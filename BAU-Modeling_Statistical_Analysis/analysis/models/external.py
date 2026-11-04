#!/usr/bin/env python3
"""Fit UKB routine-blood models and transport them to the external cohorts."""

from __future__ import annotations

import argparse
import json
import random
import shutil
from pathlib import Path
from typing import Any

import joblib
import numpy as np
import pandas as pd

from bau.calibration import apply_calibration_bundle, fit_calibration_bundle
from bau.modeling import (
    fit_and_predict,
    make_inner_folds,
    set_thread_count,
    split_calibration_sample,
    tune_model,
    tuning_sample,
)
from bau.paths import CODE_ROOT as ROOT
from bau.paths import DATA_ROOT, RESULTS_ROOT

CONFIG = ROOT / "config/external.json"
UKB_DATA = DATA_ROOT / "ukb/modeling/routine_blood.feather"
EXTERNAL_DATA = DATA_ROOT / "external"
RESULTS = RESULTS_ROOT / "models/external"
CALIBRATION = {
    "mean_crossfit_folds": 5,
    "min_sex_n": 50,
    "min_variance_factor": 0.0625,
    "max_variance_factor": 16.0,
    "phenotype_lowess_frac": 0.35,
}


def transform(values: pd.Series, rule: str) -> pd.Series:
    x = pd.to_numeric(values, errors="coerce")
    functions = {
        "identity": lambda z: z,
        "x_times_10": lambda z: z * 10,
        "x_div_10": lambda z: z / 10,
        "bun_to_urea": lambda z: z * 0.3571,
        "calcium_mgdl_to_mmoll": lambda z: z * 0.2495,
        "creatinine_mgdl_to_umoll": lambda z: z * 88.4,
        "bilirubin_mgdl_to_umoll": lambda z: z * 17.104,
        "urate_mgdl_to_umoll": lambda z: z * 59.48,
        "cholesterol_mgdl_to_mmoll": lambda z: z / 38.67,
        "triglycerides_mgdl_to_mmoll": lambda z: z / 88.57,
        "glucose_mgdl_to_mmoll": lambda z: z / 18.018,
        "hba1c_ngsp_to_ifcc": lambda z: (z - 2.152) / 0.09148,
    }
    return functions[rule](x)


def load_external(name: str, settings: dict[str, Any]) -> pd.DataFrame:
    specification = settings["cohorts"][name]
    if name == "charls":
        source = pd.read_feather(EXTERNAL_DATA / "CHARLS/candidate_qc.feather")
        eligible = (
            source["age"].between(*settings["analysis_age"])
            & source["sex"].isin(["female", "male"])
            & source["blood_eligible_initial_lt20"].fillna(False)
        )
        source = source.loc[eligible].copy().reset_index(drop=True)
        missing_column = "blood_missing_n_initial"
    else:
        source = pd.read_feather(EXTERNAL_DATA / name.upper() / "analysis_wide.feather")
        eligible = source["age"].between(*settings["analysis_age"])
        source = source.loc[eligible].copy().reset_index(drop=True)
        missing_column = "blood_missing_n"

    output = pd.DataFrame(
        {
            "eid": source["participant_id"].astype(str),
            "Age": pd.to_numeric(source["age"], errors="coerce"),
            "Sex": source["sex"].astype(str).str.lower(),
            "Baseline_unhealthy": (
                source["baseline_health_status"].astype(str) != "healthy"
            ).astype("float32"),
            "AssessmentCenter": "external",
            "SplitRole": "external_test",
        }
    )
    output["Sex_male"] = (output["Sex"] == "male").astype("float32")
    output["eid_int"] = -(np.arange(len(output), dtype="int64") + 1)
    output["BloodMissingN"] = source[missing_column].astype(int)
    output["BloodCompleteCase"] = output["BloodMissingN"].eq(0).astype("int8")
    for feature in specification["features"]:
        ukb_name = feature["ukb_feature"]
        source_name = feature["source_field"]
        if ukb_name in source:
            values = source[ukb_name]
        else:
            values = transform(source[source_name], feature["transform"])
        output[ukb_name] = pd.to_numeric(values, errors="coerce").astype("float32")

    expected = int(specification["expected_predictions"])
    feature_names = [item["ukb_feature"] for item in specification["features"]]
    missing_fraction = output[feature_names].isna().mean(axis=1)
    if len(output) != expected:
        raise ValueError(f"Expected {expected} {name} participants, found {len(output)}")
    if output["eid"].duplicated().any() or not (missing_fraction < 0.20).all():
        raise ValueError(f"Invalid {name} analysis cohort")
    return output


def load_ukb(features: list[str]) -> pd.DataFrame:
    columns = [
        "eid",
        "Age",
        "Sex",
        "Sex_male",
        "AssessmentCenter",
        "Baseline_unhealthy",
        *features,
    ]
    frame = pd.read_feather(UKB_DATA, columns=columns)
    if len(frame) != 453_163 or frame["eid"].duplicated().any():
        raise ValueError("Unexpected UKB routine-blood cohort")
    frame["eid_int"] = pd.to_numeric(frame["eid"], errors="raise").astype("int64")
    frame["Sex"] = frame["Sex"].astype(str).str.lower()
    frame["AssessmentCenter"] = frame["AssessmentCenter"].astype(int).astype(str)
    frame["Sex_male"] = frame["Sex_male"].astype("float32")
    frame["Baseline_unhealthy"] = frame["Baseline_unhealthy"].astype("float32")
    return frame


def write_json(value: Any, path: Path) -> None:
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def run(args: argparse.Namespace) -> None:
    settings = json.loads(CONFIG.read_text(encoding="utf-8"))
    specification = settings["cohorts"][args.cohort]
    features = [item["ukb_feature"] for item in specification["features"]]
    if len(features) != len(set(features)):
        raise ValueError("Feature names must be unique")

    external = load_external(args.cohort, settings)
    ukb = load_ukb(features)
    healthy = ukb.loc[ukb["Baseline_unhealthy"].eq(0)].copy().reset_index(drop=True)
    seed = int(specification["seed"])
    model_fit, calibration_reference = split_calibration_sample(
        healthy, float(settings["calibration_fraction"]), seed + 5000
    )
    model_fit["SplitRole"] = "model_fit_healthy"
    calibration_reference["SplitRole"] = "calibration_reference_healthy"

    run_dir = RESULTS / args.cohort
    if run_dir.exists():
        if not args.overwrite:
            raise FileExistsError(run_dir)
        shutil.rmtree(run_dir)
    model_dir = run_dir / "model"
    model_dir.mkdir(parents=True)

    set_thread_count(args.threads)
    random.seed(seed)
    np.random.seed(seed)
    model_features = [*features, "Sex_male"]
    if args.use_frozen_parameters:
        selected_path = CONFIG.parent / "external_parameters.json"
        parameters = json.loads(selected_path.read_text(encoding="utf-8"))[args.cohort]
    else:
        sample = tuning_sample(model_fit, int(settings["tuning_rows"]), seed)
        folds = make_inner_folds(sample, int(settings["tuning_folds"]), seed)
        parameters, _ = tune_model(
            "ngboost",
            sample,
            model_features,
            folds,
            int(settings["n_trials"]),
            seed,
            args.threads,
            run_dir,
        )
    write_json(parameters, run_dir / "parameters.json")

    scoring = pd.concat([calibration_reference, external], ignore_index=True)
    raw, _ = fit_and_predict(
        "ngboost",
        model_fit,
        scoring,
        model_features,
        parameters,
        seed + 10_000,
        model_dir,
    )
    calibration_raw = raw.loc[raw["SplitRole"].eq("calibration_reference_healthy")].reset_index(
        drop=True
    )
    external_raw = raw.loc[raw["SplitRole"].eq("external_test")].reset_index(drop=True)
    bundle, _ = fit_calibration_bundle(calibration_raw, CALIBRATION, seed + 20_000)
    predictions = apply_calibration_bundle(external_raw, bundle)
    predictions = predictions.merge(
        external[["eid_int", "eid", "BloodMissingN", "BloodCompleteCase"]],
        on="eid_int",
        how="left",
        validate="one_to_one",
    )
    predictions["Z_AgeGap_Calibrated"] = predictions["BAG"]
    predictions["Z_Uncertainty_Calibrated"] = predictions["BAU"]
    predictions["BAG_UKB"] = predictions["BAG"]
    predictions["BAU_UKB"] = predictions["BAU"]
    predictions["Cohort"] = args.cohort.upper()
    joblib.dump(bundle, model_dir / "calibration_bundle.joblib")
    predictions.to_feather(run_dir / "predictions.feather")
    predictions.to_csv(run_dir / "predictions.csv", index=False)
    write_json(
        {
            "cohort": args.cohort,
            "features": features,
            "n_ukb_model_fit": len(model_fit),
            "n_ukb_calibration": len(calibration_reference),
            "n_external": len(predictions),
            "parameter_source": "frozen_fit" if args.use_frozen_parameters else "nested_tuning",
        },
        run_dir / "summary.json",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cohort", choices=("nhanes", "charls", "chns"))
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument(
        "--use-frozen-parameters",
        action="store_true",
        help="Refit with the recorded parameter values",
    )
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
