#!/usr/bin/env python3
"""Fit one held-centre age model or collect the five out-of-fold predictions."""

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
    point_metrics,
    set_thread_count,
    split_calibration_sample,
    tune_model,
    tuning_sample,
)
from bau.paths import CODE_ROOT as ROOT
from bau.paths import DATA_ROOT, RESULTS_ROOT

CONFIG_DIR = ROOT / "config"
RESULTS_DIR = RESULTS_ROOT / "models"
REQUIRED_METADATA = {"eid", "Age", "Sex", "Sex_male", "AssessmentCenter", "Baseline_unhealthy"}


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(value: Any, path: Path) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def read_features(path: Path) -> list[str]:
    features = [
        line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()
    ]
    if len(features) != len(set(features)):
        raise ValueError(f"Repeated feature name in {path}")
    return features


def load_modality(name: str, config: dict[str, Any]) -> tuple[pd.DataFrame, list[str]]:
    specification = config["modalities"][name]
    candidate_path = DATA_ROOT / Path(specification["candidate"]).relative_to("data")
    feature_path = ROOT / specification["features"]
    features = read_features(feature_path)
    frame = pd.read_feather(candidate_path)
    missing = (REQUIRED_METADATA | set(features)).difference(frame.columns)
    if missing:
        raise ValueError(f"Missing columns in {candidate_path}: {sorted(missing)}")
    if len(frame) != int(specification["expected_rows"]):
        raise ValueError(f"Unexpected row count for {name}: {len(frame)}")
    if len(features) != int(specification["expected_features"]):
        raise ValueError(f"Unexpected feature count for {name}: {len(features)}")
    forbidden = {"Age", "Baseline_unhealthy", "AssessmentCenter", "eid", "eid_int"}
    if forbidden.intersection(features):
        raise ValueError(f"Metadata in feature list: {sorted(forbidden.intersection(features))}")
    if frame["eid"].duplicated().any():
        raise ValueError(f"Repeated participant ID in {name}")
    frame["eid_int"] = pd.to_numeric(frame["eid"], errors="raise").astype("int64")
    frame["Sex"] = frame["Sex"].astype(str).str.lower()
    frame["Sex_male"] = frame["Sex_male"].astype(np.float32)
    frame["Baseline_unhealthy"] = frame["Baseline_unhealthy"].astype(np.float32)
    frame["AssessmentCenter"] = frame["AssessmentCenter"].astype(int).astype(str)
    return frame, features + ["Sex_male"]


def run_fold(args: argparse.Namespace) -> None:
    config = read_json(CONFIG_DIR / "ukb.json")
    model_config = config["models"][args.model]
    modality_config = config["modalities"][args.modality]
    if args.model not in modality_config["tuning_rows"]:
        raise ValueError(f"{args.model} was not fitted for {args.modality}")

    frame, features = load_modality(args.modality, config)
    centre_map = read_json(CONFIG_DIR / "center_folds.json")["outer_center_to_fold"]
    frame["OuterFold"] = frame["AssessmentCenter"].map(centre_map)
    if frame["OuterFold"].isna().any():
        missing_centres = sorted(frame.loc[frame["OuterFold"].isna(), "AssessmentCenter"].unique())
        raise ValueError(f"Centres without an outer fold: {missing_centres}")
    frame["OuterFold"] = frame["OuterFold"].astype(int)

    outer_test = frame.loc[frame["OuterFold"] == args.fold].copy().reset_index(drop=True)
    development = (
        frame.loc[(frame["OuterFold"] != args.fold) & (frame["Baseline_unhealthy"] == 0)]
        .copy()
        .reset_index(drop=True)
    )
    split_seed = int(modality_config.get("seed", config["seed"])) + args.fold * 100
    model_index = list(config["models"]).index(args.model)
    model_seed = split_seed + model_index * 10
    model_fit, calibration_reference = split_calibration_sample(
        development, float(config["calibration_fraction"]), split_seed + 5000
    )
    model_fit["SplitRole"] = "model_fit_healthy"
    calibration_reference["SplitRole"] = "calibration_reference_healthy"
    outer_test["SplitRole"] = "outer_test"

    run_dir = RESULTS_DIR / args.modality / args.model / f"outer_fold_{args.fold}"
    if run_dir.exists():
        if not args.overwrite:
            raise FileExistsError(f"Output already exists: {run_dir}")
        shutil.rmtree(run_dir)
    model_dir = run_dir / "models"
    model_dir.mkdir(parents=True)

    set_thread_count(args.threads)
    random.seed(model_seed)
    np.random.seed(model_seed)
    if args.use_frozen_parameters:
        selected = read_json(CONFIG_DIR / "ukb_parameters.json")
        parameters = selected[args.modality][args.model][str(args.fold)]
    else:
        sample = tuning_sample(
            model_fit, int(modality_config["tuning_rows"][args.model]), model_seed
        )
        folds = make_inner_folds(sample, int(config["inner_folds"]), model_seed)
        parameters, _ = tune_model(
            args.model,
            sample,
            features,
            folds,
            int(model_config["trials"]),
            model_seed,
            args.threads,
            run_dir,
        )
    write_json(parameters, run_dir / "parameters.json")

    scoring = pd.concat([calibration_reference, outer_test], ignore_index=True)
    raw, _ = fit_and_predict(
        args.model, model_fit, scoring, features, parameters, model_seed + 10000, model_dir
    )
    calibration_raw = raw.loc[raw["SplitRole"] == "calibration_reference_healthy"].copy()
    outer_raw = raw.loc[raw["SplitRole"] == "outer_test"].copy()
    bundle, calibration_predictions = fit_calibration_bundle(
        calibration_raw.reset_index(drop=True), config["calibration"], model_seed + 20000
    )
    predictions = apply_calibration_bundle(outer_raw.reset_index(drop=True), bundle)
    for output in (calibration_predictions, predictions):
        output["OuterFold"] = args.fold
        output["Modality"] = args.modality
        output["Model"] = args.model
        output["ModelRole"] = model_config["role"]
    joblib.dump(bundle, model_dir / "calibration_bundle.joblib")
    calibration_predictions.to_feather(run_dir / "calibration_predictions.feather")
    predictions.to_feather(run_dir / "outer_test_predictions.feather")
    pd.concat(
        [
            model_fit[["eid_int"]].assign(role="model_fit"),
            calibration_reference[["eid_int"]].assign(role="calibration_reference"),
        ]
    ).to_csv(run_dir / "development_roles.csv", index=False)
    summary = {
        "modality": args.modality,
        "model": args.model,
        "fold": args.fold,
        "model_role": model_config["role"],
        "parameter_source": "frozen_fit" if args.use_frozen_parameters else "nested_tuning",
        "n_model_fit": len(model_fit),
        "n_calibration": len(calibration_reference),
        "n_outer_test": len(outer_test),
        "outer_centres": sorted(outer_test["AssessmentCenter"].unique()),
        "metrics": point_metrics(predictions),
    }
    write_json(summary, run_dir / "summary.json")


def aggregate(args: argparse.Namespace) -> None:
    config = read_json(CONFIG_DIR / "ukb.json")
    expected = int(config["modalities"][args.modality]["expected_rows"])
    model_root = RESULTS_DIR / args.modality / args.model
    parts = [
        pd.read_feather(model_root / f"outer_fold_{fold}/outer_test_predictions.feather")
        for fold in range(1, int(config["outer_folds"]) + 1)
    ]
    predictions = pd.concat(parts, ignore_index=True).sort_values("eid_int").reset_index(drop=True)
    if len(predictions) != expected or predictions["eid_int"].duplicated().any():
        raise ValueError("Expected one out-of-fold prediction per participant")
    predictions["Z_AgeGap_Calibrated"] = predictions["BAG"]
    predictions["Z_Uncertainty_Calibrated"] = predictions["BAU"]
    predictions.to_feather(model_root / "all_participant_predictions.feather")
    predictions.to_csv(model_root / "all_participant_predictions.csv", index=False)
    write_json(
        {"n": len(predictions), "outer_folds": sorted(predictions["OuterFold"].unique().tolist())},
        model_root / "summary.json",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    command = parser.add_subparsers(dest="command", required=True)
    fold = command.add_parser("run-fold")
    fold.add_argument("--modality", choices=("nmr", "protein", "blood"), required=True)
    fold.add_argument("--model", choices=("ngboost", "lightgbmlss", "xgboostlss"), required=True)
    fold.add_argument("--fold", type=int, choices=range(1, 6), required=True)
    fold.add_argument("--threads", type=int, default=8)
    fold.add_argument("--overwrite", action="store_true")
    fold.add_argument(
        "--use-frozen-parameters",
        action="store_true",
        help="Refit with the recorded parameter values",
    )
    fold.set_defaults(function=run_fold)
    collect = command.add_parser("aggregate")
    collect.add_argument("--modality", choices=("nmr", "protein", "blood"), required=True)
    collect.add_argument("--model", choices=("ngboost", "lightgbmlss", "xgboostlss"), required=True)
    collect.set_defaults(function=aggregate)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.function(args)


if __name__ == "__main__":
    main()
