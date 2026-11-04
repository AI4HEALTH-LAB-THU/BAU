#!/usr/bin/env python3
"""Score the repeat NMR assessment with each participant's baseline fold model."""

from __future__ import annotations

import json
import os
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

from bau.calibration import apply_calibration_bundle
from bau.modeling import predict_distribution
from bau.paths import CODE_ROOT as ROOT
from bau.paths import DATA_ROOT as PREPARED_ROOT
from bau.paths import RESULTS_ROOT

UKB_ROOT = Path(os.environ.get("UKB_DATA_ROOT", "/mnt/project/UKB"))


MODEL_ROOT = RESULTS_ROOT / "models/nmr/ngboost"
PREDICTIONS = MODEL_ROOT / "all_participant_predictions.csv"
CANDIDATE = PREPARED_ROOT / "ukb/modeling/nmr.feather"
FEATURES = ROOT / "config/features/nmr.txt"
REPEAT_NMR = UKB_ROOT / "ProcessedData/03-BiologicalSamples/NMR_instance_1.feather"
REPEAT_VISIT = UKB_ROOT / "RawData/02-AssessmentCenter/01-Recruitment/recruitment_1.csv"
OUTPUT = RESULTS_ROOT / "phenotypes/repeat_pairs.feather"


def score(frame: pd.DataFrame, fold: int) -> pd.DataFrame:
    model_dir = MODEL_ROOT / f"outer_fold_{fold}/models"
    saved = joblib.load(model_dir / "final_model.joblib")
    calibration = joblib.load(model_dir / "calibration_bundle.joblib")
    matrix = frame[saved["features"]].to_numpy(dtype="float32")
    matrix = saved["preprocessor"].transform(matrix)
    predicted = predict_distribution("ngboost", saved["model"], matrix)
    raw = pd.DataFrame(predicted)
    raw["Age"] = frame["repeat_age"].to_numpy(float)
    raw["Sex_male"] = frame["Sex_male"].to_numpy(float)
    return apply_calibration_bundle(raw, calibration)


def main() -> None:
    features = [line.strip() for line in FEATURES.read_text().splitlines() if line.strip()]
    if len(features) != 168 or len(features) != len(set(features)):
        raise ValueError("The NMR feature list must contain 168 unique measures")
    baseline = pd.read_csv(PREDICTIONS)
    if len(baseline) != 272_672:
        raise ValueError("Unexpected baseline NMR predictions")
    metadata = pd.read_feather(
        CANDIDATE,
        columns=[
            "eid",
            "Date_birth",
            "Date_attend",
            "Sex",
            "Sex_male",
            "AssessmentCenter",
        ],
    )
    repeat = pd.read_feather(REPEAT_NMR, columns=["eid", *features])
    visit = pd.read_csv(
        REPEAT_VISIT,
        usecols=[
            "eid",
            "Date of attending assessment centre",
            "UK Biobank assessment centre",
        ],
    ).rename(
        columns={
            "Date of attending assessment centre": "Repeat_date",
            "UK Biobank assessment centre": "Repeat_assessment_center",
        }
    )
    for frame in (metadata, repeat, visit):
        frame["eid"] = pd.to_numeric(frame["eid"], errors="raise").astype("int64")

    baseline = baseline.merge(
        metadata,
        left_on="eid_int",
        right_on="eid",
        validate="one_to_one",
        suffixes=("", "_candidate"),
    )
    linked = repeat.merge(visit, on="eid", validate="one_to_one").merge(
        baseline, on="eid", validate="one_to_one"
    )
    linked["Repeat_date"] = pd.to_datetime(linked["Repeat_date"], errors="coerce")
    linked["Date_birth"] = pd.to_datetime(linked["Date_birth"], errors="coerce")
    linked["Date_attend"] = pd.to_datetime(linked["Date_attend"], errors="coerce")
    linked["repeat_age"] = ((linked["Repeat_date"] - linked["Date_birth"]).dt.days / 365.25).round(
        1
    )
    linked["visit_interval_years"] = (
        linked["Repeat_date"] - linked["Date_attend"]
    ).dt.days / 365.25
    linked["repeat_feature_missing_fraction"] = linked[features].isna().mean(axis=1)
    linked = linked.loc[
        linked["repeat_age"].between(40, 70)
        & linked["visit_interval_years"].gt(0)
        & linked["repeat_feature_missing_fraction"].lt(0.20)
    ].copy()
    linked["outer_fold"] = linked["OuterFold"].astype(int)

    parts = []
    for fold in range(1, 6):
        frame = linked.loc[linked["outer_fold"].eq(fold)].copy()
        calibrated = score(frame, fold)
        parts.append(
            pd.DataFrame(
                {
                    "eid": frame["eid"].to_numpy("int64"),
                    "baseline_age": frame["Age"].to_numpy(float),
                    "repeat_age": frame["repeat_age"].to_numpy(float),
                    "visit_interval_years": frame["visit_interval_years"].to_numpy(float),
                    "sex": frame["Sex"].astype(str).str.lower().to_numpy(),
                    "sex_male": frame["Sex_male"].to_numpy(float),
                    "baseline_center": frame["AssessmentCenter"].astype(str).to_numpy(),
                    "repeat_center": frame["Repeat_assessment_center"].to_numpy(),
                    "outer_fold": fold,
                    "repeat_feature_missing_fraction": frame[
                        "repeat_feature_missing_fraction"
                    ].to_numpy(float),
                    "baseline_bag": frame["BAG"].to_numpy(float),
                    "baseline_bau": frame["BAU"].to_numpy(float),
                    "repeat_bag": calibrated["BAG"].to_numpy(float),
                    "repeat_bau": calibrated["BAU"].to_numpy(float),
                }
            )
        )
    output = pd.concat(parts, ignore_index=True).sort_values("eid")
    numeric = output.select_dtypes(include="number")
    if len(output) != 13_620:
        raise ValueError(f"Expected 13620 repeat pairs, found {len(output)}")
    if not np.isfinite(numeric.to_numpy(float)).all():
        raise ValueError("Repeat scores contain non-finite values")
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    output.reset_index(drop=True).to_feather(OUTPUT)
    print(json.dumps({"participants": len(output), "folds": 5}, indent=2))


if __name__ == "__main__":
    main()
