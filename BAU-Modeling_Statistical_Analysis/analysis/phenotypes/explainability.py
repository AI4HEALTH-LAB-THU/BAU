"""Centre-held-out factor models and total explained variance for Figure S2e."""

from __future__ import annotations

import json
import math
import os
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from lightgbm import LGBMRegressor
from scipy import stats
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder

from bau.paths import CODE_ROOT as ROOT
from bau.paths import DATA_ROOT as PREPARED_ROOT
from bau.paths import RESULTS_ROOT

UKB_ROOT = Path(os.environ.get("UKB_DATA_ROOT", "/mnt/project/UKB"))
DATA_ROOT = PREPARED_ROOT / "ukb"
OUTPUT = RESULTS_ROOT / "phenotypes/explainability"

SEED = 20260827
N_BOOT = 2000

COVARIATE_PATH = DATA_ROOT / "support/covariates.feather"
HISTORY_PATH = DATA_ROOT / "support/disease_history.feather"
ANTHRO_PATH = UKB_ROOT / "ProcessedData/02-AssessmentCenter/05-PhysicalMeasures/anthropometry_0.csv"
BP_PATH = UKB_ROOT / "ProcessedData/02-AssessmentCenter/05-PhysicalMeasures/blood_pressure_0.csv"
DRUG_RDS = UKB_ROOT / "RawData/02-AssessmentCenter/02-Touchscreen/drug_0.rds"

BASE_NUM = ["age", "sex_male", "bag"]
TECH_NUM = ["assay_missing_fraction", "attendance_days", "season_sin", "season_cos"]
LIFE_NUM = [
    "current_smoker",
    "daily_alcohol",
    "daily_healthy_food",
    "frequent_processed_meat",
    "frequent_red_meat",
    "frequent_salt",
    "mvpa_log1p",
]
LIFE_CAT = ["income", "education", "employment"]
MED_NUM = ["cholesterol_lowering", "blood_pressure_medication", "insulin"]
CLINICAL_NUM = [
    "hypertension_history",
    "diabetes_history",
    "bmi",
    "waist_hip_ratio",
    "systolic_bp",
    "diastolic_bp",
]

BLOCKS = {
    "Technical": (TECH_NUM, []),
    "Lifestyle and socioeconomic": (LIFE_NUM, LIFE_CAT),
    "Medication": (MED_NUM, []),
    "Baseline clinical state": (CLINICAL_NUM, []),
    "All measured factors": (
        TECH_NUM + LIFE_NUM + MED_NUM + CLINICAL_NUM,
        LIFE_CAT,
    ),
}


def build_drug_table() -> pd.DataFrame:
    cache = OUTPUT / "medication_fields.feather"
    if not cache.exists():
        expression = (
            "x <- as.data.frame(readRDS(" + json.dumps(str(DRUG_RDS)) + ")); "
            "arrow::write_feather(x, " + json.dumps(str(cache)) + ", compression='zstd')"
        )
        subprocess.run(["Rscript", "-e", expression], check=True)
    raw = pd.read_feather(cache)
    value_columns = [column for column in raw.columns if column != "eid"]
    values = raw[value_columns].to_numpy(dtype=float)
    recognised = np.isin(values, [-7, 1, 2, 3, 4, 5]).any(axis=1)
    out = pd.DataFrame({"eid": raw["eid"].astype("int64")})
    for code, name in (
        (1, "cholesterol_lowering"),
        (2, "blood_pressure_medication"),
        (3, "insulin"),
    ):
        indicator = np.full(len(raw), np.nan, dtype=float)
        indicator[recognised] = (values[recognised] == code).any(axis=1).astype(float)
        out[name] = indicator
    if out["eid"].duplicated().any():
        raise ValueError("Duplicated eid in baseline medication table")
    return out


def load_supporting_tables() -> dict[str, pd.DataFrame]:
    covariate_columns = [
        "eid",
        "Income",
        "Employment",
        "Education",
        "Current smoker",
        "Daily alcohol intake",
        "Daily healthy food",
        "Frequently processed meat",
        "Frequently red meat",
        "Frequently salt intake",
        "Daily moderate to vigorous physical activity",
    ]
    covariates = pd.read_feather(COVARIATE_PATH, columns=covariate_columns).rename(
        columns={
            "Income": "income",
            "Employment": "employment",
            "Education": "education",
            "Current smoker": "current_smoker",
            "Daily alcohol intake": "daily_alcohol",
            "Daily healthy food": "daily_healthy_food",
            "Frequently processed meat": "frequent_processed_meat",
            "Frequently red meat": "frequent_red_meat",
            "Frequently salt intake": "frequent_salt",
            "Daily moderate to vigorous physical activity": "mvpa",
        }
    )
    covariates["eid"] = covariates["eid"].astype("int64")
    for column in LIFE_NUM[:-1]:
        covariates[column] = pd.to_numeric(
            covariates[column].astype("string"), errors="coerce"
        ).astype(float)
    covariates["mvpa"] = pd.to_numeric(covariates["mvpa"], errors="coerce")
    covariates["mvpa_log1p"] = np.log1p(covariates["mvpa"].clip(lower=0))
    covariates = covariates.drop(columns="mvpa")
    for column in LIFE_CAT:
        covariates[column] = covariates[column].astype("string")

    histories = pd.read_feather(
        HISTORY_PATH,
        columns=["eid", "hypertension_history", "diabetes_mellitus_history"],
    ).rename(columns={"diabetes_mellitus_history": "diabetes_history"})
    histories["eid"] = histories["eid"].astype("int64")

    anthropometry = pd.read_csv(
        ANTHRO_PATH,
        usecols=["eid", "Body mass index (BMI)", "Waist circumference", "Hip circumference"],
    ).rename(columns={"Body mass index (BMI)": "bmi"})
    anthropometry["eid"] = anthropometry["eid"].astype("int64")
    anthropometry["waist_hip_ratio"] = (
        anthropometry["Waist circumference"] / anthropometry["Hip circumference"]
    )
    anthropometry.loc[~anthropometry["bmi"].between(10, 80), "bmi"] = np.nan
    anthropometry.loc[~anthropometry["waist_hip_ratio"].between(0.3, 2.0), "waist_hip_ratio"] = (
        np.nan
    )
    anthropometry = anthropometry[["eid", "bmi", "waist_hip_ratio"]]

    blood_pressure = pd.read_csv(
        BP_PATH,
        usecols=["eid", "Systolic blood pressure", "Diastolic blood pressure"],
    ).rename(
        columns={
            "Systolic blood pressure": "systolic_bp",
            "Diastolic blood pressure": "diastolic_bp",
        }
    )
    blood_pressure["eid"] = blood_pressure["eid"].astype("int64")
    blood_pressure.loc[~blood_pressure["systolic_bp"].between(60, 300), "systolic_bp"] = np.nan
    blood_pressure.loc[~blood_pressure["diastolic_bp"].between(30, 200), "diastolic_bp"] = np.nan

    return {
        "covariates": covariates,
        "histories": histories,
        "anthropometry": anthropometry,
        "blood_pressure": blood_pressure,
        "medication": build_drug_table(),
    }


def prediction_paths() -> dict[str, Path]:
    return {
        modality: RESULTS_ROOT / f"models/{modality}/ngboost/all_participant_predictions.csv"
        for modality in ("nmr", "protein")
    }


def load_modality(
    modality: str, prediction_path: Path, support: dict[str, pd.DataFrame]
) -> pd.DataFrame:
    prediction = pd.read_csv(
        prediction_path,
        usecols=[
            "eid_int",
            "Age",
            "Sex_male",
            "AssessmentCenter",
            "OuterFold",
            "Z_AgeGap_Calibrated",
            "Z_Uncertainty_Calibrated",
        ],
    ).rename(
        columns={
            "eid_int": "eid",
            "Age": "age",
            "Sex_male": "sex_male",
            "AssessmentCenter": "centre",
            "OuterFold": "outer_fold",
            "Z_AgeGap_Calibrated": "bag",
            "Z_Uncertainty_Calibrated": "bau",
        }
    )
    prediction["eid"] = prediction["eid"].astype("int64")
    if prediction.isna().any().any():
        raise ValueError(f"{modality}: invalid prediction table")

    candidate_path = DATA_ROOT / f"modeling/{modality}.feather"
    feature_order = (ROOT / f"config/features/{modality}.txt").read_text().splitlines()
    feature_order = [feature for feature in feature_order if feature]
    candidate = pd.read_feather(candidate_path, columns=["eid", "Date_attend", *feature_order])
    candidate["eid"] = candidate["eid"].astype("int64")
    candidate["assay_missing_fraction"] = candidate[feature_order].isna().mean(axis=1)
    candidate = candidate[["eid", "Date_attend", "assay_missing_fraction"]]

    frame = prediction.merge(candidate, on="eid", how="left", validate="one_to_one")
    if frame[["Date_attend", "assay_missing_fraction"]].isna().any().any():
        raise ValueError(f"{modality}: incomplete candidate linkage")
    attendance = pd.to_datetime(frame.pop("Date_attend"), errors="raise")
    frame["attendance_days"] = (attendance - pd.Timestamp("2006-01-01")).dt.days.astype(float)
    day_of_year = attendance.dt.dayofyear.astype(float)
    frame["season_sin"] = np.sin(2 * np.pi * day_of_year / 365.25)
    frame["season_cos"] = np.cos(2 * np.pi * day_of_year / 365.25)

    for table in support.values():
        frame = frame.merge(table, on="eid", how="left", validate="one_to_one")
    if (
        frame[["income", "education", "employment", "hypertension_history", "diabetes_history"]]
        .isna()
        .any()
        .any()
    ):
        raise ValueError(f"{modality}: unexpected missing supporting categories/history")
    return frame


def make_pipeline(numeric: list[str], categorical: list[str], seed: int) -> Pipeline:
    transformers = [("numeric", "passthrough", numeric)]
    if categorical:
        transformers.append(
            (
                "categorical",
                OneHotEncoder(handle_unknown="ignore"),
                categorical,
            )
        )
    preprocess = ColumnTransformer(transformers, remainder="drop", sparse_threshold=1.0)
    model = LGBMRegressor(
        objective="regression",
        n_estimators=300,
        learning_rate=0.03,
        num_leaves=15,
        min_child_samples=200,
        reg_lambda=1.0,
        reg_alpha=0.0,
        subsample=1.0,
        colsample_bytree=1.0,
        random_state=seed,
        n_jobs=int(os.environ.get("BAU_THREADS", "8")),
        verbosity=-1,
        deterministic=True,
        force_col_wise=True,
    )
    return Pipeline([("preprocess", preprocess), ("lightgbm", model)])


def cross_fit(
    frame: pd.DataFrame, numeric: list[str], categorical: list[str], seed: int
) -> np.ndarray:
    predictions = np.full(len(frame), np.nan, dtype=float)
    for fold in sorted(frame["outer_fold"].unique()):
        train = frame["outer_fold"].to_numpy() != fold
        test = ~train
        pipe = make_pipeline(numeric, categorical, seed + int(fold))
        pipe.fit(frame.loc[train, numeric + categorical], frame.loc[train, "bau"])
        predictions[test] = pipe.predict(frame.loc[test, numeric + categorical])
    if not np.isfinite(predictions).all():
        raise ValueError("Cross-fitting did not produce finite predictions")
    return predictions


def r2(y: np.ndarray, prediction: np.ndarray) -> float:
    sst = np.square(y - y.mean()).sum()
    return 1.0 - np.square(y - prediction).sum() / sst


def clustered_r2_ci(frame, prediction, rng):
    y = frame["bau"].to_numpy()
    work = pd.DataFrame(
        {"centre": frame["centre"], "y": y, "y2": y * y, "sse": np.square(y - prediction)}
    )
    values = (
        work.groupby("centre", sort=True)
        .agg(n=("y", "size"), sum_y=("y", "sum"), sum_y2=("y2", "sum"), sse=("sse", "sum"))
        .to_numpy(float)
    )
    sampled = values[rng.integers(0, len(values), size=(N_BOOT, len(values)))].sum(axis=1)
    sst = sampled[:, 2] - sampled[:, 1] ** 2 / sampled[:, 0]
    draws = 1 - sampled[:, 3] / sst
    return tuple(np.quantile(draws, [0.025, 0.975]))


def leave_one_centre_ci(
    frame: pd.DataFrame, residual: np.ndarray, metric: str
) -> tuple[float, float, float]:
    raw = frame["bau"].to_numpy()
    centres = frame["centre"].to_numpy()

    def estimate(mask: np.ndarray) -> float:
        return float(stats.spearmanr(raw[mask], residual[mask]).statistic)

    point = estimate(np.ones(len(frame), dtype=bool))
    unique = np.unique(centres)
    jack = np.array([estimate(centres != centre) for centre in unique])
    mean_jack = jack.mean()
    se = math.sqrt((len(unique) - 1) / len(unique) * np.square(jack - mean_jack).sum())
    critical = stats.t.ppf(0.975, df=len(unique) - 1)
    lower, upper = point - critical * se, point + critical * se
    lower, upper = max(-1.0, lower), min(1.0, upper)
    return point, lower, upper


def analyse_modality(
    modality: str, frame: pd.DataFrame
) -> tuple[list[dict], list[dict], list[dict]]:
    y = frame["bau"].to_numpy()
    modality_index = 0 if modality == "nmr" else 1
    seed = SEED + 100 * modality_index
    base_prediction = cross_fit(frame, BASE_NUM, [], seed)
    predictions = {"Age, sex and BAG": base_prediction}
    for block, (numeric, categorical) in BLOCKS.items():
        predictions[block] = cross_fit(frame, BASE_NUM + numeric, categorical, seed)

    base_r2 = r2(y, base_prediction)
    main_rows = []
    fold_rows = []
    for order, block in enumerate(predictions):
        prediction = predictions[block]
        block_r2 = r2(y, prediction)
        rng = np.random.default_rng(SEED + 1000 * modality_index + (order - 1 if order else 99))
        ci_low, ci_high = clustered_r2_ci(frame, prediction, rng)
        main_rows.append(
            {
                "modality": modality,
                "block": block,
                "block_order": order,
                "n": len(frame),
                "centres": frame["centre"].nunique(),
                "reference_oof_r2": base_r2,
                "oof_r2": block_r2,
                "ci_low": ci_low,
                "ci_high": ci_high,
                "bootstrap_replicates": N_BOOT,
            }
        )
        for fold in sorted(frame["outer_fold"].unique()):
            test = frame["outer_fold"].to_numpy() == fold
            fold_base = r2(y[test], base_prediction[test])
            fold_block = r2(y[test], prediction[test])
            fold_rows.append(
                {
                    "modality": modality,
                    "block": block,
                    "outer_fold": int(fold),
                    "n": int(test.sum()),
                    "test_centres": int(frame.loc[test, "centre"].nunique()),
                    "reference_oof_r2": fold_base,
                    "oof_r2": fold_block,
                }
            )

    combined_prediction = predictions["All measured factors"]
    residual = y - combined_prediction
    identity_rows = []
    for metric in ("spearman_r",):
        estimate, ci_low, ci_high = leave_one_centre_ci(frame, residual, metric)
        identity_rows.append(
            {
                "modality": modality,
                "metric": metric,
                "estimate": estimate,
                "ci_low": ci_low,
                "ci_high": ci_high,
                "ci_method": "leave-one-centre-out jackknife",
            }
        )

    cache = frame[["eid", "centre", "outer_fold", "bau"]].copy()
    for order, (block, prediction) in enumerate(predictions.items()):
        cache[f"prediction_{order}"] = prediction
    cache["residual_bau"] = residual
    cache.to_feather(OUTPUT / f"{modality}_oof_predictions.feather")
    return main_rows, identity_rows, fold_rows


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    support = load_supporting_tables()
    paths = prediction_paths()
    all_main, all_identity, all_fold = [], [], []
    cohort_rows = []
    missing_rows = []
    for modality in ("nmr", "protein"):
        print(f"[{modality}] loading and linking", flush=True)
        frame = load_modality(modality, paths[modality], support)
        cohort_rows.append(
            {
                "modality": modality,
                "participants": len(frame),
                "assessment_centres": frame["centre"].nunique(),
                "outer_folds": frame["outer_fold"].nunique(),
            }
        )
        for column in TECH_NUM + LIFE_NUM + LIFE_CAT + MED_NUM + CLINICAL_NUM:
            missing_rows.append(
                {
                    "modality": modality,
                    "variable": column,
                    "missing_n": int(frame[column].isna().sum()),
                    "missing_fraction": float(frame[column].isna().mean()),
                }
            )
        print(f"[{modality}] cross-fitting {len(frame):,} participants", flush=True)
        main_rows, identity_rows, fold_rows = analyse_modality(modality, frame)
        all_main.extend(main_rows)
        all_identity.extend(identity_rows)
        all_fold.extend(fold_rows)

    pd.DataFrame(all_main).to_csv(OUTPUT / "total_oof_r2.tsv", sep="\t", index=False)
    pd.DataFrame(all_identity).to_csv(OUTPUT / "residual_identity.tsv", sep="\t", index=False)
    pd.DataFrame(all_fold).to_csv(OUTPUT / "fold_metrics.tsv", sep="\t", index=False)
    pd.DataFrame(cohort_rows).to_csv(OUTPUT / "cohort_summary.tsv", sep="\t", index=False)
    pd.DataFrame(missing_rows).to_csv(OUTPUT / "covariate_missingness.tsv", sep="\t", index=False)

    feature_rows = []
    for block, (numeric, categorical) in BLOCKS.items():
        for feature in numeric:
            feature_rows.append({"block": block, "variable": feature, "type": "numeric"})
        for feature in categorical:
            feature_rows.append({"block": block, "variable": feature, "type": "categorical"})
    pd.DataFrame(feature_rows).to_csv(OUTPUT / "feature_blocks.tsv", sep="\t", index=False)

    main_table = pd.DataFrame(all_main)
    identity_table = pd.DataFrame(all_identity)
    qc = {
        "seed": SEED,
        "bootstrap_replicates": N_BOOT,
        "auxiliary_model": "LightGBM",
        "modalities": cohort_rows,
        "all_oof_predictions_finite": True,
        "total_oof_r2": main_table.to_dict(orient="records"),
        "residual_identity": identity_table.to_dict(orient="records"),
    }
    (OUTPUT / "summary.json").write_text(json.dumps(qc, indent=2) + "\n")
    print(main_table.to_string(index=False), flush=True)
    print(identity_table.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
