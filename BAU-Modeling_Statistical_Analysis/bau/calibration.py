"""Calibration of the predicted age distribution and the BAG and BAU scales."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from sklearn.model_selection import StratifiedKFold
from statsmodels.nonparametric.smoothers_lowess import lowess

EPS = 1e-12


def _stratification_labels(age: np.ndarray, sex: np.ndarray, n_folds: int) -> np.ndarray:
    frame = pd.DataFrame({"age": age, "sex": sex.astype(int)})
    bins = pd.qcut(frame["age"], q=min(10, frame["age"].nunique()), labels=False, duplicates="drop")
    combined = bins.astype(str) + "_" + frame["sex"].astype(str)
    if int(combined.value_counts().min()) >= n_folds:
        return combined.to_numpy()
    return bins.astype(str).to_numpy()


def _fit_affine(mu: np.ndarray, age: np.ndarray) -> dict[str, float]:
    design = np.column_stack([np.ones(len(mu)), mu])
    coefficients, *_ = np.linalg.lstsq(design, age, rcond=None)
    intercept, slope = map(float, coefficients)
    if not np.isfinite([intercept, slope]).all() or slope <= 0:
        raise RuntimeError(
            f"Invalid monotone mean calibration: intercept={intercept}, slope={slope}"
        )
    return {"intercept": intercept, "slope": slope, "n": int(len(mu))}


def fit_mean_calibrator(
    mu: np.ndarray,
    age: np.ndarray,
    sex: np.ndarray,
    min_sex_n: int,
) -> dict[str, Any]:
    finite = np.isfinite(mu) & np.isfinite(age) & np.isfinite(sex)
    if int(finite.sum()) < max(20, min_sex_n):
        raise RuntimeError("Insufficient finite rows for mean calibration")
    model: dict[str, Any] = {
        "kind": "sex_specific_affine_age_on_predicted_mean",
        "application_inputs": ["BA_Raw", "Sex_male"],
        "pooled": _fit_affine(mu[finite], age[finite]),
        "sex_models": {},
    }
    for sex_value in (0, 1):
        mask = finite & (sex.astype(int) == sex_value)
        if int(mask.sum()) >= min_sex_n:
            model["sex_models"][str(sex_value)] = _fit_affine(mu[mask], age[mask])
    return model


def predict_mean_calibrator(model: dict[str, Any], mu: np.ndarray, sex: np.ndarray) -> np.ndarray:
    mu = np.asarray(mu, dtype=float)
    sex = np.asarray(sex, dtype=float)
    pooled = model["pooled"]
    output = float(pooled["intercept"]) + float(pooled["slope"]) * mu
    for sex_value in (0, 1):
        spec = model["sex_models"].get(str(sex_value))
        if spec is None:
            continue
        mask = sex.astype(int) == sex_value
        output[mask] = float(spec["intercept"]) + float(spec["slope"]) * mu[mask]
    return output


def crossfit_mean_predictions(
    raw: pd.DataFrame,
    n_folds: int,
    seed: int,
    min_sex_n: int,
) -> tuple[np.ndarray, list[dict[str, Any]]]:
    age = raw["Age"].to_numpy(float)
    sex = raw["Sex_male"].to_numpy(float)
    mu = raw["BA_Raw"].to_numpy(float)
    labels = _stratification_labels(age, sex, n_folds)
    splitter = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=seed)
    prediction = np.full(len(raw), np.nan, dtype=float)
    records: list[dict[str, Any]] = []
    for fold_id, (train_idx, valid_idx) in enumerate(
        splitter.split(np.zeros(len(raw)), labels), start=1
    ):
        model = fit_mean_calibrator(mu[train_idx], age[train_idx], sex[train_idx], min_sex_n)
        prediction[valid_idx] = predict_mean_calibrator(model, mu[valid_idx], sex[valid_idx])
        records.append(
            {"fold": fold_id, "n_train": int(len(train_idx)), "n_valid": int(len(valid_idx))}
        )
    if not np.isfinite(prediction).all():
        raise RuntimeError("Cross-fitted mean calibration predictions are incomplete")
    return prediction, records


def fit_variance_calibrator(
    residual: np.ndarray,
    variance: np.ndarray,
    sex: np.ndarray,
    min_sex_n: int,
    min_variance_factor: float,
    max_variance_factor: float,
) -> dict[str, Any]:
    residual = np.asarray(residual, dtype=float)
    variance = np.maximum(np.asarray(variance, dtype=float), EPS)
    sex = np.asarray(sex, dtype=float)
    finite = np.isfinite(residual) & np.isfinite(variance) & np.isfinite(sex)

    def factor(mask: np.ndarray) -> float:
        value = float(np.mean((residual[mask] ** 2) / variance[mask]))
        return float(np.clip(value, min_variance_factor, max_variance_factor))

    model: dict[str, Any] = {
        "kind": "sex_specific_global_variance_multiplier",
        "application_inputs": ["Var_Raw", "Sex_male"],
        "pooled_factor": factor(finite),
        "sex_factors": {},
        "min_variance_factor": float(min_variance_factor),
        "max_variance_factor": float(max_variance_factor),
    }
    for sex_value in (0, 1):
        mask = finite & (sex.astype(int) == sex_value)
        if int(mask.sum()) >= min_sex_n:
            model["sex_factors"][str(sex_value)] = factor(mask)
    return model


def predict_variance_calibrator(
    model: dict[str, Any], variance: np.ndarray, sex: np.ndarray
) -> np.ndarray:
    variance = np.maximum(np.asarray(variance, dtype=float), EPS)
    sex = np.asarray(sex, dtype=float)
    factors = np.full(len(variance), float(model["pooled_factor"]), dtype=float)
    for sex_value in (0, 1):
        mask = sex.astype(int) == sex_value
        factors[mask] = float(model["sex_factors"].get(str(sex_value), model["pooled_factor"]))
    return variance * factors


def apply_predictive_calibration(raw: pd.DataFrame, calibrator: dict[str, Any]) -> pd.DataFrame:
    """Apply the fitted mean and variance corrections."""
    out = raw.copy()
    mu = out["BA_Raw"].to_numpy(float)
    sex = out["Sex_male"].to_numpy(float)
    variance = out["Var_Raw"].to_numpy(float)
    mu_cal = predict_mean_calibrator(calibrator["mean"], mu, sex)
    var_cal = predict_variance_calibrator(calibrator["variance"], variance, sex)
    out["BA_Calibrated"] = mu_cal
    out["Var_Calibrated"] = np.maximum(var_cal, EPS)
    out["Sigma_Calibrated"] = np.sqrt(out["Var_Calibrated"].to_numpy(float))
    out["LogVar_Calibrated"] = np.log(out["Var_Calibrated"].to_numpy(float))
    return out


def _fit_sex_lowess(
    age: np.ndarray,
    sex: np.ndarray,
    values: np.ndarray,
    frac: float,
    min_sex_n: int,
) -> dict[str, Any]:
    finite = np.isfinite(age) & np.isfinite(sex) & np.isfinite(values)
    model: dict[str, Any] = {
        "fallback": float(np.mean(values[finite])),
        "frac": float(frac),
        "sex_models": {},
    }
    for sex_value in (0, 1):
        mask = finite & (sex.astype(int) == sex_value)
        if int(mask.sum()) < min_sex_n:
            continue
        fitted = lowess(values[mask], age[mask], frac=frac, return_sorted=True)
        curve = (
            pd.DataFrame({"age": fitted[:, 0], "value": fitted[:, 1]})
            .groupby("age", as_index=False)["value"]
            .mean()
        )
        model["sex_models"][str(sex_value)] = {
            "age": curve["age"].to_numpy(float),
            "value": curve["value"].to_numpy(float),
        }
    return model


def _predict_sex_lowess(model: dict[str, Any], age: np.ndarray, sex: np.ndarray) -> np.ndarray:
    prediction = np.full(len(age), float(model["fallback"]), dtype=float)
    for sex_value in (0, 1):
        spec = model["sex_models"].get(str(sex_value))
        mask = sex.astype(int) == sex_value
        if spec is None or not mask.any():
            continue
        function = interp1d(
            spec["age"],
            spec["value"],
            bounds_error=False,
            fill_value=(spec["value"][0], spec["value"][-1]),
        )
        prediction[mask] = function(age[mask])
    return prediction


def fit_phenotype_standardizer(
    calibrated_reference: pd.DataFrame,
    lowess_frac: float,
    min_sex_n: int,
) -> dict[str, Any]:
    age = calibrated_reference["Age"].to_numpy(float)
    sex = calibrated_reference["Sex_male"].to_numpy(float)
    result: dict[str, Any] = {
        "kind": "sex_specific_age_lowess_healthy_reference_z",
        "metrics": {},
    }
    values_by_name = {
        "bag": calibrated_reference["BA_Calibrated"].to_numpy(float) - age,
        "bau": calibrated_reference["LogVar_Calibrated"].to_numpy(float),
    }
    for name, values in values_by_name.items():
        baseline = _fit_sex_lowess(age, sex, values, lowess_frac, min_sex_n)
        residual = values - _predict_sex_lowess(baseline, age, sex)
        pooled_sd = float(np.std(residual, ddof=0)) or 1.0
        sex_mean: dict[str, float] = {}
        sex_sd: dict[str, float] = {}
        for sex_value in (0, 1):
            mask = sex.astype(int) == sex_value
            sex_mean[str(sex_value)] = float(np.mean(residual[mask])) if mask.any() else 0.0
            sex_sd[str(sex_value)] = (
                float(np.std(residual[mask], ddof=0)) if mask.any() else pooled_sd
            )
            if sex_sd[str(sex_value)] <= 0:
                sex_sd[str(sex_value)] = pooled_sd
        result["metrics"][name] = {
            "baseline": baseline,
            "sex_mean": sex_mean,
            "sex_sd": sex_sd,
            "pooled_sd": pooled_sd,
        }
    return result


def apply_phenotype_standardizer(df: pd.DataFrame, standardizer: dict[str, Any]) -> pd.DataFrame:
    out = df.copy()
    age = out["Age"].to_numpy(float)
    sex = out["Sex_male"].to_numpy(float)
    values_by_name = {
        "bag": out["BA_Calibrated"].to_numpy(float) - age,
        "bau": out["LogVar_Calibrated"].to_numpy(float),
    }
    output_names = {"bag": "BAG", "bau": "BAU"}
    for name, values in values_by_name.items():
        spec = standardizer["metrics"][name]
        baseline = _predict_sex_lowess(spec["baseline"], age, sex)
        center = np.array([float(spec["sex_mean"].get(str(int(x)), 0.0)) for x in sex])
        scale = np.array([float(spec["sex_sd"].get(str(int(x)), spec["pooled_sd"])) for x in sex])
        residual = values - baseline - center
        prefix = name.upper()
        out[f"{prefix}_AgeSexBaseline"] = baseline
        out[f"{prefix}_Residual"] = residual
        out[output_names[name]] = residual / scale
    out["AgeGap_Calibrated"] = out["BA_Calibrated"].to_numpy(float) - age
    out["ZError_Calibrated"] = (age - out["BA_Calibrated"].to_numpy(float)) / out[
        "Sigma_Calibrated"
    ].to_numpy(float)
    return out


def fit_calibration_bundle(
    raw_reference: pd.DataFrame,
    config: dict[str, Any],
    seed: int,
) -> tuple[dict[str, Any], pd.DataFrame]:
    age = raw_reference["Age"].to_numpy(float)
    sex = raw_reference["Sex_male"].to_numpy(float)
    mu = raw_reference["BA_Raw"].to_numpy(float)
    variance = raw_reference["Var_Raw"].to_numpy(float)
    crossfit_mu, crossfit_records = crossfit_mean_predictions(
        raw_reference,
        int(config["mean_crossfit_folds"]),
        seed,
        int(config["min_sex_n"]),
    )
    # Estimate variance scaling from cross-fitted mean residuals.
    crossfit_residual = age - crossfit_mu
    mean_model = fit_mean_calibrator(mu, age, sex, int(config["min_sex_n"]))
    variance_model = fit_variance_calibrator(
        crossfit_residual,
        variance,
        sex,
        int(config["min_sex_n"]),
        float(config["min_variance_factor"]),
        float(config["max_variance_factor"]),
    )
    predictive = {
        "kind": "affine_mean_and_global_variance_scale",
        "application_inputs": ["BA_Raw", "Var_Raw", "Sex_male"],
        "mean": mean_model,
        "variance": variance_model,
    }
    calibrated = apply_predictive_calibration(raw_reference, predictive)
    standardizer = fit_phenotype_standardizer(
        calibrated,
        float(config["phenotype_lowess_frac"]),
        int(config["min_sex_n"]),
    )
    calibrated = apply_phenotype_standardizer(calibrated, standardizer)
    bundle = {
        "schema_version": 3,
        "fit_population": "independent_outer_development_healthy_calibration_reference",
        "predictive_calibrator": predictive,
        "phenotype_standardizer": standardizer,
        "crossfit_mean_records": crossfit_records,
        "crossfit_residual_rmse": float(np.sqrt(np.mean(crossfit_residual**2))),
        "config": dict(config),
    }
    return bundle, calibrated


def apply_calibration_bundle(raw: pd.DataFrame, bundle: dict[str, Any]) -> pd.DataFrame:
    calibrated = apply_predictive_calibration(raw, bundle["predictive_calibrator"])
    return apply_phenotype_standardizer(calibrated, bundle["phenotype_standardizer"])
