from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Any

import joblib
import numpy as np
import optuna
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.metrics import mean_absolute_error, mean_squared_error
from sklearn.model_selection import StratifiedKFold, StratifiedShuffleSplit
from sklearn.preprocessing import StandardScaler

MIN_SIGMA = 1e-6
MIN_VARIANCE = 1e-12
FAILED_TRIAL_NLL = 1e12
LOGGER = logging.getLogger(__name__)


class FoldPreprocessor:
    """Median imputation and scaling learned within a training split."""

    def __init__(self) -> None:
        self.imputer = SimpleImputer(strategy="median")
        self.scaler = StandardScaler()

    def fit(self, values: np.ndarray) -> "FoldPreprocessor":
        imputed = self.imputer.fit_transform(values)
        self.scaler.fit(imputed)
        return self

    def transform(self, values: np.ndarray) -> np.ndarray:
        imputed = self.imputer.transform(values)
        return self.scaler.transform(imputed).astype(np.float32)


def set_thread_count(n_threads: int) -> None:
    for name in (
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ.setdefault(name, "1")
    os.environ["LIGHTGBM_NUM_THREADS"] = str(n_threads)


def stratification_labels(frame: pd.DataFrame, min_count: int) -> pd.Series:
    age_bin = pd.qcut(
        frame["Age"], q=min(10, frame["Age"].nunique()), labels=False, duplicates="drop"
    )
    candidates = (
        age_bin.astype(str)
        + "_"
        + frame["Sex_male"].astype(int).astype(str)
        + "_"
        + frame["Baseline_unhealthy"].astype(int).astype(str),
        age_bin.astype(str) + "_" + frame["Sex_male"].astype(int).astype(str),
        age_bin.astype(str),
    )
    for labels in candidates:
        if labels.value_counts().min() >= min_count:
            return labels.reset_index(drop=True)
    return candidates[-1].reset_index(drop=True)


def make_inner_folds(frame: pd.DataFrame, n_folds: int, seed: int) -> np.ndarray:
    labels = stratification_labels(frame, n_folds)
    splitter = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=seed)
    folds = np.empty(len(frame), dtype=np.int16)
    for fold, (_, valid) in enumerate(splitter.split(np.zeros(len(frame)), labels), start=1):
        folds[valid] = fold
    return folds


def tuning_sample(frame: pd.DataFrame, max_rows: int, seed: int) -> pd.DataFrame:
    if max_rows <= 0 or len(frame) <= max_rows:
        return frame.copy().reset_index(drop=True)
    labels = stratification_labels(frame, 2)
    splitter = StratifiedShuffleSplit(n_splits=1, train_size=max_rows, random_state=seed)
    keep, _ = next(splitter.split(np.zeros(len(frame)), labels))
    return frame.iloc[np.sort(keep)].copy().reset_index(drop=True)


def split_calibration_sample(
    development: pd.DataFrame, fraction: float, seed: int
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Reserve a stratified healthy reference before preprocessing and tuning."""
    age_bin = pd.qcut(
        development["Age"],
        q=min(10, development["Age"].nunique()),
        labels=False,
        duplicates="drop",
    )
    candidates = (
        development["AssessmentCenter"].astype(str)
        + "_"
        + age_bin.astype(str)
        + "_"
        + development["Sex_male"].astype(int).astype(str),
        age_bin.astype(str) + "_" + development["Sex_male"].astype(int).astype(str),
        age_bin.astype(str),
    )
    labels = next(
        (value for value in candidates if value.value_counts().min() >= 2),
        candidates[-1],
    )
    splitter = StratifiedShuffleSplit(n_splits=1, test_size=fraction, random_state=seed)
    fit_index, calibration_index = next(splitter.split(np.zeros(len(development)), labels))
    model_fit = development.iloc[np.sort(fit_index)].copy().reset_index(drop=True)
    calibration = development.iloc[np.sort(calibration_index)].copy().reset_index(drop=True)
    return model_fit, calibration


def make_xy(frame: pd.DataFrame, features: list[str]) -> tuple[np.ndarray, np.ndarray]:
    x = frame[features].to_numpy(dtype=np.float32)
    y = frame["Age"].to_numpy(dtype=np.float32)
    return x, y


def normal_nll(y: np.ndarray, mean: np.ndarray, sigma: np.ndarray) -> float:
    sigma = np.maximum(np.asarray(sigma, dtype=float), MIN_SIGMA)
    residual = np.asarray(y, dtype=float) - np.asarray(mean, dtype=float)
    return float(np.mean(0.5 * np.log(2 * np.pi * sigma**2) + residual**2 / (2 * sigma**2)))


def _distribution(mean: np.ndarray, sigma: np.ndarray) -> dict[str, np.ndarray]:
    sigma = np.maximum(np.asarray(sigma, dtype=np.float32), MIN_SIGMA)
    variance = np.maximum(sigma**2, MIN_VARIANCE).astype(np.float32)
    return {
        "BA_Raw": np.asarray(mean, dtype=np.float32),
        "Sigma_Raw": sigma,
        "Var_Raw": variance,
        "LogVar_Raw": np.log(variance).astype(np.float32),
    }


def default_parameters(model_name: str, n_threads: int) -> dict[str, Any]:
    if model_name == "ngboost":
        return {
            "n_estimators": 250,
            "learning_rate": 0.02,
            "minibatch_frac": 0.8,
            "col_sample": 0.8,
            "max_depth": 3,
            "min_samples_leaf": 20,
            "min_samples_split": 40,
            "natural_gradient": True,
        }
    if model_name == "xgboostlss":
        return {
            "eta": 0.03,
            "max_depth": 4,
            "min_child_weight": 5.0,
            "gamma": 0.01,
            "subsample": 0.85,
            "colsample_bytree": 0.85,
            "lambda": 1.0,
            "alpha": 0.1,
            "num_boost_round": 450,
            "tree_method": "hist",
            "nthread": n_threads,
            "verbosity": 0,
        }
    if model_name == "lightgbmlss":
        return {
            "learning_rate": 0.03,
            "max_depth": 5,
            "num_leaves": 32,
            "min_data_in_leaf": 20,
            "min_gain_to_split": 0.01,
            "bagging_fraction": 0.85,
            "bagging_freq": 1,
            "feature_fraction": 0.85,
            "lambda_l1": 0.1,
            "lambda_l2": 1.0,
            "num_boost_round": 450,
            "boosting": "gbdt",
            "num_threads": n_threads,
            "verbosity": -1,
            "feature_pre_filter": False,
        }
    raise ValueError(model_name)


def suggest_parameters(model_name: str, trial: optuna.Trial, n_threads: int) -> dict[str, Any]:
    if model_name == "ngboost":
        return {
            "n_estimators": trial.suggest_int("n_estimators", 120, 350),
            "learning_rate": trial.suggest_float("learning_rate", 0.005, 0.08, log=True),
            "minibatch_frac": trial.suggest_float("minibatch_frac", 0.55, 1.0),
            "col_sample": trial.suggest_float("col_sample", 0.55, 1.0),
            "max_depth": trial.suggest_int("max_depth", 2, 5),
            "min_samples_leaf": trial.suggest_int("min_samples_leaf", 5, 80),
            "min_samples_split": trial.suggest_int("min_samples_split", 10, 160),
            "natural_gradient": trial.suggest_categorical("natural_gradient", [True, False]),
        }
    if model_name == "xgboostlss":
        return {
            "eta": trial.suggest_float("eta", 0.005, 0.12, log=True),
            "max_depth": trial.suggest_int("max_depth", 2, 7),
            "min_child_weight": trial.suggest_float("min_child_weight", 1e-2, 100, log=True),
            "gamma": trial.suggest_float("gamma", 1e-8, 20, log=True),
            "subsample": trial.suggest_float("subsample", 0.6, 1.0),
            "colsample_bytree": trial.suggest_float("colsample_bytree", 0.5, 1.0),
            "lambda": trial.suggest_float("lambda", 1e-4, 100, log=True),
            "alpha": trial.suggest_float("alpha", 1e-4, 100, log=True),
            "num_boost_round": trial.suggest_int("num_boost_round", 200, 700),
            "tree_method": "hist",
            "nthread": n_threads,
            "verbosity": 0,
        }
    if model_name == "lightgbmlss":
        return {
            "learning_rate": trial.suggest_float("learning_rate", 0.005, 0.12, log=True),
            "max_depth": trial.suggest_int("max_depth", 2, 8),
            "num_leaves": trial.suggest_int("num_leaves", 8, 128),
            "min_data_in_leaf": trial.suggest_int("min_data_in_leaf", 5, 500),
            "min_gain_to_split": trial.suggest_float("min_gain_to_split", 1e-8, 10, log=True),
            "bagging_fraction": trial.suggest_float("bagging_fraction", 0.6, 1.0),
            "bagging_freq": trial.suggest_int("bagging_freq", 0, 5),
            "feature_fraction": trial.suggest_float("feature_fraction", 0.5, 1.0),
            "lambda_l1": trial.suggest_float("lambda_l1", 1e-4, 100, log=True),
            "lambda_l2": trial.suggest_float("lambda_l2", 1e-4, 100, log=True),
            "num_boost_round": trial.suggest_int("num_boost_round", 200, 700),
            "boosting": "gbdt",
            "num_threads": n_threads,
            "verbosity": -1,
            "feature_pre_filter": False,
        }
    raise ValueError(model_name)


def _search_parameters(model_name: str, n_threads: int) -> dict[str, Any]:
    fixed = {"tree_method", "nthread", "verbosity", "boosting", "num_threads", "feature_pre_filter"}
    return {k: v for k, v in default_parameters(model_name, n_threads).items() if k not in fixed}


def fit_distribution_model(
    model_name: str,
    x: np.ndarray,
    y: np.ndarray,
    parameters: dict[str, Any],
    seed: int,
):
    if model_name == "ngboost":
        from ngboost import NGBRegressor
        from ngboost.distns import Normal
        from ngboost.scores import MLE
        from sklearn.tree import DecisionTreeRegressor

        learner = DecisionTreeRegressor(
            max_depth=int(parameters["max_depth"]),
            min_samples_leaf=int(parameters["min_samples_leaf"]),
            min_samples_split=int(parameters["min_samples_split"]),
            random_state=seed,
        )
        model = NGBRegressor(
            Dist=Normal,
            Score=MLE,
            Base=learner,
            n_estimators=int(parameters["n_estimators"]),
            learning_rate=float(parameters["learning_rate"]),
            minibatch_frac=float(parameters["minibatch_frac"]),
            col_sample=float(parameters["col_sample"]),
            natural_gradient=bool(parameters["natural_gradient"]),
            random_state=seed,
            verbose=False,
        )
        return model.fit(x, y)

    if model_name == "xgboostlss":
        import xgboost as xgb
        from xgboostlss.distributions.Gaussian import Gaussian
        from xgboostlss.model import XGBoostLSS

        model = XGBoostLSS(Gaussian(stabilization="None", response_fn="exp", loss_fn="nll"))
        train = xgb.DMatrix(x, label=y)
        train_parameters = {k: v for k, v in parameters.items() if k != "num_boost_round"}
        train_parameters["seed"] = seed
        model.train(train_parameters, train, num_boost_round=int(parameters["num_boost_round"]))
        return model

    if model_name == "lightgbmlss":
        import lightgbm as lgb
        from lightgbmlss.distributions.Gaussian import Gaussian
        from lightgbmlss.model import LightGBMLSS

        model = LightGBMLSS(
            Gaussian(stabilization="None", response_fn="exp", loss_fn="nll", initialize=True)
        )
        train = lgb.Dataset(pd.DataFrame(x), label=y, free_raw_data=False)
        train_parameters = {k: v for k, v in parameters.items() if k != "num_boost_round"}
        train_parameters["seed"] = seed
        train_parameters["feature_pre_filter"] = False
        model.train(train_parameters, train, num_boost_round=int(parameters["num_boost_round"]))
        return model

    raise ValueError(model_name)


def predict_distribution(model_name: str, model: Any, x: np.ndarray) -> dict[str, np.ndarray]:
    if model_name == "ngboost":
        predicted = model.pred_dist(x)
        return _distribution(predicted.loc, predicted.scale)
    if model_name == "xgboostlss":
        import xgboost as xgb

        predicted = model.predict(xgb.DMatrix(x), pred_type="parameters")
    elif model_name == "lightgbmlss":
        predicted = model.predict(pd.DataFrame(x), pred_type="parameters")
    else:
        raise ValueError(model_name)
    if not isinstance(predicted, pd.DataFrame):
        predicted = pd.DataFrame(predicted)
    mean_name = next(name for name in ("loc", "mu", "mean") if name in predicted)
    sigma_name = next(name for name in ("scale", "sigma", "sd") if name in predicted)
    return _distribution(predicted[mean_name].to_numpy(), predicted[sigma_name].to_numpy())


def tune_model(
    model_name: str,
    frame: pd.DataFrame,
    features: list[str],
    folds: np.ndarray,
    n_trials: int,
    seed: int,
    n_threads: int,
    run_dir: Path,
) -> tuple[dict[str, Any], pd.DataFrame]:
    if n_trials == 0:
        return default_parameters(model_name, n_threads), pd.DataFrame()
    x, y = make_xy(frame, features)
    study = optuna.create_study(direction="minimize", sampler=optuna.samplers.TPESampler(seed=seed))
    study.enqueue_trial(_search_parameters(model_name, n_threads))

    def objective(trial: optuna.Trial) -> float:
        parameters = suggest_parameters(model_name, trial, n_threads)
        scores: list[float] = []
        for fold in sorted(np.unique(folds)):
            try:
                train_index = np.flatnonzero(folds != fold)
                valid_index = np.flatnonzero(folds == fold)
                preprocessor = FoldPreprocessor().fit(x[train_index])
                model = fit_distribution_model(
                    model_name,
                    preprocessor.transform(x[train_index]),
                    y[train_index],
                    parameters,
                    seed + trial.number * 1000 + int(fold),
                )
                predicted = predict_distribution(
                    model_name, model, preprocessor.transform(x[valid_index])
                )
                score = normal_nll(y[valid_index], predicted["BA_Raw"], predicted["Sigma_Raw"])
                if not np.isfinite(score):
                    raise ValueError("Non-finite validation NLL")
            except Exception as exc:
                failure = f"fold {int(fold)}: {type(exc).__name__}: {str(exc)[:500]}"
                trial.set_user_attr("failure", failure)
                LOGGER.warning("%s trial %d failed on %s", model_name, trial.number, failure)
                return FAILED_TRIAL_NLL
            scores.append(score)
        return float(np.mean(scores))

    def save_trials(study_: optuna.Study, _: optuna.trial.FrozenTrial) -> None:
        study_.trials_dataframe().to_csv(run_dir / "tuning_trials.csv", index=False)

    study.optimize(objective, n_trials=n_trials, gc_after_trial=True, callbacks=[save_trials])
    trials = study.trials_dataframe()
    successful = [
        trial
        for trial in study.trials
        if trial.state == optuna.trial.TrialState.COMPLETE
        and not trial.user_attrs.get("failure")
        and trial.value is not None
        and np.isfinite(trial.value)
    ]
    if not successful:
        raise RuntimeError(f"No successful {model_name} tuning trial")
    best = min(successful, key=lambda trial: trial.value)
    parameters = {**default_parameters(model_name, n_threads), **best.params}
    return parameters, trials


def fit_and_predict(
    model_name: str,
    training: pd.DataFrame,
    scoring: pd.DataFrame,
    features: list[str],
    parameters: dict[str, Any],
    seed: int,
    model_dir: Path,
) -> tuple[pd.DataFrame, FoldPreprocessor]:
    x_train, y_train = make_xy(training, features)
    x_score, _ = make_xy(scoring, features)
    preprocessor = FoldPreprocessor().fit(x_train)
    model = fit_distribution_model(
        model_name, preprocessor.transform(x_train), y_train, parameters, seed
    )
    predicted = predict_distribution(model_name, model, preprocessor.transform(x_score))
    columns = [
        "eid_int",
        "Age",
        "Sex",
        "Sex_male",
        "Baseline_unhealthy",
        "AssessmentCenter",
        "SplitRole",
    ]
    output = scoring[columns].copy()
    for name, values in predicted.items():
        output[name] = values
    output["AgeGap_Raw"] = output["BA_Raw"] - output["Age"]
    model_dir.mkdir(parents=True, exist_ok=True)
    joblib.dump(
        {
            "model": model,
            "preprocessor": preprocessor,
            "features": features,
            "parameters": parameters,
            "model_name": model_name,
        },
        model_dir / "final_model.joblib",
    )
    return output, preprocessor


def point_metrics(frame: pd.DataFrame, mean_column: str = "BA_Calibrated") -> dict[str, float]:
    age = frame["Age"].to_numpy(float)
    predicted = frame[mean_column].to_numpy(float)
    return {
        "n": int(len(frame)),
        "mae": float(mean_absolute_error(age, predicted)),
        "rmse": float(mean_squared_error(age, predicted) ** 0.5),
        "pearson_r": float(np.corrcoef(age, predicted)[0, 1]),
    }
