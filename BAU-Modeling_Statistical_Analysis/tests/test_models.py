import unittest

import numpy as np
import pandas as pd

from bau.calibration import (
    apply_calibration_bundle,
    apply_predictive_calibration,
    fit_calibration_bundle,
)
from bau.modeling import (
    FoldPreprocessor,
    default_parameters,
    fit_distribution_model,
    make_inner_folds,
    normal_nll,
    predict_distribution,
    split_calibration_sample,
)


class CalibrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        rng = np.random.default_rng(19)
        n = 1200
        sex = rng.integers(0, 2, n)
        age = rng.uniform(40, 70, n)
        sigma = np.exp(rng.normal(0.9, 0.2, n))
        cls.raw = pd.DataFrame(
            {
                "Age": age,
                "Sex_male": sex,
                "BA_Raw": 0.9 * age + 5 + rng.normal(0, sigma),
                "Var_Raw": sigma**2,
            }
        )
        config = dict(
            mean_crossfit_folds=5,
            min_sex_n=20,
            min_variance_factor=0.01,
            max_variance_factor=100,
            phenotype_lowess_frac=0.35,
        )
        cls.bundle, cls.reference = fit_calibration_bundle(cls.raw, config, 7)

    def test_reference_scale(self):
        for _, group in self.reference.groupby("Sex_male"):
            for term in ("BAG", "BAU"):
                self.assertAlmostEqual(group[term].mean(), 0, places=12)
                self.assertAlmostEqual(group[term].std(ddof=0), 1, places=12)

    def test_predictive_calibration_does_not_use_age(self):
        shifted = self.raw.copy()
        shifted["Age"] += 8
        calibrator = self.bundle["predictive_calibrator"]
        a = apply_predictive_calibration(self.raw, calibrator)
        b = apply_predictive_calibration(shifted, calibrator)
        np.testing.assert_array_equal(a.BA_Calibrated, b.BA_Calibrated)
        np.testing.assert_array_equal(a.Var_Calibrated, b.Var_Calibrated)

    def test_application_is_independent_of_batch(self):
        whole = apply_calibration_bundle(self.raw, self.bundle)
        single = apply_calibration_bundle(self.raw.iloc[[47]], self.bundle)
        np.testing.assert_allclose(
            whole.iloc[47][["BAG", "BAU"]].to_numpy(float),
            single.iloc[0][["BAG", "BAU"]].to_numpy(float),
        )

    def test_fitted_reference_reproduces(self):
        applied = apply_calibration_bundle(self.raw, self.bundle)
        np.testing.assert_allclose(applied[["BAG", "BAU"]], self.reference[["BAG", "BAU"]])


class ModelingTests(unittest.TestCase):
    def test_preprocessing_uses_training_rows(self):
        train = np.array([[1, 3], [3, np.nan], [5, 7]], dtype=float)
        p = FoldPreprocessor().fit(train)
        before = p.imputer.statistics_.copy(), p.scaler.mean_.copy()
        p.transform(np.array([[900, -900], [np.nan, 0]]))
        np.testing.assert_array_equal(p.imputer.statistics_, before[0])
        np.testing.assert_array_equal(p.scaler.mean_, before[1])
        np.testing.assert_allclose(p.transform(train).mean(axis=0), 0, atol=1e-6)

    def test_disjoint_calibration_split(self):
        rng = np.random.default_rng(23)
        frame = pd.DataFrame(
            {
                "eid": np.arange(2000),
                "Age": rng.uniform(40, 70, 2000),
                "Sex_male": rng.integers(0, 2, 2000),
                "AssessmentCenter": rng.integers(0, 4, 2000),
                "Baseline_unhealthy": 0,
            }
        )
        fit, reference = split_calibration_sample(frame, 0.2, 31)
        self.assertEqual(len(reference), 400)
        self.assertFalse(set(fit.eid) & set(reference.eid))
        self.assertEqual(set(fit.eid) | set(reference.eid), set(frame.eid))
        self.assertEqual(set(make_inner_folds(fit, 5, 13)), set(range(1, 6)))

    def test_gaussian_likelihood(self):
        self.assertAlmostEqual(
            normal_nll(np.array([0.0]), np.array([0.0]), np.array([1.0])), 0.5 * np.log(2 * np.pi)
        )

    def test_distribution_models(self):
        rng = np.random.default_rng(8)
        x = rng.normal(size=(160, 5)).astype(np.float32)
        y = (55 + 3 * x[:, 0] + rng.normal(size=160)).astype(np.float32)
        for name in ("ngboost", "lightgbmlss", "xgboostlss"):
            with self.subTest(model=name):
                params = default_parameters(name, 1)
                params["n_estimators" if name == "ngboost" else "num_boost_round"] = 4
                model = fit_distribution_model(name, x, y, params, 9)
                predicted = predict_distribution(name, model, x[:7])
                self.assertEqual(predicted["BA_Raw"].shape, (7,))
                self.assertTrue(all(np.isfinite(a).all() for a in predicted.values()))
                self.assertTrue((predicted["Var_Raw"] > 0).all())
                np.testing.assert_allclose(
                    predicted["Sigma_Raw"] ** 2, predicted["Var_Raw"], rtol=1e-6
                )


if __name__ == "__main__":
    unittest.main()
