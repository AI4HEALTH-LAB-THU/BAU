"""Failed parameter candidates leave the remaining tuning trials available."""

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd

from bau.modeling import FAILED_TRIAL_NLL, tune_model


class TuningTests(unittest.TestCase):
    def setUp(self):
        self.frame = pd.DataFrame({"measurement": np.arange(12), "Age": np.arange(12) + 50})
        self.folds = np.tile([1, 2], 6)
        self.prediction = {"BA_Raw": np.full(6, 55.0), "Sigma_Raw": np.ones(6)}

    def tune(self, directory, n_trials=2):
        return tune_model(
            "ngboost", self.frame, ["measurement"], self.folds, n_trials, 19, 1, Path(directory)
        )

    def test_failed_candidate_preserves_later_trials(self):
        for failure in (RuntimeError("candidate failed"), None):
            with (
                self.subTest(failure=type(failure).__name__),
                tempfile.TemporaryDirectory() as directory,
            ):
                if failure is None:
                    predictions = [
                        {"BA_Raw": np.full(6, np.nan), "Sigma_Raw": np.ones(6)},
                        self.prediction,
                        self.prediction,
                    ]
                    fit_effect = None
                else:
                    predictions = [self.prediction, self.prediction]
                    fit_effect = [failure, object(), object()]
                with (
                    patch("bau.modeling.fit_distribution_model", side_effect=fit_effect),
                    patch("bau.modeling.predict_distribution", side_effect=predictions),
                ):
                    parameters, trials = self.tune(directory)
                self.assertEqual(len(trials), 2)
                self.assertEqual(trials.loc[0, "value"], FAILED_TRIAL_NLL)
                self.assertLess(trials.loc[1, "value"], FAILED_TRIAL_NLL)
                self.assertIn("max_depth", parameters)
                self.assertTrue((Path(directory) / "tuning_trials.csv").is_file())

    def test_failed_candidate_cannot_win_against_large_finite_nll(self):
        prediction = {"BA_Raw": np.full(6, 1000.0), "Sigma_Raw": np.full(6, 1e-6)}
        with (
            tempfile.TemporaryDirectory() as directory,
            patch(
                "bau.modeling.fit_distribution_model",
                side_effect=[RuntimeError("candidate failed"), object(), object()],
            ),
            patch("bau.modeling.predict_distribution", side_effect=[prediction, prediction]),
        ):
            parameters, trials = self.tune(directory)
        self.assertGreater(trials.loc[1, "value"], FAILED_TRIAL_NLL)
        self.assertTrue(pd.notna(trials.loc[0, "user_attrs_failure"]))
        for column in trials.columns[trials.columns.str.startswith("params_")]:
            self.assertEqual(parameters[column.removeprefix("params_")], trials.loc[1, column])

    def test_all_failed_candidates_report_failure(self):
        with (
            tempfile.TemporaryDirectory() as directory,
            patch(
                "bau.modeling.fit_distribution_model", side_effect=RuntimeError("candidate failed")
            ),
        ):
            with self.assertRaisesRegex(RuntimeError, "No successful ngboost tuning trial"):
                self.tune(directory)


if __name__ == "__main__":
    unittest.main()
