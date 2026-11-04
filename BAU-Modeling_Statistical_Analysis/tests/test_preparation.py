"""Cohort eligibility and participant alignment on synthetic measurements."""

import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd

from analysis.prepare import supporting as support
from analysis.prepare import ukb as cohorts
from bau.paths import RESULTS_ROOT


class PreparationTests(unittest.TestCase):
    def setUp(self):
        self.population = pd.DataFrame(
            {
                "eid": np.arange(1, 11),
                "Age": 55.0,
                "Sex": "female",
                "Sex_male": 0,
                "AssessmentCenter": 11012,
                "Date_birth": "1950-01-01",
                "Date_attend": "2005-01-01",
            }
        )
        self.features = [f"measurement_{i}" for i in range(10)]
        self.assay = pd.DataFrame(np.ones((10, 10)), columns=self.features)
        self.assay.insert(0, "eid", self.population.eid)
        self.history = pd.DataFrame({"eid": self.population.eid, "Baseline_unhealthy": 0})
        test_root = Path(os.environ.get("BAU_TEST_ROOT", RESULTS_ROOT / "tests"))
        test_root.mkdir(parents=True, exist_ok=True)
        self.temporary = tempfile.TemporaryDirectory(dir=test_root)
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)

    def filter_assay(self):
        return cohorts.attach_metadata_and_filter(
            "nmr", self.population, self.assay, self.features, self.history, 10
        )

    def test_row_missingness_threshold(self):
        self.assay.loc[0, self.features[:2]] = np.nan
        result, _, _, statistics = self.filter_assay()
        self.assertEqual(result.eid.tolist(), list(range(2, 11)))
        self.assertEqual(statistics["excluded_participant_missing_ge_20pct"], 1)

    def test_feature_missingness_threshold(self):
        self.assay.loc[:1, self.features[0]] = np.nan
        result, _, features, _ = self.filter_assay()
        self.assertEqual(len(result), 10)
        self.assertNotIn(self.features[0], features)
        self.assertEqual(features, self.features[1:])

    def test_age_filter_precedes_rounding(self):
        source = self.directory / "population.csv"
        pd.DataFrame(
            {
                "eid": [1, 2],
                "Date of birth": ["1950-01-01"] * 2,
                "Date of attending assessment centre": ["1989-12-31", "1990-01-01"],
                "UK Biobank assessment centre": [11012] * 2,
                "Age when attended assessment centre": [40, 40],
                "Sex": [0, 0],
            }
        ).to_csv(source, index=False)
        with patch.object(cohorts, "POPULATION_PATH", source):
            result, _ = cohorts.load_population()
        self.assertEqual(result.eid.tolist(), [2])
        self.assertAlmostEqual(result.Age.iloc[0], 40.0)

    def test_source_and_age_eligible_counts(self):
        _, _, _, statistics = cohorts.attach_metadata_and_filter(
            "nmr", self.population, self.assay, self.features, self.history, 12
        )
        self.assertEqual(statistics["source_rows_with_any_measurement"], 12)
        self.assertEqual(statistics["age_eligible_with_any_measurement"], 10)

    def test_duplicate_assay_identifiers(self):
        self.assay.loc[1, "eid"] = 1
        with self.assertRaises(ValueError):
            self.filter_assay()

    def test_support_table_order(self):
        source = self.directory / "source.feather"
        target = self.directory / "aligned.feather"
        pd.DataFrame({"eid": [1, 2, 3], "value": [0.1, 0.2, 0.3]}).to_feather(source)
        support.rewrite_feather(source, target, np.array([1, 2, 3], dtype=np.int64))
        self.assertEqual(pd.read_feather(target).eid.tolist(), [1, 2, 3])
        with self.assertRaisesRegex(ValueError, "order"):
            support.rewrite_feather(source, target, np.array([3, 2, 1], dtype=np.int64))


if __name__ == "__main__":
    unittest.main()
