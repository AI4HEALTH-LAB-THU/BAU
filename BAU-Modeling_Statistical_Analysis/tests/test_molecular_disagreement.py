"""Synthetic checks for the observed-biomarker likelihood loss."""

import unittest

import numpy as np

from analysis.phenotypes.molecular_disagreement import AGE_GRID, disagreement


class MolecularDisagreementTests(unittest.TestCase):
    def test_shared_peak_and_conflict(self):
        common = np.stack([-((AGE_GRID - 55) ** 2), -((AGE_GRID - 55) ** 2)])[None]
        conflict = np.stack([-((AGE_GRID - 45) ** 2), -((AGE_GRID - 65) ** 2)])[None]
        self.assertEqual(disagreement(common)[0][0], 0)
        self.assertGreater(disagreement(conflict)[0][0], 0)
        self.assertEqual(AGE_GRID[disagreement(conflict)[1][0]], 55)
        expected = conflict.max(axis=2).mean(axis=1) - conflict.mean(axis=1).max(axis=1)
        np.testing.assert_allclose(disagreement(conflict)[0], expected)

    def test_missing_and_vertical_offsets(self):
        curves = np.stack([-((AGE_GRID - 45) ** 2), -((AGE_GRID - 65) ** 2)])[None]
        with_missing = np.concatenate([curves, np.full((1, 1, len(AGE_GRID)), np.nan)], axis=1)
        np.testing.assert_allclose(disagreement(with_missing)[0], disagreement(curves)[0])
        np.testing.assert_allclose(
            disagreement(curves + np.array([3, -10])[None, :, None])[0], disagreement(curves)[0]
        )
        with self.assertRaises(ValueError):
            disagreement(np.full((1, 2, len(AGE_GRID)), np.nan))


if __name__ == "__main__":
    unittest.main()
