#!/usr/bin/env python3
"""Prepare the ten proteomic organ panels used for the domain clocks."""

from __future__ import annotations

import argparse
import json

import pandas as pd

from analysis.prepare import ukb as cohorts
from bau.paths import CODE_ROOT as ROOT
from bau.paths import DATA_ROOT, RESULTS_ROOT

FEATURE_DIR = ROOT / "config/organ_features"
OUTPUT_DIR = DATA_ROOT / "ukb/domain_clocks"
REPORT_DIR = RESULTS_ROOT / "preparation/ukb"
EXPECTED_ROWS = {
    "artery": 43_501,
    "brain": 43_500,
    "heart": 38_089,
    "immune": 43_777,
    "intestine": 43_645,
    "kidney": 43_396,
    "liver": 43_661,
    "lung": 49_866,
    "muscle": 40_817,
    "pancreas": 43_613,
}


def read_features(panel: str) -> list[str]:
    path = FEATURE_DIR / f"{panel}.txt"
    features = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    if not features or len(features) != len(set(features)):
        raise ValueError(f"Invalid feature list for {panel}")
    return features


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    existing = [OUTPUT_DIR / f"{panel}.feather" for panel in EXPECTED_ROWS]
    if not args.overwrite and any(path.exists() for path in existing):
        raise FileExistsError("Organ-panel data already exist")

    population, _ = cohorts.load_population()
    history = cohorts.load_history()
    feature_sets = {panel: read_features(panel) for panel in EXPECTED_ROWS}
    union = sorted(set().union(*feature_sets.values()))
    raw = cohorts.read_numeric_csv(cohorts.PROTEIN_PATH, ["eid", *union])

    records = []
    for panel, features in feature_sets.items():
        cohort, profile, retained, statistics = cohorts.attach_metadata_and_filter(
            panel,
            population,
            raw[["eid", *features]],
            features,
            history,
            int(raw[features].notna().any(axis=1).sum()),
        )
        if retained != features:
            raise ValueError(f"Feature missingness changed the {panel} panel")
        if len(cohort) != EXPECTED_ROWS[panel]:
            raise ValueError(
                f"Expected {EXPECTED_ROWS[panel]} {panel} participants, found {len(cohort)}"
            )
        cohorts.atomic_feather(cohort, OUTPUT_DIR / f"{panel}.feather")
        profile.to_csv(REPORT_DIR / f"{panel}_features.tsv", sep="\t", index=False)
        records.append(statistics)

    pd.DataFrame(records).to_csv(REPORT_DIR / "organ_panel_flow.tsv", sep="\t", index=False)
    print(json.dumps({row["modality"]: row["final_participants"] for row in records}, indent=2))


if __name__ == "__main__":
    main()
