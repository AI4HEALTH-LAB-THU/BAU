"""Partition the healthy reference for proteomic organ clocks."""

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

from bau.modeling import split_calibration_sample
from bau.paths import CODE_ROOT as ROOT
from bau.paths import RESULTS_ROOT

UKB = Path(os.environ.get("UKB_DATA_ROOT", "/mnt/project/UKB"))
REFERENCE = Path(os.environ.get("DOMAIN_REFERENCE_ROOT", UKB / "Derived/DomainClockReference"))


def main():
    config = json.loads((ROOT / "config/organ_clocks.json").read_text())
    folds = json.loads((ROOT / "config/center_folds.json").read_text())
    centres = pd.read_csv(
        UKB / "ProcessedData/02-AssessmentCenter/01-Recruitment/basic_population_info_0.csv",
        usecols=["eid", "UK Biobank assessment centre"],
    ).rename(columns={"UK Biobank assessment centre": "AssessmentCenter"})
    centres["eid"] = centres.eid.astype("int64")
    # Row order determines the stratified partitions.
    x = pd.read_feather(REFERENCE / "protein_healthy.feather", columns=["eid", "Age", "Sex"])
    x["eid"] = x.eid.astype("int64")
    x["Sex"] = x.Sex.astype(str).str.strip().str.lower()
    x = x.loc[x.Age.between(37, 73) & x.Sex.isin(["female", "male"])].copy()
    x = x.merge(centres, on="eid", how="left", sort=False, validate="one_to_one")
    x["AssessmentCenter"] = x.AssessmentCenter.astype(int).astype(str)
    x["Sex_male"] = x.Sex.eq("male").astype(np.float32)
    x["eid_int"] = x.eid
    outer = x.AssessmentCenter.map(folds["outer_center_to_fold"])
    if len(x) != 22731 or outer.isna().any():
        raise ValueError("Unexpected protein domain-reference population")
    for fold in range(1, 6):
        development = x.loc[outer.ne(fold)].copy().reset_index(drop=True)
        fit, calibration = split_calibration_sample(
            development, 0.2, config["reference_seed"] + 100 * fold + 5000
        )
        roles = pd.concat(
            [
                fit[["eid_int"]].assign(role="model_fit"),
                calibration[["eid_int"]].assign(role="calibration_reference"),
            ],
            ignore_index=True,
        )
        output = RESULTS_ROOT / "phenotypes/domain_reference" / "protein"
        output.mkdir(parents=True, exist_ok=True)
        roles.to_csv(output / f"fold_{fold}.csv", index=False)


if __name__ == "__main__":
    main()
