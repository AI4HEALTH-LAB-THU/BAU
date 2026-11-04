#!/usr/bin/env python3
"""Build participant-aligned UK Biobank covariate and outcome tables."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.feather as feather
import pyarrow.ipc as ipc

from bau.paths import CODE_ROOT as ROOT
from bau.paths import DATA_ROOT

UKB_ROOT = Path(os.environ.get("UKB_DATA_ROOT", "/mnt/project/UKB"))
H7_ROOT = Path(os.environ.get("UKB_H7_ROOT", "/mnt/project/UKB/Derived/Composite_Outcomes_H7_H7MM"))
CONFIG = ROOT / "config/supporting.json"
OUTPUT = DATA_ROOT / "ukb/support"

HISTORY_COLUMNS = [
    "angina_pectoris_history",
    "myocardial_infarction_history",
    "other_ischaemic_heart_diseases_history",
    "stroke_history",
    "chronic_obstructive_pulmonary_disease_history",
    "diabetes_mellitus_history",
    "nephrotic_syndrome_history",
    "renal_failure_history",
    "parkinsonism_history",
    "dementia_history",
    "depressive_episode_history",
    "cancer_history",
    "hypertension_history",
    "hypertensive_heart_and_renal_diseases_history",
    "heart_failure_history",
    "respiratory_failure_history",
    "hepatic_failure_history",
    "fibrosis_and_cirrhosis_of_liver_history",
    "multiple_sclerosis_and_other_demyelinating_diseases_history",
    "schizophrenia_and_schizoaffective_disorders_history",
    "rheumatoid_arthritis_history",
    "other_arthritis_history",
    "osteoporosis_history",
]


def read_eids(path: Path) -> np.ndarray:
    values = pd.read_feather(path, columns=["eid"])["eid"]
    eids = pd.to_numeric(values, errors="raise").round().astype("int64").to_numpy()
    if len(np.unique(eids)) != len(eids):
        raise ValueError(f"Duplicated participant identifier in {path}")
    return eids


def write_feather(frame: pd.DataFrame, path: Path) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    frame.to_feather(temporary, compression="zstd")
    os.replace(temporary, path)


def rewrite_feather(source: Path, target: Path, eids: np.ndarray) -> None:
    temporary = target.with_suffix(target.suffix + ".tmp")
    with source.open("rb") as source_handle:
        reader = ipc.open_file(source_handle)
        eid_column = reader.schema.get_field_index("eid")
        if eid_column < 0:
            raise ValueError(f"Missing eid column in {source}")
        fields = list(reader.schema)
        fields[eid_column] = pa.field("eid", pa.int64())
        schema = pa.schema(fields)
        options = ipc.IpcWriteOptions(compression="zstd")
        offset = 0
        with (
            temporary.open("wb") as target_handle,
            ipc.new_file(target_handle, schema, options=options) as writer,
        ):
            for batch_number in range(reader.num_record_batches):
                batch = reader.get_batch(batch_number)
                observed = np.rint(batch.column(eid_column).to_numpy(zero_copy_only=False)).astype(
                    "int64"
                )
                expected = eids[offset : offset + batch.num_rows]
                if not np.array_equal(observed, expected):
                    raise ValueError(f"Participant order differs in {source.name}")
                columns = list(batch.columns)
                columns[eid_column] = pa.array(observed, type=pa.int64())
                writer.write_batch(pa.record_batch(columns, schema=schema))
                offset += batch.num_rows
    if offset != len(eids):
        raise ValueError(f"Participant count differs in {source.name}")
    os.replace(temporary, target)


def build_participant_index(
    population_path: Path, history: pd.DataFrame, candidate_dir: Path
) -> pd.DataFrame:
    population = pd.read_csv(
        population_path,
        usecols=[
            "eid",
            "Date of birth",
            "Date of attending assessment centre",
            "UK Biobank assessment centre",
            "Sex",
        ],
    )
    population["eid"] = pd.to_numeric(population["eid"], errors="raise").astype("int64")
    if population["eid"].duplicated().any():
        raise ValueError("Duplicated participant identifier in population table")
    out = history[["eid", "Baseline_unhealthy"]].merge(
        population, on="eid", how="left", validate="one_to_one"
    )
    birth = pd.to_datetime(out.pop("Date of birth"), errors="coerce")
    visit = pd.to_datetime(out.pop("Date of attending assessment centre"), errors="coerce")
    out["Age"] = ((visit - birth).dt.days / 365.2425).round(1).astype("float32")
    out["Sex"] = out["Sex"].map({0: "female", 1: "male"})
    out["Sex_male"] = (out["Sex"] == "male").astype("int8")
    out["AssessmentCenter"] = pd.to_numeric(
        out.pop("UK Biobank assessment centre"), errors="coerce"
    ).astype("int32")
    out["Date_birth"] = birth.dt.strftime("%Y-%m-%d")
    out["Date_attend"] = visit.dt.strftime("%Y-%m-%d")
    required = ["Age", "Sex", "AssessmentCenter", "Date_birth", "Date_attend"]
    if out[required].isna().any().any():
        raise ValueError("Incomplete population metadata")

    candidates = {
        "in_nmr": candidate_dir / "nmr.feather",
        "in_protein": candidate_dir / "protein.feather",
        "in_clinical": candidate_dir / "routine_blood.feather",
    }
    for column, path in candidates.items():
        ids = set(read_eids(path))
        out[column] = out["eid"].isin(ids).astype("int8")
    out["modality_count"] = out[list(candidates)].sum(axis=1).astype("int8")
    return out[
        [
            "eid",
            "Age",
            "Sex",
            "Sex_male",
            "AssessmentCenter",
            "Date_birth",
            "Date_attend",
            "Baseline_unhealthy",
            *candidates,
            "modality_count",
        ]
    ]


def check_h7(path: Path, expected_eids: np.ndarray, columns: list[str]) -> None:
    table = feather.read_table(path, columns=columns)
    frame = table.to_pandas()
    observed = pd.to_numeric(frame["eid"], errors="raise").astype("int64").to_numpy()
    if not np.array_equal(observed, expected_eids):
        raise ValueError(f"Participant order differs in {path.name}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    settings = json.loads(CONFIG.read_text(encoding="utf-8"))
    sources = {name: UKB_ROOT / path for name, path in settings["sources"].items()}
    h7 = {name: H7_ROOT / path for name, path in settings["h7_files"].items()}
    OUTPUT.mkdir(parents=True, exist_ok=True)
    targets = {
        "participants": OUTPUT / "participants.feather",
        "covariates": OUTPUT / "covariates.feather",
        "outcomes_primary": OUTPUT / "outcomes_primary.feather",
        "disease_history": OUTPUT / "disease_history.feather",
        "h7_events": OUTPUT / "h7_events.feather",
        "h7_baseline": OUTPUT / "h7_baseline.feather",
    }
    if not args.overwrite:
        existing = [path for path in targets.values() if path.exists()]
        if existing:
            raise FileExistsError(f"Supporting tables already exist: {existing}")

    history = pd.read_feather(sources["disease_history"], columns=["eid", *HISTORY_COLUMNS])
    history["eid"] = pd.to_numeric(history["eid"], errors="raise").astype("int64")
    if not history[HISTORY_COLUMNS].isin([0, 1]).all().all():
        raise ValueError("Disease history values must be binary")
    history["Baseline_unhealthy"] = history[HISTORY_COLUMNS].max(axis=1).astype("int8")
    eids = history["eid"].to_numpy()
    expected_n = int(settings["expected_participants"])
    if len(eids) != expected_n or len(np.unique(eids)) != expected_n:
        raise ValueError(f"Expected {expected_n} participants, found {len(eids)}")

    participants = build_participant_index(
        sources["population"], history, DATA_ROOT / "ukb/modeling"
    )
    write_feather(participants, targets["participants"])
    write_feather(history, targets["disease_history"])
    rewrite_feather(sources["outcomes_primary"], targets["outcomes_primary"], eids)

    export_script = Path(__file__).with_name("covariates.R")
    subprocess.run(
        [
            "Rscript",
            str(export_script),
            str(targets["participants"]),
            str(sources["covariates"]),
            str(targets["covariates"]),
        ],
        check=True,
    )
    check_h7(
        h7["events"],
        eids,
        [
            "eid",
            "h7_healthspan_diagnose",
            "h7_healthspan_duration",
            "h7_multimorbidity_diagnose",
            "h7_multimorbidity_duration",
        ],
    )
    check_h7(h7["baseline"], eids, ["eid", "baseline_h7_domain_count"])
    shutil.copyfile(h7["events"], targets["h7_events"])
    shutil.copyfile(h7["baseline"], targets["h7_baseline"])

    summary = {
        "participants": len(participants),
        "healthy": int((participants["Baseline_unhealthy"] == 0).sum()),
        "nmr": int(participants["in_nmr"].sum()),
        "proteomic": int(participants["in_protein"].sum()),
        "routine_blood": int(participants["in_clinical"].sum()),
    }
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
