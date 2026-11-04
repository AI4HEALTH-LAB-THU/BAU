#!/usr/bin/env python3
"""Prepare the NMR, proteomic and routine blood modeling cohorts."""

from __future__ import annotations

import argparse
import csv
import os
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from bau.paths import DATA_ROOT, RESULTS_ROOT

UKB_ROOT = Path(os.environ.get("UKB_DATA_ROOT", "/mnt/project/UKB"))
PROCESSED = UKB_ROOT / "ProcessedData"
POPULATION_PATH = PROCESSED / "01-PopulationCharacteristics/basic_population_info_0.csv"
NMR_PATH = PROCESSED / "03-BiologicalSamples/NMR_instance_original_0.csv"
PROTEIN_PATH = PROCESSED / "03-BiologicalSamples/ukb_protein_full.csv"
HISTORY_PATH = UKB_ROOT / "ResearchData/Aging_Indicator/dat_ukb_aging_indicator_medhis_0.feather"

DATA_DIR = DATA_ROOT / "ukb/modeling"
REPORT_DIR = RESULTS_ROOT / "preparation/ukb"

MISSING_THRESHOLD = 0.20
FEATURE_MISSINGNESS_EXCEPTIONS = {
    "clinical": {"Lipoprotein A": 0.2005},
}
AGE_MIN = 40.0
AGE_MAX = 70.0
METADATA_COLUMNS = [
    "eid",
    "Age",
    "Sex",
    "Sex_male",
    "AssessmentCenter",
    "Date_birth",
    "Date_attend",
    "Baseline_unhealthy",
]

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

PROTEIN_FEATURES = [
    "ACRV1",
    "ACTA2",
    "ADA",
    "ADAMTS13",
    "ADAMTS15",
    "ADAMTS16",
    "ADAMTSL5",
    "ADGRG1",
    "AFP",
    "AGER",
    "AGRP",
    "AHNAK2",
    "ANGPT2",
    "BAG3",
    "BCAN",
    "BGLAP",
    "BOC",
    "BSG",
    "C19ORF12",
    "C1QL2",
    "CA14",
    "CA4",
    "CALB1",
    "CCDC80",
    "CCL28",
    "CCN5",
    "CD1C",
    "CD248",
    "CD8A",
    "CD93",
    "CDCP1",
    "CDH2",
    "CDH3",
    "CDHR2",
    "CDON",
    "CELSR2",
    "CFHR5",
    "CHGB",
    "CHIT1",
    "CHRDL1",
    "CHRDL2",
    "CKAP4",
    "CLEC14A",
    "CNTN5",
    "COL15A1",
    "COL6A3",
    "COL9A1",
    "CR2",
    "CRH",
    "CRTAC1",
    "CRYBB2",
    "CSPG5",
    "CST1",
    "CST5",
    "CTHRC1",
    "CTSF",
    "CTSV",
    "CXADR",
    "CXCL12",
    "CXCL14",
    "CXCL17",
    "CXCL9",
    "CYB5R2",
    "CYTL1",
    "DCBLD2",
    "DCN",
    "DIPK2B",
    "DKK3",
    "DKKL1",
    "DLK1",
    "DMP1",
    "DPEP2",
    "DPT",
    "EDA2R",
    "EDDM3B",
    "EDIL3",
    "EFEMP1",
    "EFHD1",
    "EGFR",
    "ELN",
    "ENAH",
    "ENG",
    "ENO3",
    "ENPP2",
    "ENPP5",
    "ERBB4",
    "FABP4",
    "FAM3B",
    "FAP",
    "FAS",
    "FASLG",
    "FBLN2",
    "FCRL2",
    "FGF5",
    "FSHB",
    "FSTL1",
    "GAS6",
    "GDF15",
    "GFAP",
    "GFRAL",
    "GHRL",
    "GIP",
    "GIPC2",
    "GP2",
    "GZMB",
    "HAVCR1",
    "HMCN2",
    "HSD11B1",
    "IGDCC4",
    "IL17D",
    "IL5RA",
    "IL7R",
    "INSL3",
    "ITGAV",
    "ITGB5",
    "ITGBL1",
    "KIF22",
    "KIT",
    "KLK14",
    "KLK3",
    "KLK4",
    "KLK7",
    "KLK8",
    "KLRF1",
    "L1CAM",
    "LACRT",
    "LECT2",
    "LEG1",
    "LHB",
    "LMOD1",
    "LPO",
    "LTBP2",
    "LYPD3",
    "MAMDC4",
    "MATN3",
    "MEP1B",
    "MEPE",
    "MERTK",
    "MFGE8",
    "MLN",
    "MMP12",
    "MOG",
    "MXRA8",
    "NCAN",
    "NEFL",
    "NME3",
    "NOTCH3",
    "NPL",
    "NPTX2",
    "NTF3",
    "NTF4",
    "NTPROBNP",
    "ODAM",
    "PAEP",
    "PAMR1",
    "PINLYP",
    "PKD1",
    "PLAT",
    "PODXL2",
    "POMC",
    "PRELP",
    "PRL",
    "PRND",
    "PROK1",
    "PSPN",
    "PTGDS",
    "PTN",
    "PTPRM",
    "PTPRN2",
    "PTPRR",
    "PTPRZ1",
    "REN",
    "RET",
    "RGMA",
    "RGMB",
    "RLN2",
    "ROBO1",
    "RRM2",
    "SCARF2",
    "SCG2",
    "SCG3",
    "SCGB1A1",
    "SDK2",
    "SEPTIN3",
    "SOD2",
    "SORCS2",
    "SOST",
    "SPINK1",
    "SPON2",
    "SPRR3",
    "SRPX",
    "SUSD2",
    "SUSD5",
    "TFF1",
    "THBS2",
    "TNFRSF11B",
    "TNFRSF13B",
    "TNFSF13",
    "TNXB",
    "TSPAN1",
    "WFDC2",
    "WIF1",
    "WNT9A",
    "XCL1",
]


@dataclass(frozen=True)
class ClinicalSpec:
    field_id: str
    name: str
    source_file: str
    source_columns: tuple[str, ...]
    rule: str = "direct"
    lower: float | None = 0.0
    upper: float | None = None


def c(
    field_id: str,
    name: str,
    source_file: str,
    *source_columns: str,
    rule: str = "direct",
    lower: float | None = 0.0,
    upper: float | None = None,
) -> ClinicalSpec:
    return ClinicalSpec(field_id, name, source_file, tuple(source_columns), rule, lower, upper)


CLINICAL_SPECS = [
    c(
        "Handgrip strength (average)",
        "Handgrip strength (average)",
        "anthropometry_0.csv",
        "Hand grip strength (left)",
        "Hand grip strength (right)",
        rule="row_mean",
        lower=0.1,
        upper=100,
    ),
    c(
        "48-0.0",
        "Waist circumference",
        "anthropometry_0.csv",
        "Waist circumference",
        lower=30,
        upper=250,
    ),
    c(
        "49-0.0",
        "Hip circumference",
        "anthropometry_0.csv",
        "Hip circumference",
        lower=30,
        upper=250,
    ),
    c("50-0.0", "Standing height", "anthropometry_0.csv", "Standing height", lower=100, upper=250),
    c(
        "102-0.0",
        "Pulse rate, automated reading",
        "blood_pressure_0.csv",
        "Pulse rate",
        lower=20,
        upper=250,
    ),
    c(
        "3062-0.0",
        "Forced vital capacity (FVC)",
        "spirometry_0.csv",
        "Forced vital capacity (FVC) 1",
        "Forced vital capacity (FVC) 2",
        "Forced vital capacity (FVC) 3",
        rule="max_measured_blow",
        lower=0.5,
        upper=10,
    ),
    c(
        "3063-0.0",
        "Forced expiratory volume in 1-second (FEV1)",
        "spirometry_0.csv",
        "Forced expiratory volume in 1-second (FEV1) 1",
        "Forced expiratory volume in 1-second (FEV1) 2",
        "Forced expiratory volume in 1-second (FEV1) 3",
        rule="max_measured_blow",
        lower=0.5,
        upper=8,
    ),
    c(
        "3064-0.0",
        "Peak expiratory flow (PEF)",
        "spirometry_0.csv",
        "Peak expiratory flow (PEF) 1",
        "Peak expiratory flow (PEF) 2",
        "Peak expiratory flow (PEF) 3",
        rule="max_measured_blow",
        lower=30,
        upper=1000,
    ),
    c(
        "4079-0.0",
        "Diastolic blood pressure, automated reading",
        "blood_pressure_0.csv",
        "Diastolic blood pressure",
        lower=30,
        upper=200,
    ),
    c(
        "4080-0.0",
        "Systolic blood pressure, automated reading",
        "blood_pressure_0.csv",
        "Systolic blood pressure",
        lower=60,
        upper=300,
    ),
    c(
        "21001-0.0",
        "Body mass index (BMI)",
        "anthropometry_0.csv",
        "Body mass index (BMI)",
        lower=10,
        upper=80,
    ),
    c("21002-0.0", "Weight", "anthropometry_0.csv", "Weight", lower=20, upper=300),
]

for field_id, name in [
    ("30000-0.0", "White blood cell (leukocyte) count"),
    ("30010-0.0", "Red blood cell (erythrocyte) count"),
    ("30020-0.0", "Haemoglobin concentration"),
    ("30030-0.0", "Haematocrit percentage"),
    ("30040-0.0", "Mean corpuscular volume"),
    ("30050-0.0", "Mean corpuscular haemoglobin"),
    ("30060-0.0", "Mean corpuscular haemoglobin concentration"),
    ("30070-0.0", "Red blood cell (erythrocyte) distribution width"),
    ("30080-0.0", "Platelet count"),
    ("30090-0.0", "Platelet crit"),
    ("30100-0.0", "Mean platelet (thrombocyte) volume"),
    ("30110-0.0", "Platelet distribution width"),
    ("30120-0.0", "Lymphocyte count"),
    ("30130-0.0", "Monocyte count"),
    ("30140-0.0", "Neutrophill count"),
    ("30150-0.0", "Eosinophill count"),
    ("30160-0.0", "Basophill count"),
    ("30170-0.0", "Nucleated red blood cell count"),
    ("30180-0.0", "Lymphocyte percentage"),
    ("30190-0.0", "Monocyte percentage"),
    ("30200-0.0", "Neutrophill percentage"),
    ("30210-0.0", "Eosinophill percentage"),
    ("30220-0.0", "Basophill percentage"),
    ("30230-0.0", "Nucleated red blood cell percentage"),
    ("30240-0.0", "Reticulocyte percentage"),
    ("30250-0.0", "Reticulocyte count"),
    ("30260-0.0", "Mean reticulocyte volume"),
    ("30270-0.0", "Mean sphered cell volume"),
    ("30280-0.0", "Immature reticulocyte fraction"),
    ("30290-0.0", "High light scatter reticulocyte percentage"),
    ("30300-0.0", "High light scatter reticulocyte count"),
]:
    CLINICAL_SPECS.append(c(field_id, name, "blood_count_0.csv", name))

for field_id, output_name, source_name in [
    ("30510-0.0", "Creatinine (enzymatic) in urine", "Creatinine (enzymatic) in urine"),
    ("30520-0.0", "Potassium in urine", "Potassium in urine"),
    ("30530-0.0", "Sodium in urine", "Sodium in urine"),
]:
    CLINICAL_SPECS.append(c(field_id, output_name, "urine_0.csv", source_name))

for field_id, output_name, source_name in [
    ("30600-0.0", "Albumin", "Albumin"),
    ("30610-0.0", "Alkaline phosphatase", "Alkaline phosphatase"),
    ("30620-0.0", "Alanine aminotransferase (ALT)", "Alanine aminotransferase"),
    ("30630-0.0", "Apolipoprotein A", "Apolipoprotein A"),
    ("30640-0.0", "Apolipoprotein B", "Apolipoprotein B"),
    ("30650-0.0", "Aspartate aminotransferase (AST)", "Aspartate aminotransferase"),
    ("30660-0.0", "Direct bilirubin", "Direct bilirubin"),
    ("30670-0.0", "Urea", "Urea"),
    ("30680-0.0", "Calcium", "Calcium"),
    ("30690-0.0", "Cholesterol", "Cholesterol"),
    ("30700-0.0", "Creatinine", "Creatinine"),
    ("30710-0.0", "C-reactive protein", "C-reactive protein"),
    ("30720-0.0", "Cystatin C", "Cystatin C"),
    ("30730-0.0", "Gamma glutamyltransferase", "Gamma glutamyltransferase"),
    ("30740-0.0", "Glucose", "Glucose"),
    ("30750-0.0", "Glycated haemoglobin (HbA1c)", "Glycated haemoglobin (HbA1c)"),
    ("30760-0.0", "HDL cholesterol", "HDL cholesterol"),
    ("30770-0.0", "IGF-1", "IGF-1"),
    ("30780-0.0", "LDL direct", "LDL direct"),
    ("30790-0.0", "Lipoprotein A", "Lipoprotein A"),
    ("30810-0.0", "Phosphate", "Phosphate"),
    ("30830-0.0", "SHBG", "SHBG"),
    ("30840-0.0", "Total bilirubin", "Total bilirubin"),
    ("30850-0.0", "Testosterone", "Testosterone"),
    ("30860-0.0", "Total protein", "Total protein"),
    ("30870-0.0", "Triglycerides", "Triglycerides"),
    ("30880-0.0", "Urate", "Urate"),
    ("30890-0.0", "Vitamin D", "Vitamin D"),
]:
    CLINICAL_SPECS.append(c(field_id, output_name, "blood_chemistry_0.csv", source_name))

CLINICAL_SPECS.extend(
    [
        c(
            "FEV1-FVC ratio",
            "FEV1-FVC ratio",
            "spirometry_0.csv",
            "Forced expiratory volume in 1-second (FEV1)",
            "Forced vital capacity (FVC)",
            rule="derived_from_isolated_measured_blows",
            lower=0.2,
            upper=1.2,
        ),
        c(
            "Waist-hip circumference ratio",
            "Waist-hip circumference ratio",
            "anthropometry_0.csv",
            "Waist circumference",
            "Hip circumference",
            rule="ratio",
            lower=0.3,
            upper=2.0,
        ),
        c(
            "Heel bone mineral density (average)",
            "Heel bone mineral density (average)",
            "processed_heel_instance_0.csv",
            "Heel bone mineral density (BMD)",
            lower=-2.0,
            upper=3.0,
        ),
        c(
            "Ankle spacing width (average)",
            "Ankle spacing width (average)",
            "processed_heel_instance_0.csv",
            "Ankle spacing width",
            lower=10,
            upper=300,
        ),
    ]
)


def csv_columns(path: Path) -> list[str]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return next(csv.reader(handle))


def ensure_unique_eid(df: pd.DataFrame, context: str) -> None:
    if df["eid"].isna().any():
        raise ValueError(f"{context}: missing eid")
    if df["eid"].duplicated().any():
        raise ValueError(f"{context}: duplicated eid={int(df['eid'].duplicated().sum())}")


def load_population() -> tuple[pd.DataFrame, dict[str, int | float]]:
    usecols = [
        "eid",
        "Date of birth",
        "Date of attending assessment centre",
        "UK Biobank assessment centre",
        "Age when attended assessment centre",
        "Sex",
    ]
    raw = pd.read_csv(POPULATION_PATH, usecols=usecols)
    raw["eid"] = pd.to_numeric(raw["eid"], errors="raise").astype("int64")
    ensure_unique_eid(raw, "population")
    raw["Date_birth"] = pd.to_datetime(raw["Date of birth"], errors="coerce")
    raw["Date_attend"] = pd.to_datetime(raw["Date of attending assessment centre"], errors="coerce")
    raw["Age_unrounded"] = (raw["Date_attend"] - raw["Date_birth"]).dt.total_seconds() / (
        365.2425 * 86400.0
    )
    raw["Sex"] = raw["Sex"].map({0: "female", 1: "male"})
    raw["AssessmentCenter"] = pd.to_numeric(raw["UK Biobank assessment centre"], errors="coerce")
    valid_metadata = raw[
        raw[["Date_birth", "Date_attend", "Age_unrounded", "Sex", "AssessmentCenter"]]
        .notna()
        .all(axis=1)
    ].copy()
    age_eligible = valid_metadata[
        valid_metadata["Age_unrounded"].between(AGE_MIN, AGE_MAX, inclusive="both")
    ].copy()
    age_eligible["Age"] = age_eligible["Age_unrounded"].round(1).astype("float32")
    age_eligible["Sex_male"] = (age_eligible["Sex"] == "male").astype("int8")
    age_eligible["AssessmentCenter"] = age_eligible["AssessmentCenter"].astype("int32")
    age_eligible["Date_birth"] = age_eligible["Date_birth"].dt.strftime("%Y-%m-%d")
    age_eligible["Date_attend"] = age_eligible["Date_attend"].dt.strftime("%Y-%m-%d")
    out = age_eligible[
        ["eid", "Age", "Sex", "Sex_male", "AssessmentCenter", "Date_birth", "Date_attend"]
    ].reset_index(drop=True)
    stats = {
        "population_source_rows": int(len(raw)),
        "complete_core_metadata_rows": int(len(valid_metadata)),
        "age_40_70_rows": int(len(out)),
        "age_original_disagreement_n": int(
            (
                age_eligible["Age"].round().astype(int)
                != age_eligible["Age when attended assessment centre"]
            ).sum()
        ),
        "age_difference_vs_original_mean": float(
            (age_eligible["Age"] - age_eligible["Age when attended assessment centre"]).mean()
        ),
    }
    return out, stats


def load_history() -> pd.DataFrame:
    history = pd.read_feather(HISTORY_PATH, columns=["eid", *HISTORY_COLUMNS])
    history["eid"] = pd.to_numeric(history["eid"], errors="raise").astype("int64")
    ensure_unique_eid(history, "disease history")
    invalid = ~history[HISTORY_COLUMNS].isin([0, 1])
    if invalid.any().any():
        raise ValueError("Disease-history table contains values outside {0,1}")
    history["Baseline_unhealthy"] = history[HISTORY_COLUMNS].max(axis=1).astype("int8")
    return history[["eid", "Baseline_unhealthy"]]


def read_numeric_csv(path: Path, columns: list[str]) -> pd.DataFrame:
    df = pd.read_csv(path, usecols=columns, low_memory=False)
    df["eid"] = pd.to_numeric(df["eid"], errors="raise").astype("int64")
    ensure_unique_eid(df, path.name)
    for column in columns:
        if column != "eid":
            df[column] = pd.to_numeric(df[column], errors="coerce").astype("float32")
    return df


def apply_bounds(
    series: pd.Series, lower: float | None, upper: float | None
) -> tuple[pd.Series, int]:
    out = pd.to_numeric(series, errors="coerce").replace([np.inf, -np.inf], np.nan)
    invalid = pd.Series(False, index=out.index)
    if lower is not None:
        invalid |= out < lower
    if upper is not None:
        invalid |= out > upper
    count = int((invalid & out.notna()).sum())
    return out.mask(invalid).astype("float32"), count


def source_path_for_clinical(name: str) -> Path:
    mapping = {
        "anthropometry_0.csv": PROCESSED
        / "02-AssessmentCenter/05-PhysicalMeasures/anthropometry_0.csv",
        "blood_pressure_0.csv": PROCESSED
        / "02-AssessmentCenter/05-PhysicalMeasures/blood_pressure_0.csv",
        "spirometry_0.csv": PROCESSED / "02-AssessmentCenter/05-PhysicalMeasures/spirometry_0.csv",
        "processed_heel_instance_0.csv": PROCESSED
        / "02-AssessmentCenter/05-PhysicalMeasures/processed_heel_instance_0.csv",
        "blood_count_0.csv": PROCESSED / "03-BiologicalSamples/blood_count_0.csv",
        "blood_chemistry_0.csv": PROCESSED / "03-BiologicalSamples/blood_chemistry_0.csv",
        "urine_0.csv": PROCESSED / "03-BiologicalSamples/urine_0.csv",
    }
    return mapping[name]


def attach_metadata_and_filter(
    modality: str,
    population: pd.DataFrame,
    raw_features: pd.DataFrame,
    features: list[str],
    history: pd.DataFrame,
    source_rows_with_any_measurement: int,
) -> tuple[pd.DataFrame, pd.DataFrame, list[str], dict[str, int | float]]:
    merged = population.merge(raw_features, on="eid", how="left", validate="one_to_one")
    before_row_missing = len(merged)
    age_eligible_with_any_measurement = int(merged[features].notna().any(axis=1).sum())
    row_missing = merged[features].isna().mean(axis=1)
    row_keep = row_missing < MISSING_THRESHOLD
    retained = merged.loc[row_keep].copy()
    feature_missing = retained[features].isna().mean()
    feature_keep = feature_missing < MISSING_THRESHOLD
    exception_caps = FEATURE_MISSINGNESS_EXCEPTIONS.get(modality, {})
    for feature, cap in exception_caps.items():
        feature_keep.loc[feature] = bool(feature_missing.loc[feature] <= cap)
    final_features = feature_missing.index[feature_keep].tolist()
    dropped_features = feature_missing.index[~feature_keep].tolist()
    exception_applied = {
        feature: bool(
            feature in final_features
            and feature_missing.loc[feature] >= MISSING_THRESHOLD
            and feature_missing.loc[feature] <= cap
        )
        for feature, cap in exception_caps.items()
    }
    retained = retained[[*population.columns, *final_features]]
    before_history = len(retained)
    retained = retained.merge(history, on="eid", how="inner", validate="one_to_one")
    retained = (
        retained[[*METADATA_COLUMNS, *final_features]].sort_values("eid").reset_index(drop=True)
    )
    retained["Baseline_unhealthy"] = retained["Baseline_unhealthy"].astype("int8")
    stats = {
        "modality": modality,
        "age_eligible_population": int(before_row_missing),
        "source_rows_with_any_measurement": source_rows_with_any_measurement,
        "age_eligible_with_any_measurement": age_eligible_with_any_measurement,
        "excluded_participant_missing_ge_20pct": int((~row_keep).sum()),
        "retained_after_participant_missingness": int(len(row_keep[row_keep])),
        "excluded_missing_health_metadata": int(before_history - len(retained)),
        "final_participants": int(len(retained)),
        "candidate_features": int(len(features)),
        "excluded_feature_missing_ge_20pct": int(len(dropped_features)),
        "features_retained_by_named_exception": int(sum(exception_applied.values())),
        "final_features": int(len(final_features)),
        "healthy_n": int((retained["Baseline_unhealthy"] == 0).sum()),
        "unhealthy_n": int((retained["Baseline_unhealthy"] == 1).sum()),
        "female_n": int((retained["Sex"] == "female").sum()),
        "male_n": int((retained["Sex"] == "male").sum()),
        "age_min": float(retained["Age"].min()),
        "age_max": float(retained["Age"].max()),
        "max_final_participant_missing_rate": float(
            retained[final_features].isna().mean(axis=1).max()
        ),
        "max_final_feature_missing_rate": float(retained[final_features].isna().mean().max()),
    }
    profile = pd.DataFrame(
        {
            "modality": modality,
            "feature": features,
            "missing_rate_after_participant_filter": feature_missing.reindex(features).to_numpy(),
            "missingness_rule": [
                "named_exception_le_20.05pct" if feature in exception_caps else "standard_lt_20pct"
                for feature in features
            ],
            "exception_cap": [exception_caps.get(feature, np.nan) for feature in features],
            "exception_applied": [
                int(exception_applied.get(feature, False)) for feature in features
            ],
            "include_final": [int(feature in final_features) for feature in features],
        }
    )
    return retained, profile, final_features, stats


def build_nmr(population: pd.DataFrame, history: pd.DataFrame):
    headers = csv_columns(NMR_PATH)
    features = [
        column
        for column in headers[1:]
        if "percentage" not in column.lower() and " ratio" not in column.lower()
    ]
    if len(features) != 168:
        raise ValueError(f"NMR rule must yield 168 features, observed {len(features)}")
    raw = read_numeric_csv(NMR_PATH, ["eid", *features])
    any_measurement = int(raw[features].notna().any(axis=1).sum())
    return attach_metadata_and_filter("nmr", population, raw, features, history, any_measurement)


def build_protein(population: pd.DataFrame, history: pd.DataFrame):
    if len(PROTEIN_FEATURES) != 204 or len(set(PROTEIN_FEATURES)) != 204:
        raise ValueError("Expected 204 unique proteins in the feature list")
    raw = read_numeric_csv(PROTEIN_PATH, ["eid", *PROTEIN_FEATURES])
    any_measurement = int(raw[PROTEIN_FEATURES].notna().any(axis=1).sum())
    return attach_metadata_and_filter(
        "protein", population, raw, PROTEIN_FEATURES, history, any_measurement
    )


def build_clinical(population: pd.DataFrame, history: pd.DataFrame):
    features = [spec.name for spec in CLINICAL_SPECS]
    if len(features) != 78 or len(set(features)) != 78:
        raise ValueError("Expected 78 unique clinical features")
    grouped: dict[str, set[str]] = {}
    for spec in CLINICAL_SPECS:
        if spec.rule == "derived_from_isolated_measured_blows":
            continue
        grouped.setdefault(spec.source_file, set()).update(spec.source_columns)
    source_frames: dict[str, pd.DataFrame] = {}
    for source_file, columns in grouped.items():
        source_frames[source_file] = read_numeric_csv(
            source_path_for_clinical(source_file), ["eid", *sorted(columns)]
        )

    clinical = population[["eid"]].copy()
    invalidated: dict[str, int] = {}
    calculated: dict[str, pd.Series] = {}
    for spec in CLINICAL_SPECS:
        if spec.rule == "derived_from_isolated_measured_blows":
            continue
        frame = source_frames[spec.source_file].set_index("eid")
        if spec.rule == "direct":
            values = frame[spec.source_columns[0]]
        elif spec.rule == "row_mean":
            values = frame[list(spec.source_columns)].mean(axis=1, skipna=True)
        elif spec.rule == "max_measured_blow":
            bounded = []
            invalid_blow_n = 0
            for column in spec.source_columns:
                value, invalid_n = apply_bounds(frame[column], spec.lower, spec.upper)
                invalid_blow_n += invalid_n
                bounded.append(value)
            values = pd.concat(bounded, axis=1).max(axis=1, skipna=True)
        elif spec.rule == "ratio":
            values = frame[spec.source_columns[0]] / frame[spec.source_columns[1]]
        else:
            raise ValueError(f"Unsupported clinical rule: {spec.rule}")
        bounded_values, invalid_n = apply_bounds(values, spec.lower, spec.upper)
        calculated[spec.name] = bounded_values
        invalidated[spec.name] = invalid_n + (
            invalid_blow_n if spec.rule == "max_measured_blow" else 0
        )

    fev1 = calculated["Forced expiratory volume in 1-second (FEV1)"]
    fvc = calculated["Forced vital capacity (FVC)"]
    ratio_spec = next(spec for spec in CLINICAL_SPECS if spec.name == "FEV1-FVC ratio")
    ratio, ratio_invalid = apply_bounds(fev1 / fvc, ratio_spec.lower, ratio_spec.upper)
    calculated[ratio_spec.name] = ratio
    invalidated[ratio_spec.name] = ratio_invalid

    for name, values in calculated.items():
        clinical = clinical.merge(
            values.rename(name).reset_index(), on="eid", how="left", validate="one_to_one"
        )
    any_measurement = int(clinical[features].notna().any(axis=1).sum())
    cohort, profile, final_features, stats = attach_metadata_and_filter(
        "clinical", population, clinical, features, history, any_measurement
    )
    spec_frame = pd.DataFrame(
        [
            {
                "field_id": spec.field_id,
                "feature": spec.name,
                "source_file": str(source_path_for_clinical(spec.source_file)),
                "source_columns": " | ".join(spec.source_columns),
                "construction_rule": spec.rule,
                "plausibility_lower": spec.lower,
                "plausibility_upper": spec.upper,
                "invalid_values_set_missing": invalidated.get(spec.name, 0),
            }
            for spec in CLINICAL_SPECS
        ]
    )
    profile = profile.merge(spec_frame, on="feature", how="left", validate="one_to_one")
    return cohort, profile, final_features, stats


def atomic_feather(df: pd.DataFrame, path: Path) -> None:
    temp = path.with_suffix(path.suffix + ".tmp")
    df.to_feather(temp, compression="zstd")
    os.replace(temp, path)


def write_cohort(
    modality: str,
    cohort: pd.DataFrame,
    profile: pd.DataFrame,
    features: list[str],
) -> None:
    output_name = "routine_blood" if modality == "clinical" else modality
    data_path = DATA_DIR / f"{output_name}.feather"
    feature_path = DATA_DIR / f"{output_name}_features.txt"
    atomic_feather(cohort, data_path)
    feature_path.write_text("\n".join(features) + "\n", encoding="utf-8")
    profile.to_csv(REPORT_DIR / f"{output_name}_features.tsv", sep="\t", index=False)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    output_names = ("nmr.feather", "protein.feather", "routine_blood.feather")
    existing = [DATA_DIR / name for name in output_names if (DATA_DIR / name).exists()]
    if existing and not args.overwrite:
        raise FileExistsError(f"Modeling data already exist: {existing}")

    population, _ = load_population()
    history = load_history()
    builders = {"nmr": build_nmr, "protein": build_protein, "clinical": build_clinical}
    cohort_flow = []
    for modality, builder in builders.items():
        cohort, profile, features, statistics = builder(population, history)
        write_cohort(modality, cohort, profile, features)
        print(f"{modality}: {len(cohort):,} participants, {len(features)} features")
        cohort_flow.append(statistics)
    pd.DataFrame(cohort_flow).to_csv(REPORT_DIR / "cohort_flow.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
