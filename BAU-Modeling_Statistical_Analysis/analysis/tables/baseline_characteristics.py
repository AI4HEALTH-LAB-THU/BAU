"""Summarize baseline characteristics in UKB and the external cohorts."""

import os
from pathlib import Path

import pandas as pd

from bau.paths import DATA_ROOT, RESULTS_ROOT

UKB = Path(os.environ.get("UKB_DATA_ROOT", "/mnt/project/UKB"))
OUTPUT = RESULTS_ROOT / "tables"


def normalise_id(x):
    return x.astype(str).str.strip().str.replace(r"\.0$", "", regex=True)


def summarize(frame, continuous, indicators):
    rows = []
    for name, values in continuous.items():
        x = pd.to_numeric(values, errors="coerce").dropna()
        rows.append(
            dict(
                variable=name,
                population_n=len(frame),
                observed_n=len(x),
                mean=x.mean(),
                sd=x.std(ddof=1),
                median=x.median(),
                q1=x.quantile(0.25),
                q3=x.quantile(0.75),
            )
        )
    for name, indicator in indicators.items():
        count = int(indicator.fillna(False).sum())
        rows.append(
            dict(
                variable=name,
                population_n=len(frame),
                count=count,
                percent=100 * count / len(frame),
            )
        )
    return pd.DataFrame(rows)


def ukb_summary(x):
    ethnicity = x["Ethnic background"]
    continuous = {
        "Age": x.Age,
        "BMI": x.BMI,
        "Systolic blood pressure": x["Systolic blood pressure"],
        "Daily moderate to vigorous physical activity": x[
            "Daily moderate to vigorous physical activity"
        ],
    }
    indicators = {
        "Female": x.Sex.eq("female"),
        "White": ethnicity.isin([1, 1001, 1002, 1003]),
        "Mixed": ethnicity.isin([2, 2001, 2002, 2003, 2004]),
        "Asian": ethnicity.isin([3, 3001, 3002, 3003, 3004]),
        "Black": ethnicity.isin([4, 4001, 4002, 4003]),
        "Chinese": ethnicity.eq(5),
        "Baseline healthy": x.Baseline_unhealthy.eq(0),
    }
    indicators["Other or unknown ethnicity"] = ~pd.concat(
        [indicators[k] for k in ("White", "Mixed", "Asian", "Black", "Chinese")], axis=1
    ).any(axis=1)
    for name in (
        "Current smoker",
        "Daily alcohol intake",
        "Hypertension history",
        "Diabetes history",
    ):
        indicators[name] = pd.to_numeric(x[name], errors="coerce").eq(1)
    indicators["College education"] = x.Education.astype(str).eq("college")
    indicators["Paid employment"] = x.Employment.astype(str).eq("in paid employment")
    for level in ("low", "middle", "high", "unknown"):
        indicators[f"Income {level}"] = x.Income.astype(str).eq(level)
    return summarize(x, continuous, indicators)


def read_external(cohort):
    source_dir = DATA_ROOT / "external" / cohort
    if cohort == "CHARLS":
        x = pd.read_feather(source_dir / "candidate_qc.feather")
        survival = pd.read_feather(source_dir / "survival_outcomes.feather")
        fields = ["participant_id", "survival_eligible", "time_midpoint_years", "event"]
        survival = survival[fields].copy()
        for frame in (x, survival):
            frame["participant_id"] = normalise_id(frame.participant_id)
        x = x.drop(columns=[c for c in fields[1:] if c in x]).merge(
            survival, on="participant_id", validate="one_to_one"
        )
    else:
        x = pd.read_feather(source_dir / "analysis_wide.feather")
    x["participant_id"] = normalise_id(x.participant_id)
    prediction = pd.read_feather(
        RESULTS_ROOT / "models/external" / cohort.lower() / "predictions.feather"
    )
    ids = pd.DataFrame({"participant_id": normalise_id(prediction.eid)})
    x = ids.merge(x, on="participant_id", validate="one_to_one")
    return x.loc[
        x.age.between(40, 70) & x.survival_eligible.eq(True) & x.time_midpoint_years.gt(0)
    ].copy()


def external_summary(x, cohort):
    bmi = pd.to_numeric(x.bmi, errors="coerce")
    sbp = pd.to_numeric(x.systolic_bp, errors="coerce")
    continuous = {
        "Age": x.age,
        "BMI": bmi.where(bmi.between(10, 100)),
        "Systolic blood pressure": sbp.where(sbp.between(60, 300)),
    }
    indicators = {
        "Female": x.sex.eq("female"),
        "Hypertension history": x.hypertension_status.eq("yes"),
        "Diabetes history": x.diabetes_status.eq("yes"),
        "Baseline disease": x.baseline_health_status.eq("disease"),
    }
    if cohort == "NHANES":
        for name, codes in {
            "Mexican American": [1],
            "Other Hispanic": [2],
            "White": [3],
            "Black": [4],
            "Asian": [6],
            "Other ethnicity": [5, 7],
        }.items():
            indicators[name] = x.race_ethnicity_raw.isin(codes)
        indicators["Ever smoked"] = x.smoking_status.isin(["former", "current"])
        indicators["Secondary education or higher"] = x.education_raw.isin([3, 4, 5])
        indicators["US citizen"] = x.citizenship_raw.eq(1)
        income = x.income_poverty_ratio
        indicators.update(
            {
                "Income low": income.lt(1.3),
                "Income middle": income.ge(1.3) & income.lt(3.5),
                "Income high": income.ge(3.5),
                "Income unknown": income.isna(),
            }
        )
    elif cohort == "CHARLS":
        indicators["Ever smoked"] = x.smoking_status.eq("ever")
        indicators["Secondary education or higher"] = x.education_raw.astype("string").str.match(
            r"^(6|7|8|9|10) ", na=False
        )
    else:
        indicators["Ever smoked"] = x.smoking_status.eq("ever")
        indicators["Secondary education or higher"] = x.education_raw.ge(3)
        indicators["Han ethnicity"] = x.ethnicity_raw.eq(1)
        indicators["Minority ethnicity"] = x.ethnicity_raw.gt(0) & x.ethnicity_raw.ne(1)
        indicators["Unknown ethnicity"] = x.ethnicity_raw.isna() | x.ethnicity_raw.le(0)
        for name, field in {
            "Urban": "urban_rural_raw",
            "Working": "working_raw",
            "Medical insurance": "medical_insurance_raw",
        }.items():
            indicators[name] = x[field].eq(1)
    return summarize(x, continuous, indicators)


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    support = DATA_ROOT / "ukb/support"
    participants = pd.read_feather(support / "participants.feather")
    covariates = pd.read_feather(support / "covariates.feather")
    covariates = covariates.drop(
        columns=[c for c in participants if c != "eid" and c in covariates]
    )
    ethnicity = pd.read_csv(
        UKB / "ProcessedData/02-AssessmentCenter/02-Touchscreen/sociodemographics_0.csv",
        usecols=["eid", "Ethnic background"],
    )
    data = participants.merge(covariates, on="eid", validate="one_to_one").merge(
        ethnicity, on="eid", how="left", validate="one_to_one"
    )
    cohorts = {"NMR": "in_nmr", "Proteomic": "in_protein", "Routine blood": "in_clinical"}
    table1 = []
    for name, flag in cohorts.items():
        x = data.loc[data[flag].eq(1)]
        table1.append(ukb_summary(x).assign(cohort=name))
    pd.concat(table1).to_csv(OUTPUT / "ukb_baseline.tsv", sep="\t", index=False)

    organ = pd.read_feather(
        RESULTS_ROOT / "phenotypes/organ_disagreement.feather", columns=["eid"]
    ).eid
    repeat = pd.read_feather(RESULTS_ROOT / "phenotypes/repeat_pairs.feather", columns=["eid"]).eid
    table2 = []
    for name, ids in (("Organ-age disagreement", organ), ("Repeat NMR", repeat)):
        x = data.loc[data.eid.isin(ids)]
        if len(x) != len(ids) or ids.duplicated().any():
            raise ValueError(f"Incomplete or duplicated {name} participant identifiers")
        table2.append(ukb_summary(x).assign(cohort=name))
    pd.concat(table2).to_csv(OUTPUT / "phenotype_baseline.tsv", sep="\t", index=False)

    table3, followup = [], []
    for cohort in ("NHANES", "CHARLS", "CHNS"):
        x = read_external(cohort)
        table3.append(external_summary(x, cohort).assign(cohort=cohort))
        followup.append(
            dict(
                cohort=cohort,
                n_risk=len(x),
                n_events=int(x.event.sum()),
                follow_up_q1_years=x.time_midpoint_years.quantile(0.25),
                median_follow_up_years=x.time_midpoint_years.median(),
                follow_up_q3_years=x.time_midpoint_years.quantile(0.75),
            )
        )
    pd.concat(table3).to_csv(OUTPUT / "external_baseline.tsv", sep="\t", index=False)
    pd.DataFrame(followup).to_csv(OUTPUT / "external_follow_up.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
