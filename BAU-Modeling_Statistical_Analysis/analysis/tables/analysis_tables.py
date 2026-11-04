"""Export analysis populations and manuscript variable definitions."""

import argparse
import json

import pandas as pd

from bau.paths import CODE_ROOT as ROOT
from bau.paths import RESULTS_ROOT

OUTPUT = RESULTS_ROOT / "tables"


def write(frame, name):
    frame.to_csv(OUTPUT / name, sep="\t", index=False)


def definition_tables():
    """Join manuscript definitions to the feature and analysis configuration."""
    path = ROOT / "config/table_definitions.json"
    definitions = json.loads(path.read_text())
    tables = {
        name: pd.DataFrame(table["records"], columns=table["columns"])
        for name, table in definitions.items()
    }

    features = tables["age_model_features"]
    features["modality"] = (
        features.section.str.split()
        .str[0]
        .map(
            {
                "Predictor": "shared",
                "NMR": "nmr",
                "Plasma": "protein",
                "Haematology": "blood",
                "Serum": "blood",
            }
        )
    )
    model_features = [{"modality": "shared", "order": 1, "feature": "Sex"}]
    for modality in ("nmr", "protein", "blood"):
        path = ROOT / f"config/features/{modality}.txt"
        names = [line for line in path.read_text().splitlines() if line.strip()]
        model_features.extend(
            {"modality": modality, "order": order, "feature": feature}
            for order, feature in enumerate(names, 1)
        )
    models = pd.DataFrame(model_features)
    models["Model feature"] = models.feature.replace(
        {
            "Remnant Cholesterol (Non-HDL, Non-LDL -Cholesterol)": "Remnant Cholesterol (Non-HDL, Non-LDL Cholesterol)",
            **{
                f"{cell}ll {measure}": f"{cell}l {measure}"
                for cell in ("Neutrophi", "Eosinophi", "Basophi")
                for measure in ("count", "percentage")
            },
        }
    )
    tables["age_model_features"] = models.merge(
        features, on=["modality", "Model feature"], validate="one_to_one"
    )
    if len(tables["age_model_features"]) != len(models):
        raise ValueError("Age-model feature definitions are incomplete")

    organ_rows = []
    config_dir = ROOT / "config/organ_features"
    for path in sorted(config_dir.glob("*.txt")):
        proteins = [line for line in path.read_text().splitlines() if line.strip()]
        organ_rows.append(
            {
                "Organ-age clock": path.stem.title(),
                "Proteins, n": len(proteins),
                "Gene symbols": "; ".join(proteins),
            }
        )
    tables["organ_clock_features"] = pd.DataFrame(organ_rows)

    transfer = json.loads((ROOT / "config/external.json").read_text())
    transfer_rows = [
        dict(cohort=cohort.upper(), **feature)
        for cohort, settings in transfer["cohorts"].items()
        for feature in settings["features"]
    ]
    harmonisation = tables["external_feature_harmonisation"].rename(
        columns={"UK Biobank feature name": "ukb_feature"}
    )
    harmonisation["cohort"] = harmonisation.section.str.split().str[0]
    tables["external_feature_harmonisation"] = pd.DataFrame(transfer_rows).merge(
        harmonisation, on=["cohort", "ukb_feature"], validate="one_to_one"
    )
    if len(tables["external_feature_harmonisation"]) != len(transfer_rows):
        raise ValueError("External feature definitions are incomplete")

    physiology = pd.read_csv(ROOT / "config/physiology.tsv", sep="\t")
    tables["physiological_traits"] = physiology.merge(
        tables["physiological_traits"], left_on="label", right_on="Trait", validate="one_to_one"
    )
    if len(tables["physiological_traits"]) != len(physiology):
        raise ValueError("Physiological trait definitions are incomplete")
    return tables


def analysis_populations():
    columns = [
        "n_risk",
        "n_events",
        "follow_up_q1_years",
        "median_follow_up_years",
        "follow_up_q3_years",
    ]
    baseline = pd.read_csv(RESULTS_ROOT / "outcomes/ngboost/joint_cox_estimates.tsv", sep="\t")
    baseline = baseline.loc[baseline.exposure.eq("BAG"), ["modality", "outcome", *columns]]
    repeat = pd.read_csv(RESULTS_ROOT / "outcomes/longitudinal/landmark_risk_sets.tsv", sep="\t")
    repeat["modality"] = "Repeat NMR"
    external = pd.read_csv(OUTPUT / "external_follow_up.tsv", sep="\t").rename(
        columns={"cohort": "modality"}
    )
    external["outcome"] = "all-cause mortality"
    return pd.concat([baseline, repeat, external], ignore_index=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--definitions-only",
        action="store_true",
        help="Export feature and analysis definitions without cohort results",
    )
    args = parser.parse_args()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    for name, frame in definition_tables().items():
        write(frame, f"{name}.tsv")
    definitions = json.loads((ROOT / "config/table_definitions.json").read_text())
    write(
        pd.DataFrame(
            [
                {"definition": name, "notes": definition["notes"]}
                for name, definition in definitions.items()
            ]
        ),
        "table_definition_notes.tsv",
    )
    if not args.definitions_only:
        write(analysis_populations(), "analysis_populations.tsv")
        timing = pd.read_csv(RESULTS_ROOT / "outcomes/external/mortality_estimates.tsv", sep="\t")
        write(timing.loc[timing.cohort.isin(["CHARLS", "CHNS"])], "mortality_timing.tsv")


if __name__ == "__main__":
    main()
