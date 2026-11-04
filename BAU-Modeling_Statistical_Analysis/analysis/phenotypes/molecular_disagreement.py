"""Score observed biomarker age support and select the four protein examples."""

import numpy as np
import pandas as pd

from bau.paths import RESULTS_ROOT

AGE_GRID = np.arange(40, 70.001, 0.25)


def likelihood(values, means, scales):
    return -0.5 * ((values[:, :, None] - means.T[None, :, :]) / scales[None, :, None]) ** 2


def disagreement(log_likelihood):
    observed = np.isfinite(log_likelihood[:, :, 0])
    count = observed.sum(axis=1)
    if (count == 0).any():
        raise ValueError("A molecular profile has no observed biomarkers")
    masked = np.where(observed[:, :, None], log_likelihood, -np.inf)
    maxima = masked.max(axis=2)
    maxima[~observed] = 0
    joint = np.where(observed[:, :, None], log_likelihood, 0).sum(axis=1) / count[:, None]
    peak = joint.argmax(axis=1)
    loss = maxima.sum(axis=1) / count - joint[np.arange(len(joint)), peak]
    return loss, peak


def read_model(directory, fold, sex):
    parameters = pd.read_csv(directory / f"parameters_{fold}_{sex}.tsv", sep="\t")
    means = pd.read_feather(directory / f"mean_grid_{fold}_{sex}.feather")
    return means.iloc[:, 1:].to_numpy(), parameters.sd.to_numpy()


def score_modality(modality):
    directory = RESULTS_ROOT / "phenotypes/molecular" / modality
    population = pd.read_feather(directory / "population.feather")
    measurements = pd.read_feather(directory / "measurements.feather")
    if not np.array_equal(population.eid, measurements.eid):
        raise ValueError("Molecular measurements and population IDs differ")
    values = measurements.iloc[:, 1:].to_numpy(float)
    scores = np.full(len(population), np.nan)
    common_age = np.full(len(population), np.nan)
    for fold in range(1, 6):
        for sex in (0, 1):
            means, scales = read_model(directory, fold, sex)
            indices = np.flatnonzero((population.OuterFold == fold) & (population.Sex_male == sex))
            for start in range(0, len(indices), 128):
                rows = indices[start : start + 128]
                scores[rows], peaks = disagreement(likelihood(values[rows], means, scales))
                common_age[rows] = AGE_GRID[peaks]
    if not np.isfinite(scores).all() or (scores < -1e-10).any():
        raise ValueError("Invalid molecular disagreement scores")
    population["H"] = scores
    population["common_age"] = common_age
    population.to_feather(directory / "scores.feather")
    if modality == "protein":
        select_examples(population, measurements, directory)


def select_examples(population, measurements, directory):
    features = measurements.columns[1:].to_numpy()
    values = measurements.iloc[:, 1:].to_numpy(float)
    eligible = population.loc[population.n_observed == len(features)].drop(columns="BAU")
    q05, q10, q90, q95 = eligible.H.quantile([0.05, 0.10, 0.90, 0.95])
    selected = []
    for sex in (0, 1):
        low = eligible[(eligible.Sex_male == sex) & eligible.H.between(q05, q10)]
        high = eligible[(eligible.Sex_male == sex) & eligible.H.between(q90, q95)]
        pairs = []
        for row in low.itertuples():
            matches = high[
                (high.OuterFold == row.OuterFold)
                & ((high.Age - row.Age).abs() <= 1)
                & ((high.BAG - row.BAG).abs() <= 0.15)
            ]
            cost = (
                (matches.Age - row.Age) ** 2
                + ((matches.BAG - row.BAG) / 0.15) ** 2
                + (((matches.Age + row.Age) / 2 - 55) / 5) ** 2
                + (((matches.BAG + row.BAG) / 2) / 0.5) ** 2
            )
            pairs.extend(
                (float(value), int(row.eid), int(eid)) for eid, value in zip(matches.eid, cost)
            )
        if not pairs:
            raise ValueError(f"No eligible matched protein examples for sex {sex}")
        selected.extend(min(pairs)[1:])
    examples = population.set_index("eid").loc[selected].reset_index()
    examples["example"] = [
        "Female: lower disagreement",
        "Female: higher disagreement",
        "Male: lower disagreement",
        "Male: higher disagreement",
    ]
    references = [read_model(directory, fold, sex) for fold in range(1, 6) for sex in (0, 1)]
    mean_reference = np.mean([means / scales[None, :] for means, scales in references], axis=0)
    order = np.argsort(mean_reference[-1] - mean_reference[0], kind="stable")
    rows = []
    for example in examples.itertuples():
        means, scales = read_model(directory, example.OuterFold, example.Sex_male)
        index = np.flatnonzero(population.eid == example.eid)
        ll = likelihood(values[index], means, scales)[0]
        common = int(np.argmin(np.abs(AGE_GRID - example.common_age)))
        loss = ll.max(axis=1) - ll[:, common]
        rows.append(
            pd.DataFrame(
                {
                    "example": example.example,
                    "row": np.arange(1, len(features) + 1),
                    "protein": features[order],
                    "disagreement": loss[order],
                    "mean_disagreement": example.H,
                    "BAG": example.BAG,
                    "BAU": example.BAU,
                    "age": example.Age,
                    "common_age": example.common_age,
                }
            )
        )
    pd.concat(rows, ignore_index=True).to_csv(
        directory / "protein_examples.tsv", sep="\t", index=False
    )
    examples.to_feather(directory / "examples.feather")


if __name__ == "__main__":
    for modality in ("protein", "nmr"):
        score_modality(modality)
