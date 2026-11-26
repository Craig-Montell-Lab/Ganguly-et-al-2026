"""Statistical tests on the exact archived Figure S4C PC1–PC2 coordinates.

Tests:
 * 9,999-permutation global between/within Euclidean distance ratio
 * 9,999-permutation sensitivity test, labels shuffled within hemispheres
 * all 12 two-sided Mann–Whitney U tests for group pairs (PC1 and PC2)

Anatomical hemisphere labels are read from archived metadata; never inferred
from PC1 coordinates. No side-collapsed PCA is claimed or recomputed here.
"""
import numpy as np
import pandas as pd
from itertools import combinations
from scipy.spatial.distance import pdist
from scipy.stats import mannwhitneyu
from common import OUT, SEEDS, ORIGINAL_COORDS, HEMISPHERES

N_PERMUTATIONS = 9999
RANDOM_SEED = 42


def validate_and_load():
    pc = pd.read_csv(ORIGINAL_COORDS)
    hemi = pd.read_csv(HEMISPHERES)[["bodyId_pre", "seed_group", "PC1", "PC2", "hemisphere"]]
    merged = pc.merge(hemi, on=["bodyId_pre", "seed_group"],
                      suffixes=("", "_annotated"), how="inner", validate="one_to_one")
    if len(merged) != len(pc) or len(merged) != len(hemi):
        raise ValueError("PCA and hemisphere files do not have matching neuron IDs")
    if not np.allclose(merged[["PC1", "PC2"]],
                       merged[["PC1_annotated", "PC2_annotated"]], atol=1e-10):
        raise ValueError("Archived PCA and hemisphere-labeled coordinates disagree")
    if not set(merged.hemisphere).issubset({"Left", "Right"}):
        raise ValueError("Expected anatomical Left/Right annotations")
    for group, ids in SEEDS.items():
        if set(merged.loc[merged.seed_group == group, "bodyId_pre"]) != set(ids):
            raise ValueError(f"Incomplete or incorrect group membership: {group}")
    return merged


def main():
    df = validate_and_load()
    labels = df.seed_group.to_numpy()
    side = df.hemisphere.to_numpy()
    x = df[["PC1", "PC2"]].to_numpy()
    # Condensed pair distances and pair-membership masks avoid redundant work.
    dist = pdist(x, metric="euclidean")
    i, j = np.triu_indices(len(df), 1)

    def distance_ratio(g):
        same = g[i] == g[j]
        within = dist[same].mean()
        between = dist[~same].mean()
        return between / within, within, between

    observed, within, between = distance_ratio(labels)
    rng = np.random.RandomState(RANDOM_SEED)  # matches legacy seed=42 permutation sequence
    null = np.array([distance_ratio(rng.permutation(labels))[0]
                     for _ in range(N_PERMUTATIONS)])
    # Include the observed configuration in Monte Carlo p calculation.
    p_global = (1 + np.count_nonzero(null >= observed)) / (N_PERMUTATIONS + 1)

    # Secondary sensitivity analysis: keeps anatomical hemisphere membership fixed.
    rng = np.random.RandomState(RANDOM_SEED)
    side_indices = [np.flatnonzero(side == s) for s in ("Left", "Right")]
    hemi_null = []
    for _ in range(N_PERMUTATIONS):
        perm = labels.copy()
        for ix in side_indices:
            perm[ix] = rng.permutation(perm[ix])
        hemi_null.append(distance_ratio(perm)[0])
    p_hemi = (1 + np.count_nonzero(np.asarray(hemi_null) >= observed)) / (N_PERMUTATIONS + 1)

    statistics = pd.DataFrame([
        {"test": "global_group_separation", "observed_ratio": observed,
         "mean_within_distance": within, "mean_between_distance": between,
         "n_permutations": N_PERMUTATIONS, "random_seed": RANDOM_SEED,
         "p_value": p_global, "permutation_scheme": "shuffle_group_labels_over_all_neurons"},
        {"test": "group_separation_conditioned_on_hemisphere", "observed_ratio": observed,
         "mean_within_distance": within, "mean_between_distance": between,
         "n_permutations": N_PERMUTATIONS, "random_seed": RANDOM_SEED,
         "p_value": p_hemi, "permutation_scheme": "shuffle_group_labels_within_each_hemisphere"},
    ])
    statistics.to_csv(OUT / "Figure_S4C_permutation_tests.csv", index=False)
    print(statistics.to_string(index=False))

    results = []
    for axis in ("PC1", "PC2"):
        for g1, g2 in combinations(SEEDS, 2):
            a = df.loc[df.seed_group == g1, axis].to_numpy()
            b = df.loc[df.seed_group == g2, axis].to_numpy()
            u, p = mannwhitneyu(a, b, alternative="two-sided")
            results.append({"PC": axis, "group1": g1, "group2": g2,
                            "n1": len(a), "n2": len(b), "U_statistic": u,
                            "p_value_uncorrected": p})
    mw = pd.DataFrame(results)
    mw.to_csv(OUT / "Figure_S4C_Mann_Whitney_tests.csv", index=False)
    print("\nMann–Whitney comparisons (uncorrected two-sided p-values):")
    print(mw.to_string(index=False))
    print("\nAnatomical hemisphere counts:")
    print(pd.crosstab(df.seed_group, df.hemisphere).to_string())
    print("\nInterpretation: the hemisphere-stratified test examines whether group")
    print("separation exceeds that expected after preserving hemisphere composition.")
    print("It does NOT prove that hemisphere has no influence on the PCA.")


if __name__ == "__main__":
    main()
