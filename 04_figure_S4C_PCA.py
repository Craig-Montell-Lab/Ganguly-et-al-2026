"""Reproduce Figure S4C PCA; colors = neuron group, markers = hemisphere.

The ROI-record weight >= 2 threshold is applied before pair aggregation,
matching the author's ORIGINAL PCA code and archived PCA coordinates.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from sklearn.decomposition import PCA
from sklearn.preprocessing import normalize
from common import OUT, SEEDS, HEMISPHERES, ORIGINAL_COORDS, load_connectivity

MIN_WEIGHT = 2  # threshold for each ROI-level adjacency record, not pair total


def build_pca():
    df = load_connectivity()
    filtered = df.loc[df.weight >= MIN_WEIGHT].copy()
    mat = filtered.pivot_table(index=["bodyId_pre", "seed_group"],
                               columns="bodyId_post", values="weight",
                               aggfunc="sum", fill_value=0)
    expected = {(bid, group) for group, ids in SEEDS.items() for bid in ids}
    missing = expected - set(mat.index)
    if missing:
        raise ValueError(f"{len(missing)} seeds have no retained connections: {missing}")
    x = normalize(mat.values, norm="l1", axis=1)
    pca = PCA(n_components=2, random_state=42)
    z = pca.fit_transform(x)
    coords = pd.DataFrame({"bodyId_pre": mat.index.get_level_values("bodyId_pre"),
                           "seed_group": mat.index.get_level_values("seed_group"),
                           "PC1": z[:, 0], "PC2": z[:, 1]})
    return coords, pca


def main():
    coords, pca = build_pca()
    # Cross-check against the archived figure coordinates (up to numerical tolerance).
    if ORIGINAL_COORDS.exists():
        archive = pd.read_csv(ORIGINAL_COORDS)
        check = coords.merge(archive, on=["bodyId_pre", "seed_group"],
                             how="outer", validate="one_to_one", indicator=True,
                             suffixes=("_new", "_archive"))
        if (check._merge != "both").any() or not np.allclose(
            check[["PC1_new", "PC2_new"]], check[["PC1_archive", "PC2_archive"]], atol=1e-8):
            raise ValueError("Computed coordinates do not match archived Figure S4C coordinates")
        print("Verified PCA coordinates against archived original (tolerance 1e-8).")
    if not HEMISPHERES.exists():
        raise FileNotFoundError("Anatomical hemisphere reference CSV is required for Figure S4C")
    metadata = pd.read_csv(HEMISPHERES)[["bodyId_pre", "seed_group", "hemisphere"]]
    metadata = metadata.drop_duplicates(["bodyId_pre", "seed_group"])
    coords = coords.merge(metadata, on=["bodyId_pre", "seed_group"],
                          how="left", validate="one_to_one")
    if coords.hemisphere.isna().any() or not set(coords.hemisphere).issubset({"Left", "Right"}):
        raise ValueError("Missing or unexpected hemisphere assignments; check input metadata")
    coords.to_csv(OUT / "Figure_S4C_PCA_coordinates_with_hemisphere.csv", index=False)
    palette = dict(zip(SEEDS, plt.get_cmap("tab10").colors[:len(SEEDS)]))
    marker = {"Left": "s", "Right": "^"}
    fig, ax = plt.subplots(figsize=(7, 6))
    for group in SEEDS:
        for side in ("Left", "Right"):
            sub = coords[(coords.seed_group == group) & (coords.hemisphere == side)]
            ax.scatter(sub.PC1, sub.PC2, color=palette[group], marker=marker[side],
                       s=70, alpha=0.85, edgecolors="black", linewidths=0.35)
    ax.set_xlabel(f"PC1 ({100*pca.explained_variance_ratio_[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({100*pca.explained_variance_ratio_[1]:.1f}%)")
    ax.set_title("PCA of LB3a and LB3d subgroup downstream connectivity")
    group_handles = [Line2D([0], [0], marker="o", color="none", label=g,
                            markerfacecolor=palette[g], markeredgecolor="black", markersize=7)
                     for g in SEEDS]
    side_handles = [Line2D([0], [0], marker=marker[s], color="none", label=s,
                           markerfacecolor="gray", markeredgecolor="black", markersize=7)
                    for s in ("Left", "Right")]
    lg = ax.legend(handles=group_handles, title="Neuron group", frameon=False,
                   loc="upper left", bbox_to_anchor=(1.01, 1))
    ax.add_artist(lg)
    ax.legend(handles=side_handles, title="Hemisphere", frameon=False,
              loc="upper left", bbox_to_anchor=(1.01, .5))
    fig.tight_layout()
    fig.savefig(OUT / "Figure_S4C_PCA_hemisphere.pdf", bbox_inches="tight")
    fig.savefig(OUT / "Figure_S4C_PCA_hemisphere.png", dpi=600, bbox_inches="tight")
    plt.close(fig)
    print("PCA explained variance:", np.round(100*pca.explained_variance_ratio_, 4))
    print(pd.crosstab(coords.seed_group, coords.hemisphere))


if __name__ == "__main__":
    main()
