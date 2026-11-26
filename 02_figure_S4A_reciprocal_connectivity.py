"""Generate LB3a ↔ LB3d_glut direct connectivity heatmaps (Figure S4A).

IMPORTANT: sums ROI-specific records to obtain total synapses per neuron pair;
this corrects original pivot_table() mean aggregation. Shows all input neuron IDs,
and exports full matrices; figures display connected rows and columns as in the original script.
"""
import pandas as pd
import matplotlib.pyplot as plt
from common import OUT, SEEDS, load_connectivity, aggregate_pairs


def plot_matrix(pairs, pre_group, post_group, stem):
    source_ids, target_ids = SEEDS[pre_group], SEEDS[post_group]
    matrix = (pairs.pivot_table(index="bodyId_pre", columns="bodyId_post",
                                values="weight", aggfunc="sum", fill_value=0)
              .reindex(index=source_ids, columns=target_ids, fill_value=0))
    matrix.to_csv(OUT / f"{stem}_all_seeds_matrix.csv")
    # Match the original figure: plot only rows and columns with detected links.
    plotted = matrix.loc[matrix.sum(axis=1).gt(0), matrix.sum(axis=0).gt(0)]
    plotted.to_csv(OUT / f"{stem}_plotted_matrix.csv")
    fig, ax = plt.subplots(figsize=(3.4, 2.6))
    im = ax.imshow(plotted.values, aspect="auto", cmap="viridis", vmin=0)
    ax.set_xticks(range(plotted.shape[1]), labels=[str(i) for i in plotted.columns], rotation=90)
    ax.set_yticks(range(plotted.shape[0]), labels=[str(i) for i in plotted.index])
    ax.tick_params(labelsize=7)
    ax.set_xlabel(f"{post_group} bodyId")
    ax.set_ylabel(f"{pre_group} bodyId")
    ax.set_title(f"{pre_group} → {post_group}")
    fig.colorbar(im, ax=ax, label="Synaptic weight")
    fig.tight_layout()
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.png", dpi=600, bbox_inches="tight")
    plt.close(fig)
    return matrix


def main():
    df = load_connectivity()
    a = df[df.seed_group.eq("LB3a")]
    d = df[df.seed_group.eq("LB3d_glut")]
    fwd = aggregate_pairs(a[a.bodyId_post.isin(SEEDS["LB3d_glut"])])
    rev = aggregate_pairs(d[d.bodyId_post.isin(SEEDS["LB3a"])])
    fwd.to_csv(OUT / "LB3a_to_LB3d_glut_direct_pairs.csv", index=False)
    rev.to_csv(OUT / "LB3d_glut_to_LB3a_direct_pairs.csv", index=False)
    plot_matrix(fwd, "LB3a", "LB3d_glut", "Figure_S4A_top_LB3a_to_LB3d_glut")
    plot_matrix(rev, "LB3d_glut", "LB3a", "Figure_S4A_bottom_LB3d_glut_to_LB3a")
    fwd_ids = set(map(tuple, fwd[["bodyId_pre", "bodyId_post"]].to_numpy()))
    rev_ids = set(map(tuple, rev[["bodyId_post", "bodyId_pre"]].to_numpy()))
    reciprocal = fwd_ids & rev_ids
    summary = pd.DataFrame([
        {"direction": "LB3a_to_LB3d_glut", "connected_pairs": len(fwd), "total_weight": int(fwd.weight.sum())},
        {"direction": "LB3d_glut_to_LB3a", "connected_pairs": len(rev), "total_weight": int(rev.weight.sum())},
    ])
    summary.to_csv(OUT / "Figure_S4A_connectivity_summary.csv", index=False)
    pd.DataFrame(sorted(reciprocal), columns=["LB3a_bodyId", "LB3d_glut_bodyId"]).to_csv(
        OUT / "Figure_S4A_bidirectional_pairs.csv", index=False)
    print(summary.to_string(index=False))
    print("Unique bidirectionally connected neuron pairs:", len(reciprocal))


if __name__ == "__main__":
    main()
