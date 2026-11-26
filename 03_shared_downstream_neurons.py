"""Identify downstream neurons with direct input from BOTH LB3a and LB3d_glut."""
import pandas as pd
from common import OUT, load_connectivity


def main():
    df = load_connectivity()
    df = df[df.seed_group.isin(["LB3a", "LB3d_glut"])].copy()
    summary = (df.groupby(["bodyId_post", "seed_group"])
               .agg(total_weight=("weight", "sum"),
                    seed_count=("bodyId_pre", "nunique"))
               .reset_index())
    hubs = summary.pivot(index="bodyId_post", columns="seed_group",
                         values=["total_weight", "seed_count"]).fillna(0)
    hubs.columns = [f"{metric}_{group}" for metric, group in hubs.columns]
    hubs = hubs.reset_index()
    hubs = hubs[(hubs.total_weight_LB3a > 0) &
                (hubs.total_weight_LB3d_glut > 0)].copy()
    hubs["combined_weight"] = hubs.total_weight_LB3a + hubs.total_weight_LB3d_glut
    hubs["balanced_input_score"] = hubs[["total_weight_LB3a", "total_weight_LB3d_glut"]].min(axis=1)
    hubs["total_seed_count"] = hubs.seed_count_LB3a + hubs.seed_count_LB3d_glut
    hubs = hubs.sort_values(["balanced_input_score", "combined_weight", "total_seed_count"],
                            ascending=False)
    hubs.to_csv(OUT / "shared_hubs_LB3a_LB3d_glut.csv", index=False)
    print("Shared downstream neurons receiving direct input from both groups:", len(hubs))
    print("Saved:", OUT / "shared_hubs_LB3a_LB3d_glut.csv")


if __name__ == "__main__":
    main()
