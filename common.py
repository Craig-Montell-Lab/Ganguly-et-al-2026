"""Shared settings and validation for the PNAS connectome analysis scripts."""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)
DATASET = "male-cns:v0.9"
SEEDS = {
    "LB3a": [136183, 140619, 141663, 160435, 166190, 167663, 178913,
             183061, 187776, 199308, 203234, 374701, 518542,
             388541892, 456838775, 553738738, 786482749],
    "LB3d_ACh": [33339, 93290, 125405, 129171, 158537, 178321,
                 199196, 204726, 512552, 516218, 916953],
    "LB3d_glut": [130693, 140851, 201075, 236475, 261048, 516216, 550065],
    "LB3d_unannotated": [210344, 312102, 317616, 147619274,
                         245892505, 422510463, 581710868, 885642755],
}
RAW = DATA / "LB3a_LB3d_all_subgroups_raw_downstream_connectivity.csv"
HEMISPHERES = DATA / "LB3a_LB3d_PCA_hemisphere_data.csv"
ORIGINAL_COORDS = DATA / "LB3a_LB3d_subgroups_PCA_coordinates.csv"


def load_connectivity(path=RAW):
    df = pd.read_csv(path, dtype={"bodyId_pre": "int64", "bodyId_post": "int64"})
    required = {"bodyId_pre", "bodyId_post", "roi", "weight", "seed_group"}
    if not required.issubset(df.columns):
        raise ValueError(f"Missing raw connectivity columns: {required-set(df.columns)}")
    if df["weight"].isna().any() or (df["weight"] < 0).any():
        raise ValueError("Connectivity weights must be nonnegative and nonmissing")
    for grp, ids in SEEDS.items():
        actual = set(df.loc[df.seed_group.eq(grp), "bodyId_pre"])
        if not actual.issubset(set(ids)):
            raise ValueError(f"Unexpected seed IDs for {grp}: {actual-set(ids)}")
    return df


def aggregate_pairs(df):
    """Sum ROI-specific records for each directed neuron pair."""
    return (df.groupby(["bodyId_pre", "bodyId_post"], as_index=False,
                       sort=True)["weight"].sum())
