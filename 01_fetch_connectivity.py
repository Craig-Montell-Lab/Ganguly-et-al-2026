"""Optional: fetch the original downstream adjacency records from neuPrint.

The archived data/data/...csv is supplied, so an internet query is not required
for reproducing the analysis. Use this only to refresh the archived raw data.
"""
import os
from common import DATASET, RAW, SEEDS


def main():
    try:
        from neuprint import Client, fetch_adjacencies
        import pandas as pd
    except ImportError as exc:
        raise SystemExit("Install neuprint-python and pandas to query neuPrint") from exc
    token = os.environ.get("NEUPRINT_TOKEN")
    if not token:
        raise SystemExit("Set NEUPRINT_TOKEN to your personal token from neuprint.janelia.org")
    client = Client("neuprint.janelia.org", dataset=DATASET, token=token)
    tables = []
    for group, body_ids in SEEDS.items():
        _, adj = fetch_adjacencies(sources=body_ids, targets=None, client=client)
        adj["seed_group"] = group
        tables.append(adj[["bodyId_pre", "bodyId_post", "roi", "weight", "seed_group"]])
        print(group, "records:", len(adj))
    all_conn = pd.concat(tables, ignore_index=True)
    all_conn.to_csv(RAW, index=False)
    print("Saved:", RAW, "records:", len(all_conn))


if __name__ == "__main__":
    main()
