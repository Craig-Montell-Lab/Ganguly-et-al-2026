# Ganguly et al 2026: connectome analysis code

This repository reproduces the LB3a/LB3d connectome analyses associated with Ganguly et al.'s manuscript titled "Water satiety requires Otopetrin-mediated osmosensing
in inhibitory taste neurons", including **Supplementary Figure S4A**, the count of **185 shared downstream neurons**, and **Supplementary Figure S4C** (PCA and statistical tests).

## Data source and original design

- Source: [neuPrint](https://neuprint.janelia.org/), dataset **`male-cns:v0.9`** (*Drosophila* male CNS).
- Seed groups: **LB3a (17)**, **LB3d_ACh (11)**, **LB3d_glut (7)**, and **LB3d_unannotated (8)**. The full body IDs are documented in `common.py`.
- Raw adjacency records, including `roi` and `weight`, are supplied in `data/LB3a_LB3d_all_subgroups_raw_downstream_connectivity.csv` (3,080 rows). These data were supplied by the authors; no neuPrint query is needed to run the downstream scripts.
- `data/LB3a_LB3d_PCA_hemisphere_data.csv` contains the author's anatomical left/right annotations merged with original PCA coordinates. These annotations were provided by the authors and are **not independently revalidated** by these offline scripts. Confirm provenance of neuPrint `somaSide`/`rootSide` annotations before submission.
- The separately supplied archived PCA coordinates are used as a reproducibility reference and as the input to the statistics script.

## Files

| Script | Purpose | Main output(s) |
| --- | --- | --- |
| `01_fetch_connectivity.py` | OPTIONAL neuPrint refetch of the raw adjacencies for all four seed groups | `data/LB3a_LB3d_all_subgroups_raw_downstream_connectivity.csv` |
| `02_figure_S4A_reciprocal_connectivity.py` | LB3a → LB3d_glut and LB3d_glut → LB3a connectivity heatmaps, summaries and reciprocal pair list | `outputs/Figure_S4A_*.pdf`, `*.png`, `*.csv` |
| `03_shared_downstream_neurons.py` | Intersect postsynaptic targets receiving direct inputs from both LB3a and LB3d_glut | `outputs/shared_hubs_LB3a_LB3d_glut.csv` |
| `04_figure_S4C_PCA.py` | Reproduce the original individual-neuron PCA; add squares for anatomical left and triangles for anatomical right | `outputs/Figure_S4C_PCA_hemisphere.pdf`, `*.png`, `*_coordinates_with_hemisphere.csv` |
| `05_figure_S4C_statistics.py` | Group separation permutations; two-sided pairwise Mann–Whitney U; hemisphere-stratified permutation sensitivity test | `outputs/Figure_S4C_permutation_tests.csv`, `Figure_S4C_Mann_Whitney_tests.csv` |
| `common.py` | Shared body IDs, paths, data validation | — |

## Installation

Python 3.10+ recommended. From this directory:

```bash
python -m pip install -r requirements.txt
```

The `neuprint-python` package is required only if you wish to refetch neuPrint records. `numpy`, `pandas`, `scipy`, `scikit-learn` and `matplotlib` are used for offline reproduction. For archival reproducibility, record the actual package versions in the GitHub/Zenodo release; this `requirements.txt` currently does **not** pin versions.

## Running the analyses

From this directory:

```bash
python 02_figure_S4A_reciprocal_connectivity.py
python 03_shared_downstream_neurons.py
python 04_figure_S4C_PCA.py
python 05_figure_S4C_statistics.py
```

All generated files are saved under `outputs/`. The supplied data allow these scripts to run **without a neuPrint token**.

### Optional raw-data retrieval

Create a personal API token at [neuprint.janelia.org](https://neuprint.janelia.org/). Set the environment variable `NEUPRINT_TOKEN` (do **not** commit your token), then run:

```bash
python 01_fetch_connectivity.py
```

The optional fetch overwrites the archived raw CSV. Back up the provided data first if you want an exact comparison with the original analysis. The result may change if the underlying dataset or service changes.

## Methods and output interpretations

### Figure S4A: reciprocal-connectivity heatmaps

Both directions are tested: LB3a → LB3d_glut and LB3d_glut → LB3a. ROI-specific adjacency records for each directed pair are **summed across ROIs**. For each direction, the exported `*_all_seeds_matrix.csv` includes all input neurons and explicit zeros for absent reported connections. To match the existing Figure S4A design, the plotted heatmaps show only connected rows/columns, with the exact plotted values separately saved as `*_plotted_matrix.csv`. The exploratory script used `pivot_table()`'s default **mean**; this has been replaced with **sum**. **In the supplied S4A data there are no duplicate directed neuron-pair ROI records, so this aggregation correction does not alter the plotted nonzero synaptic weights.** Nevertheless, compare output files with the manuscript panels before submission. “Reciprocal populations” does not imply that every individual neuron pair has bidirectional connections; a separate CSV lists confirmed bidirectional pairs.

### Shared downstream neuron count

Targets are retained when they have positive summed synaptic weight from **at least one LB3a neuron AND at least one LB3d_glut neuron**. The archived data yield **185 distinct postsynaptic body IDs**. We retain the author's `combined_weight`, `balanced_input_score`, and `total_seed_count` descriptive rankings. “Shared downstream neuron” is more precise than calling each candidate a network hub, because a centrality calculation was not performed.

### Figure S4C: individual-neuron PCA

The original analysis filter is preserved exactly: `MIN_WEIGHT = 2` is applied to **each ROI adjacency record before aggregation**. A matrix of presynaptic neurons × postsynaptic target IDs is assembled by **summing retained records per pair**. Each row is L1-normalized, followed by `PCA(n_components=2, random_state=42)` on 43 individual neurons. PC1 and PC2 explain approximately **23.7%** and **14.2%** of the variance, respectively, matching the archived PCA coordinate file to numerical precision. Colors indicate four seed groups; **squares indicate left** and **triangles right** anatomical hemisphere according to supplied metadata. Hemisphere must never be assigned by the sign or rank of PC1.

### PCA statistics

On the archived Figure S4C PC1–PC2 coordinates, the primary permutation statistic is:

`mean Euclidean distance between groups / mean Euclidean distance within groups`

All unordered neuron pairs are used. Group labels are shuffled **9,999** times using random seed **42**, preserving group sizes. Monte Carlo P-values use `(1 + number_of_permuted_statistics_at_least_observed)/(9,999 + 1)`. Separately, for each PC coordinate, all six pairwise group comparisons are tested using **two-sided Mann–Whitney U** (12 tests in all). **These P-values are uncorrected for multiple comparisons**, so this must be explicitly reported if they are cited.

An additional **hemisphere-stratified permutation sensitivity test** shuffles neuron-group labels *within* each anatomical hemisphere (also 9,999 permutations, seed 42), keeping hemisphere composition fixed. This checks for residual group separation conditional on laterality, but does not show that hemisphere has no effect. The stratified test is a **new supplemental analysis**, not the original reported test; only report it in a manuscript if you choose to incorporate its findings.





