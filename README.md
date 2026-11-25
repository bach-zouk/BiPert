# BiPert: How Far Can Linear Go? BiPert, a Unified Framework for Single-Cell Perturbation Modelling

Notebooks that recompute the results from compact data bundles. They need only numpy, scipy, pandas and matplotlib, run on a
CPU, and fit a free Google Colab session.

| Notebook | Open | What it reproduces | Data bundle | Runtime (2 CPU cores) | Peak memory |
|---|---|---|---|---|---|
| `BiPert-transfer/BiPert_transfer.ipynb` | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/bach-zouk/BiPert/blob/main/BiPert-transfer/BiPert_transfer.ipynb) | **Perturbation transfer** (BiPert, BiPert-FS): Fig. 3, the rank ablation, the JIVE comparison and the few-shot Table 1 | `transfer_data/`, 0.23 GB | about 25 min | under 10 GB |
| `BiPert-ID/BiPert_ID_identification.ipynb` | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/bach-zouk/BiPert/blob/main/BiPert-ID/BiPert_ID_identification.ipynb) | **Perturbation identification** (BiPert-ID): Tables 2 to 6, Fig. S2, Tables S2 and S3, and the intervals of Appendix F, on ExPert's full benchmark | `identification_data/`, 1.5 GB | about 7 min | about 6.5 GB |
| `cross_cell_type/BiPert_cross_cell_type.ipynb` | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/bach-zouk/BiPert/blob/main/cross_cell_type/BiPert_cross_cell_type.ipynb) | **Identifiability within a cell type limits transfer**: Table 7 and the paragraph on guides in Appendix F | `crispra_data/`, 0.34 GB | under 1 min | under 2 GB |

Each notebook defines the model in short cells, then has one section per result. Every section ends with a table
`paper value | recomputed | difference | OK`, and the last cell lists every check. 

## Data bundles

**Data coming soon.** The bundles are too large for a git repository and will be linked here. Put each folder
next to its notebook:

```
BiPert-transfer/transfer_data/
BiPert-ID/identification_data/
cross_cell_type/crispra_data/
```

## How to run

**Locally**

```
pip install -r requirements.txt
jupyter nbconvert --to notebook --execute --output my_run.ipynb BiPert_transfer.ipynb
```

or open a notebook in Jupyter and run all cells.

**Google Colab**

Copy the data folder into a folder `BiPert` on your Google Drive (`MyDrive/BiPert/<data folder>`), open the notebook in Colab, and in the configuration cell uncomment
the two lines that mount the Drive and set `DATA_DIR`.

**Quick mode.** `FAST = True` in the configuration cell of the transfer and identification notebooks runs a
labelled subset (transfer: fewer random splits, a few minutes; identification: only the four cell lines of
Table 2).

## Data sources

* Transfer benchmark (COMPASS's benchmark): Replogle et al. 2022 (K562, RPE1), Nadig et al. 2025 (HepG2, Jurkat),
  X-Atlas/Orion (HCT116, HEK293T).
  
* Identification benchmark: ExPert's benchmark of 22 datasets, rebuilt from the public cell-level data with
  ExPert's quality filters; the candidate list, folds, test cells and ExPert's released results are ExPert's.
  
* Cross-cell-type analysis: Southard et al., *Nature Genetics* 2025 (TF CRISPR activation in fibroblasts and
  RPE-1 cells), the authors' released per-guide mean profiles.

Please cite the original datasets and benchmarks when you use the bundles.

## Licence

The code is released under the MIT licence. The data bundles will carry their own terms, which follow those of the
source datasets. Please cite our work if you use the code/its findings.
