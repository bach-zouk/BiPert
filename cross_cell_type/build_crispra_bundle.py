"""build_crispra_bundle.py: packs the two per-guide mean-expression files of Southard et al. 2025 (TF CRISPRa screen in
fibroblasts and RPE-1 cells; files fibroblast_CRISPRa_mean_pop.h5ad and RPE1_CRISPRa_mean_pop.h5ad of their release)
into two small .npz files: the guides present in both cell types, on the genes present in both files.
Authors' record; a reader of the notebook does not need it. Needs h5py."""
import os, sys
import numpy as np, h5py
SRC = sys.argv[1]; OUT = sys.argv[2]; os.makedirs(OUT, exist_ok=True)


def load(fn):
    with h5py.File(f"{SRC}/{fn}", "r") as f:
        X = f["X"][:]
        rows = f["obs"][dict(f["obs"].attrs)["_index"]][:].astype(str)
        genes = np.char.upper(f["var"][dict(f["var"].attrs)["_index"]][:].astype(str))
    return X, rows, genes


Xa, ra, ga = load("fibroblast_CRISPRa_mean_pop.h5ad"); Xb, rb, gb = load("RPE1_CRISPRa_mean_pop.h5ad")
guides = np.intersect1d(ra, rb)
genes, ia, ib = np.intersect1d(ga, gb, return_indices=True)
pa = {r: i for i, r in enumerate(ra)}; pb = {r: i for i, r in enumerate(rb)}
A = Xa[[pa[g] for g in guides]][:, ia].astype(np.float32); B = Xb[[pb[g] for g in guides]][:, ib].astype(np.float32)
np.savez_compressed(f"{OUT}/fibroblast.npz", X=A, guides=guides, genes=genes)
np.savez_compressed(f"{OUT}/RPE1.npz", X=B, guides=guides, genes=genes)
print("guides per file", len(ra), len(rb), "| shared guides", len(guides), "| shared genes", len(genes), "| dtype", Xa.dtype,
      "| MB", [round(os.path.getsize(f"{OUT}/{n}.npz") / 1e6, 1) for n in ("fibroblast", "RPE1")])
