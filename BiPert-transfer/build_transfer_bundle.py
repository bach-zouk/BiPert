"""Build the compact data bundle read by linear_transfer_claims.ipynb.

Input (--src): one pseudobulk file per cell line, prep_<line>.npz, with the keys
    genes, symbols   Ensembl id and symbol per gene
    perts            knockdown names (target gene symbols)
    pb_mean          mean log(1 + x) profile per knockdown (counts normalised to 10,000 per cell), knockdowns x genes
    ctrl             log(1 + x) profiles of non-targeting control cells, cells x genes (at most ~12,000 per line)
    hv_sum, hv_sq    per-gene sum and sum of squares of the normalised counts over all cells of the line
    n_all            number of cells of the line
Output (--out): one float32 .npz per cell line, restricted to the benchmark's knockdowns and gene panel, plus
MANIFEST.json. Nothing else is needed to re-derive the perturbation-transfer numbers.

Two benchmarks, built exactly as in the analysis behind the paper:
    four_lines   K562, RPE1, HepG2, Jurkat. Knockdowns measured in all four lines; genes present in all four files;
                 panel = 2,000 most variable genes of the four lines pooled.
    orion_pair   HCT116, HEK293T. Knockdowns measured in both lines and in RPE1, HepG2 and Jurkat (2,270); panel =
                 2,000 most variable genes of the two lines pooled.
The most variable genes are ranked by the Seurat normalised dispersion (scanpy's flavor="seurat"), computed from the
pooled per-gene mean and variance of the normalised counts.

Usage: python build_transfer_bundle.py --src <folder with prep_*.npz> --out transfer_data
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

N_PANEL = 2000
BENCHMARKS = {
    "four_lines": {"lines": ["K562", "RPE1", "HepG2", "Jurkat"], "files": ["K562", "RPE1", "HepG2", "Jurkat"],
                   "also_measured_in": []},
    "orion_pair": {"lines": ["HCT116", "HEK293T"], "files": ["HCT116all", "HEK293Tall"],
                   "also_measured_in": ["RPE1", "HepG2", "Jurkat"]},
}


def seurat_hvg(mean, var, n_top):
    """Indices of the n_top genes with the highest dispersion, normalised within 20 bins of mean expression."""
    mean = np.where(mean == 0, 1e-12, mean)
    disp = var / mean
    disp = np.log(np.where(disp == 0, np.nan, disp))
    df = pd.DataFrame({"m": np.log1p(mean), "d": disp})
    df["bin"] = pd.cut(df.m, bins=20)
    g = df.groupby("bin", observed=True).d
    mu, sd = g.transform("mean"), g.transform("std")
    single = sd.isna()
    sd, mu = sd.where(~single, mu), mu.where(~single, 0)
    return np.argsort(-((df.d - mu) / sd).fillna(-np.inf).to_numpy())[:n_top]


def build(src, out, name, lines, files, also_measured_in):
    raw = {c: np.load(src / f"prep_{f}.npz", allow_pickle=True) for c, f in zip(lines, files)}
    names = {c: {k: raw[c][k].astype(str) for k in ("genes", "symbols", "perts")} for c in lines}
    genes = names[lines[0]]["genes"]
    for c in lines[1:]:
        genes = genes[np.isin(genes, names[c]["genes"])]
    perts = set.intersection(*[set(names[c]["perts"]) for c in lines])
    for f in also_measured_in:
        perts &= set(np.load(src / f"prep_{f}.npz", allow_pickle=True)["perts"].astype(str))
    perts = np.array(sorted(perts))
    col = {c: pd.Index(names[c]["genes"]).get_indexer(genes) for c in lines}
    n = sum(int(raw[c]["n_all"]) for c in lines)
    s = sum(raw[c]["hv_sum"][col[c]] for c in lines)
    q = sum(raw[c]["hv_sq"][col[c]] for c in lines)
    panel = np.sort(seurat_hvg(s / n, (q - s ** 2 / n) / (n - 1), N_PANEL))
    symbols = names[lines[0]]["symbols"][col[lines[0]]][panel]
    entries = []
    for c in lines:
        row = pd.Index(names[c]["perts"]).get_indexer(perts)
        assert (row >= 0).all() and (col[c] >= 0).all()
        pert_mean = raw[c]["pb_mean"][row][:, col[c][panel]].astype(np.float32)
        ctrl = raw[c]["ctrl"][:, col[c][panel]].astype(np.float32)
        path = out / f"{c}.npz"
        np.savez_compressed(path, pert_mean=pert_mean, ctrl=ctrl, perts=perts, genes=genes[panel], symbols=symbols)
        entries.append({"file": path.name, "line": c, "benchmark": name, "knockdowns": int(pert_mean.shape[0]),
                        "genes": int(pert_mean.shape[1]), "control_cells": int(ctrl.shape[0]),
                        "bytes": path.stat().st_size, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
        print(f"{name}/{c}: {pert_mean.shape[0]} knockdowns x {pert_mean.shape[1]} genes, {ctrl.shape[0]} control cells, "
              f"{path.stat().st_size / 1e6:.1f} MB", flush=True)
    return entries


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--src", required=True, type=Path, help="folder with the per-line prep_<line>.npz files")
    ap.add_argument("--out", default=Path("transfer_data"), type=Path, help="output folder")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    entries = []
    for name, spec in BENCHMARKS.items():
        entries += build(a.src, a.out, name, **spec)
    manifest = {"description": "Per-line knockdown mean profiles and control cells on a 2,000-gene panel; "
                               "values are log(1 + x) of counts normalised to 10,000 per cell, float32.",
                "keys": {"pert_mean": "knockdowns x genes, mean profile of the cells carrying each knockdown",
                         "ctrl": "cells x genes, non-targeting control cells",
                         "perts": "knockdown names (target gene symbols), same order in every line of a benchmark",
                         "genes": "Ensembl gene ids of the panel", "symbols": "gene symbols of the panel"},
                "files": entries, "total_bytes": sum(e["bytes"] for e in entries)}
    (a.out / "MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"total {manifest['total_bytes'] / 1e6:.1f} MB", flush=True)


if __name__ == "__main__":
    main()
