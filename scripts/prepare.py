"""Stage 0 (prepare): materialize one input sample as an .h5ad.

For a self-contained, fast, reproducible demo we derive several 'samples' by deterministically
subsampling the classic 10x pbmc3k dataset (different seed/fraction per sample, set in config). In a
real pipeline this rule would instead point at each sample's CellRanger output. Producing one file per
sample is what lets the workflow SCATTER the downstream stages across samples in parallel.

  python scripts/prepare.py --sample s1 --fraction 0.8 --seed 1 --out samples/s1.h5ad
"""
from __future__ import annotations
import argparse
import numpy as np
import scanpy as sc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--fraction", type=float, default=0.8)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    adata = sc.datasets.pbmc3k()                 # ~2700 cells x 32k genes (downloaded + cached once)
    rng = np.random.default_rng(a.seed)
    n = int(adata.n_obs * a.fraction)
    idx = rng.choice(adata.n_obs, size=n, replace=False)
    adata = adata[np.sort(idx)].copy()
    adata.obs["sample"] = a.sample
    adata.write(a.out)
    print(f"[prepare] {a.sample}: wrote {adata.n_obs} cells x {adata.n_vars} genes -> {a.out}")


if __name__ == "__main__":
    main()
