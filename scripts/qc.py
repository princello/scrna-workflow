"""Stage 1 (qc): per-cell / per-gene quality control.

Standard scRNA QC: drop cells with too few detected genes, drop genes seen in too few cells, and drop
cells with a high mitochondrial fraction (a dying-cell signature). Writes the filtered matrix plus a
one-row metrics TSV so the aggregate stage can report cohort-wide QC.

  python scripts/qc.py --input samples/s1.h5ad --output qc/s1.h5ad --metrics qc/s1_metrics.tsv \\
      --min-genes 200 --min-cells 3 --max-pct-mt 10
"""
from __future__ import annotations
import argparse
import scanpy as sc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--metrics", required=True)
    ap.add_argument("--min-genes", type=int, default=200)
    ap.add_argument("--min-cells", type=int, default=3)
    ap.add_argument("--max-pct-mt", type=float, default=10.0)
    a = ap.parse_args()

    adata = sc.read_h5ad(a.input)
    sample = str(adata.obs["sample"].iloc[0]) if "sample" in adata.obs else "NA"
    n0 = adata.n_obs

    sc.pp.filter_cells(adata, min_genes=a.min_genes)
    sc.pp.filter_genes(adata, min_cells=a.min_cells)
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None)
    adata = adata[adata.obs["pct_counts_mt"] < a.max_pct_mt].copy()

    adata.write(a.output)
    median_genes = float(adata.obs["n_genes_by_counts"].median())
    with open(a.metrics, "w") as fh:
        fh.write("sample\tcells_before\tcells_after\tgenes\tmedian_genes_per_cell\n")
        fh.write(f"{sample}\t{n0}\t{adata.n_obs}\t{adata.n_vars}\t{median_genes:.0f}\n")
    print(f"[qc] {sample}: {n0} -> {adata.n_obs} cells (median {median_genes:.0f} genes/cell)")


if __name__ == "__main__":
    main()
