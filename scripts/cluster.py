"""Stage 2 (cluster): normalize -> HVG -> PCA -> neighbors -> Leiden -> UMAP.

The canonical scanpy clustering pipeline. Writes the clustered object, a one-row summary TSV
(sample, n_cells, n_clusters), and a UMAP PNG colored by cluster.

  python scripts/cluster.py --input qc/s1.h5ad --output clustered/s1.h5ad \\
      --summary clustered/s1_summary.tsv --umap figures/s1_umap.png \\
      --n-hvg 2000 --n-pcs 30 --resolution 1.0
"""
from __future__ import annotations
import argparse
import matplotlib
matplotlib.use("Agg")
import scanpy as sc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--umap", required=True)
    ap.add_argument("--n-hvg", type=int, default=2000)
    ap.add_argument("--n-pcs", type=int, default=30)
    ap.add_argument("--resolution", type=float, default=1.0)
    a = ap.parse_args()

    adata = sc.read_h5ad(a.input)
    sample = str(adata.obs["sample"].iloc[0]) if "sample" in adata.obs else "NA"

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=a.n_hvg)
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable].copy()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=a.n_pcs)
    sc.pp.neighbors(adata, n_pcs=a.n_pcs)
    sc.tl.leiden(adata, resolution=a.resolution, flavor="igraph", n_iterations=2, directed=False)
    sc.tl.umap(adata)

    n_clusters = adata.obs["leiden"].nunique()
    sc.pl.umap(adata, color="leiden", title=f"{sample} ({adata.n_obs} cells, {n_clusters} clusters)",
               show=False, save=None)
    import matplotlib.pyplot as plt
    plt.savefig(a.umap, dpi=150, bbox_inches="tight")
    plt.close()

    adata.write(a.output)
    with open(a.summary, "w") as fh:
        fh.write("sample\tn_cells\tn_clusters\n")
        fh.write(f"{sample}\t{adata.n_obs}\t{n_clusters}\n")
    print(f"[cluster] {sample}: {adata.n_obs} cells -> {n_clusters} clusters")


if __name__ == "__main__":
    main()
