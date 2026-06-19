"""Stage 3 (markers): rank marker genes per cluster.

Wilcoxon rank-sum test (scanpy's rank_genes_groups) for each Leiden cluster vs the rest; writes the
top-N markers per cluster as a tidy TSV.

  python scripts/markers.py --input clustered/s1.h5ad --output markers/s1.tsv --top 10
"""
from __future__ import annotations
import argparse
import scanpy as sc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--top", type=int, default=10)
    a = ap.parse_args()

    adata = sc.read_h5ad(a.input)
    sample = str(adata.obs["sample"].iloc[0]) if "sample" in adata.obs else "NA"
    # rank_genes_groups wants log-normalized values; those live in .raw after cluster.py
    use = adata.raw.to_adata() if adata.raw is not None else adata
    use.obs["leiden"] = adata.obs["leiden"].values
    sc.tl.rank_genes_groups(use, "leiden", method="wilcoxon")

    with open(a.output, "w") as fh:
        fh.write("sample\tcluster\trank\tgene\tlog2fc\tpval_adj\n")
        for cl in use.obs["leiden"].cat.categories:
            df = sc.get.rank_genes_groups_df(use, group=cl).head(a.top).reset_index(drop=True)
            for r, row in df.iterrows():
                fh.write(f"{sample}\t{cl}\t{r+1}\t{row['names']}\t{row['logfoldchanges']:.2f}\t"
                         f"{row['pvals_adj']:.2e}\n")
    print(f"[markers] {sample}: top {a.top} markers x {len(use.obs['leiden'].cat.categories)} clusters -> {a.output}")


if __name__ == "__main__":
    main()
