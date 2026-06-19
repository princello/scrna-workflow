"""Stage 4 (aggregate / GATHER): combine all per-sample results into one cohort report.

This is the fan-in: it takes every sample's QC-metrics and cluster-summary TSVs and produces a single
cohort summary table plus a bar plot. This stage is why a workflow engine earns its keep — it can't run
until ALL samples have finished, and the engine tracks that dependency for you.

  python scripts/aggregate.py --qc qc/s1_metrics.tsv qc/s2_metrics.tsv ... \\
      --summaries clustered/s1_summary.tsv ... --out-tsv results/cohort_summary.tsv \\
      --out-plot results/cohort_summary.png
"""
from __future__ import annotations
import argparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--qc", nargs="+", required=True)
    ap.add_argument("--summaries", nargs="+", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-plot", required=True)
    a = ap.parse_args()

    qc = pd.concat([pd.read_csv(f, sep="\t") for f in a.qc], ignore_index=True)
    cl = pd.concat([pd.read_csv(f, sep="\t") for f in a.summaries], ignore_index=True)
    cohort = qc.merge(cl[["sample", "n_clusters"]], on="sample").sort_values("sample")
    cohort.to_csv(a.out_tsv, sep="\t", index=False)
    print("[aggregate] cohort summary:")
    print(cohort.to_string(index=False))

    fig, axes = plt.subplots(1, 3, figsize=(11, 3.2))
    for ax, col, title in zip(
        axes,
        ["cells_after", "median_genes_per_cell", "n_clusters"],
        ["cells passing QC", "median genes / cell", "Leiden clusters"],
    ):
        ax.bar(cohort["sample"], cohort[col], color="#4C78A8")
        ax.set_title(title)
        ax.set_xlabel("sample")
    fig.suptitle("scRNA cohort summary")
    fig.tight_layout()
    fig.savefig(a.out_plot, dpi=150, bbox_inches="tight")
    print(f"[aggregate] wrote {a.out_tsv} and {a.out_plot}")


if __name__ == "__main__":
    main()
