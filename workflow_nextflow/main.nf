#!/usr/bin/env nextflow
// scRNA cohort workflow — Nextflow port (DSL2).
//
// Nextflow is DATAFLOW-based: processes are connected by channels, and a process runs as soon as its
// inputs are available. Emitting a per-sample channel SCATTERS prepare/qc/cluster/markers across
// samples automatically; `.collect()` on the qc-metrics and cluster-summary channels is the GATHER —
// it blocks AGGREGATE until every sample has finished. Same scripts as the Snakemake port; only this
// orchestration layer differs.
//
//   JAVA_HOME=../.venv/lib/jvm PATH=../.venv/bin:$PATH \
//     nextflow run workflow_nextflow/main.nf -with-dag workflow_nextflow/dag.png

nextflow.enable.dsl = 2


process PREPARE {
    tag "$sample"
    publishDir "${params.outdir}/samples", mode: 'copy'
    input:
        tuple val(sample), val(fraction), val(seed)
    output:
        tuple val(sample), path("${sample}.h5ad")
    script:
        """
        python ${params.scripts}/prepare.py --sample ${sample} \
            --fraction ${fraction} --seed ${seed} --out ${sample}.h5ad
        """
}


process QC {
    tag "$sample"
    publishDir "${params.outdir}/qc", mode: 'copy'
    input:
        tuple val(sample), path(h5ad)
    output:
        tuple val(sample), path("${sample}_qc.h5ad"), emit: h5ad
        path "${sample}_metrics.tsv",                 emit: metrics
    script:
        """
        python ${params.scripts}/qc.py --input ${h5ad} --output ${sample}_qc.h5ad \
            --metrics ${sample}_metrics.tsv --min-genes ${params.min_genes} \
            --min-cells ${params.min_cells} --max-pct-mt ${params.max_pct_mt}
        """
}


process CLUSTER {
    tag "$sample"
    publishDir "${params.outdir}/clustered", mode: 'copy', pattern: "*.{h5ad,tsv}"
    publishDir "${params.outdir}/figures",   mode: 'copy', pattern: "*.png"
    input:
        tuple val(sample), path(h5ad)
    output:
        tuple val(sample), path("${sample}_clustered.h5ad"), emit: h5ad
        path "${sample}_summary.tsv",                        emit: summary
        path "${sample}_umap.png"
    script:
        """
        python ${params.scripts}/cluster.py --input ${h5ad} \
            --output ${sample}_clustered.h5ad --summary ${sample}_summary.tsv \
            --umap ${sample}_umap.png --n-hvg ${params.n_hvg} \
            --n-pcs ${params.n_pcs} --resolution ${params.resolution}
        """
}


process MARKERS {
    tag "$sample"
    publishDir "${params.outdir}/markers", mode: 'copy'
    input:
        tuple val(sample), path(h5ad)
    output:
        path "${sample}_markers.tsv"
    script:
        """
        python ${params.scripts}/markers.py --input ${h5ad} \
            --output ${sample}_markers.tsv --top ${params.top_markers}
        """
}


process AGGREGATE {
    publishDir "${params.outdir}", mode: 'copy'
    input:
        path metrics
        path summaries
    output:
        path "cohort_summary.tsv"
        path "cohort_summary.png"
    script:
        """
        python ${params.scripts}/aggregate.py --qc ${metrics} --summaries ${summaries} \
            --out-tsv cohort_summary.tsv --out-plot cohort_summary.png
        """
}


workflow {
    samples_ch = Channel.fromList(params.samples)        // scatter over samples
    prepared   = PREPARE(samples_ch)
    qcd        = QC(prepared)
    clustered  = CLUSTER(qcd.h5ad)
    MARKERS(clustered.h5ad)
    AGGREGATE(qcd.metrics.collect(), clustered.summary.collect())   // gather: waits for all samples
}
