# APEX Snakemake Pipeline

This repository provides a minimal Snakemake workflow for running APEX on multiple samples in parallel and automatically compiling results into cohort-level outputs. The pipeline is designed for local execution or deployment on Slurm-based HPC systems.

The workflow:
1. Computes cohort-level QC metrics
2. Runs APEX independently per sample
3. Merges per-sample outputs into a single gene-by-sample expression matrix

The design prioritizes simplicity, transparency, and portability.

## Requirements

Software:
- R >= 4.4.0
- Snakemake
- APEX R package: https://github.com/Baca-Lab/APEX
- Standard Unix utilities (bash, awk)

Optional:
- Slurm (for cluster execution)

## Input

Manifest file

The pipeline is driven by a single tab-separated manifest file with four columns:

sample_id	frag_k4	frag_k27	frag_k36
SAMPLE_001	/path/to/K4.bed	NA	/path/to/K36.bed
SAMPLE_002	NA	NA	/path/to/K36.bed

Notes:
- sample_id must be unique
- Use NA for missing histone marks
- File paths must be accessible from compute nodes
- APEX-compatible BED files can be generated with the SNAPIE pipeline: https://github.com/prc992/SNAPIE

## Configuration

Edit config.yaml to point to your manifest and output directory:

manifest: manifests/manifest.tsv
outdir: results

## Running the pipeline

Dry run (recommended):

snakemake -n

Local execution:

snakemake -j 4

Slurm execution:

snakemake --cluster "sbatch -c 1 --mem 32G -t 00:30:00" -j 100

Each sample is processed as an independent job.

## Outputs

After completion, the output directory contains:

results/
├── APEX_QC_metrics.tsv
├── APEX_expression_matrix.tsv
└── per_sample/
    ├── SAMPLE_001.apex.tsv
    ├── SAMPLE_002.apex.tsv
    └── ...

QC metrics

APEX_QC_metrics.tsv includes enrichment scores and fragment counts for each sample and histone mark.

Expression matrix

APEX_expression_matrix.tsv contains inferred gene expression values (rows = genes, columns = samples) and can be analyzed using standard transcriptomic workflows.

## Notes

- QC is computed once per cohort prior to expression inference.
- APEX automatically handles missing histone marks specified as NA.
- The pipeline avoids containers and heavy orchestration to remain portable across environments.

## Citation

If you use this pipeline, please cite the APEX manuscript and relevant cfChIP-seq methodology papers.

## Contact

For questions or issues, please contact the APEX developers or open an issue on GitHub.
