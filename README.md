<p align="center">
  <img src="man/figures/APEX_logo_v2.png" alt="APEX Logo" width="750">
</p>
<h1 align="center">APEX – Version 4.0.0 Released </h1>

<p align="center">
  <img src="https://img.shields.io/badge/status-active-brightgreen?style=flat-square">
  <img src="https://img.shields.io/badge/R->=4.4.0-blue?style=flat-square&logo=r">
  <!-- <img src="https://zenodo.org/badge/<ZENODO_ID>.svg"> -->
  <img src="https://img.shields.io/badge/docs-online-blue?style=flat-square&logo=readthedocs">
</p>

---
**APEX (Associating Plasma Epigenomics with eXpression)** infers genome-wide tumor gene expression from plasma cfChIP-seq by integrating positional coverage and fragment-derived features from histone mark–enriched circulating chromatin (e.g., **H3K4me3**, **H3K36me3**) using pretrained models.  

## Installation

APEX depends on several Bioconductor packages that are not always automatically resolved by CRAN-based installers.  
We therefore recommend installing Bioconductor dependencies explicitly prior to installing APEX.

### Step 1: Install Bioconductor dependencies

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_packages <- c(
  "BiocGenerics", "XVector", "Biostrings", "BiocIO",
  "BSgenome.Hsapiens.UCSC.hg19", "GSVA", "limma",
  "rtracklayer", "S4Vectors", "GenomeInfoDb",
  "IRanges", "GenomicRanges", "BiocStyle"
)

BiocManager::install(bioc_packages, ask = FALSE, update = TRUE)
```

### Step 2: Install APEX from GitHub
We recommend installing APEX directly from GitHub using the 'remotes' package:

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github(
  "Baca-Lab/APEX",
  build_vignettes = TRUE,
  dependencies = TRUE
)
```

## Load apex
After installation, load the `apex` package using the `library` function in R:  

```r
library(apex)
```

---

# Minimal workflow

Below we provide a brief overview of the capabilites and functions available through the APEX R package. For a detailed example to work through, please refer to the tutorial vignette at:  

```r
browseVignettes("apex")
```

## Required input

APEX expects fragment-level BED-like files that include fragment coordinates and fragment-derived covariates (GC content, fragment length, end motifs). We recommend generating these files using the SNAP Nextflow pipeline, which is publicly available at [SNAP pipeline](https://github.com/prc992/SNAP).

**Important**: Please use the hg19 reference genome as APEX currently uses hg19 coordinates.  

Each fragment file contains:  
1.  chromosome  
2.	start  
3.	end  
4.	strand (optional)  
5.	fragment length  
6.	fragment GC content  
7.	5′ end motif (read 1)  
8.	5′ end motif (read 2)  

## Prepare a manifest

A manifest is a data.frame with one row per sample and paths to fragment files for each histone mark, along with a grouping variable for downstream comparisons.

```r
manifest <- data.frame(
  sample_id = c("S1", "S2"),
  frag_k4   = c("S1_K4_4NMER_bp_motif.bed", "S2_K4_4NMER_bp_motif.bed"),
  frag_k27  = NA,
  frag_k36  = c("S1_K36_4NMER_bp_motif.bed", "S2_K36_4NMER_bp_motif.bed"),
  group     = c("Responder", "NonResponder")
)
```

## Quality control
Before feature extraction and expression inference, we recommend assessing cfChIP-seq library quality using histone mark–specific enrichment metrics. `apex_qc()` reports two complementary measures per sample and mark:  
  
	•	Fragment number: total uniquely mapped fragments (proxy for library complexity and sequencing depth)  
	•	Enrichment score: signal-to-noise metric comparing normalized coverage over expected on-target versus off-target genomic regions  
	
When `plotQC = TRUE`, the function will generate boxplots and QC pass summaries for enrichment scores and fragment counts using mark-specific thresholds.

```r
qc <- apex_qc(manifest = manifest, plotQC = TRUE)
```

**Recommended QC thresholds**
These thresholds were used during model training and benchmarking and serve as practical guidelines (not strict cutoffs):  
	•	H3K4me3: enrichment > 7 and > 1 million fragments  
	•	H3K27ac: enrichment > 7 and > 1 million fragments  
	•	H3K36me3: enrichment > 2 and > 2 million fragments  

Samples below these thresholds may still be informative but should be interpreted with caution.  

## Infer gene expression
Once samples pass basic QC, APEX extracts epigenomic and fragmentomic features and infers genome-wide gene expression using pretrained models.

APEX provides multiple pretrained models corresponding to the chromatin immunoprecipitation data available for a given sample:

- **H3K4me3 only**
- **H3K36me3 only**
- **H3K4me3 + H3K36me3** *(recommended use case)*
- **H3K4me3 + H3K36me3 + H3K27ac**

APEX selects the appropriate model based on the fragment files supplied.

> **Note:**  
> H3K27ac alone is not provided as a standalone model, as it performed poorly in model evaluation when used in isolation.

### Single-sample analysis

```r
apex_single <- apex(
  frag_file_k4  = "/PATH/TO/H3K4me3/FRAGMENT/FILE",
  frag_file_k36 = "/PATH/TO/H3K36me3/FRAGMENT/FILE"
)
```

### Cohort-level analysis
```r
apex_mat <- apex_batch(manifest = manifest)
```

The output is a `genes × samples` matrix analogous to bulk RNA-seq expression data.


## Differential gene expression analysis

Because APEX outputs inferred expression in a familiar matrix format, results can be analyzed using standard transcriptomic workflows. `apex_diff()` performs a limma-based differential analysis to estimate log₂ fold changes and moderated statistics between groups.At least **three samples per group** are recommended to ensure stable variance estimation.

```{r}
de <- apex_diff(apex_mat, group = manifest$group)

apex_volcano_plot(
  de,
  highlight_genes = c("NECTIN4", "TACSTD2", "ERBB2"),
  title = "Responders vs non-responders"
)
```

This plot summarizes gene-level differential expression, with the option to highlight genes of interest.

## Gene set analysis

In addition to gene-level inference, APEX supports pathway- and program-level analyses. The function `apex_geneset_score()` computes gene set activity scores from APEX-inferred expression using either mean expression (`method = "mean"`) or single-sample gene set enrichment analysis (`method = "ssgsea"`), which estimates relative enrichment for each gene set independently per sample.

Gene sets may be provided directly (`geneset_source = "custom"`) or drawn from the MSigDB collection (`geneset_source = "msigdb"`), with optional specification of a particular collection or subcollection using `msigdb_collection` or `msigdb_subcollection`. The `min_genes` parameter ensures that only gene sets with sufficient representation in the inferred expression matrix are scored. See `?apex_geneset_score` and `?msigdbr::msigdbr` for additional details.

### Example: Custom gene sets
```{r}
gs_scores_custom <- apex_geneset_score(
  apex_mat,
  method = "ssgsea",
  geneset_source = "custom",
  genesets = list(
    EMT     = c("VIM", "FN1", "ZEB1", "TWIST1"),
    LUMINAL = c("UPK1A", "GATA3", "PPARG"),
    BASAL   = c("KRT5", "KRT14", "TP63")
  )
)
```

### Example: MSigDB Hallmark gene sets
```{r}
gs_scores_hallmark <- apex_geneset_score(
  apex_mat,
  method = "ssgsea",
  geneset_source = "msigdb",
  msigdb_collection = "H"
)
```

Gene set scores can be analyzed similarly to gene-level data, including differential testing (`apex_geneset_diff()`) and volcano-style visualization (`apex_geneset_volcano_plot()`).

---

# Nominating expression-based cancer targets

The `apex_rank_targets()` function computes gene-wise expression percentiles for user-supplied APEX-inferred expression data using either a pan-cancer or cancer-specific reference. Higher percentiles indicate genes that are unusually highly expressed relative to a heterogeneous cohort of tumors profiled by plasma cfChIP-seq.


```{r}
ranked <- apex_rank_targets(
  apex_mat,
  top_n = 5,
  plot = TRUE
)
```

---

# Snakemake workflow for scalable execution

APEX includes a reference Snakemake workflow for running analyses across multiple samples in parallel. The workflow supports local execution and is compatible with Slurm-based HPC systems for large-scale runs.

The Snakemake workflow, along with scripts and usage instructions, is provided in:

`inst/workflows/snakemake/`

A dedicated README in that directory describes required inputs, configuration, and commands for execution.

---

# Documentation

A detailed vignette with example data and workflow can be accessed here:  
```r
browseVignettes("apex")
```

Other references:  
	•	SNAPIE pipeline: https://github.com/prc992/SNAPIE
	•	Manuscript: in preparation  
	
---

# Citation

If you use APEX in your work, please cite the accompanying manuscript (details pending).  

---

