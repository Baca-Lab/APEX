<p align="center">
  <img src="man/figures/apex_gif-ezgif.com-optimize.gif" alt="APEX Logo" width="750">
</p>

<h1 align="center">APEX – Version 3.2.3 Released </h1>

<p align="center">
  <img src="https://img.shields.io/badge/status-active-brightgreen?style=flat-square">
  <img src="https://img.shields.io/badge/R->=4.4.0-blue?style=flat-square&logo=r">
  <!-- <img src="https://zenodo.org/badge/<ZENODO_ID>.svg"> -->
  <img src="https://img.shields.io/badge/docs-online-blue?style=flat-square&logo=readthedocs">
</p>

---

**APEX (Associating Peripheral Epigenomics with eXpression)** infers genome-wide tumor gene expression from plasma cfChIP-seq by integrating positional coverage and fragment-derived features from histone mark–enriched circulating chromatin (e.g., **H3K4me3**, **H3K36me3**) using pretrained models.  

## Installation
We recommend installing the APEX package using the `remotes` package from the R console. If you do not have `remotes` installed, you can install it by copying and pasting the following code in the R console:  

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("Baca-Lab/APEX")
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
  frag_k4   = c("S1.K4.bed", "S2.K4.bed"),
  frag_k27  = NA,
  frag_k36  = c("S1.K36.bed", "S2.K36.bed"),
  group     = c("Responder", "NonResponder")
)
```

## Quality control
Before feature extraction and expression inference, we recommend assessing cfChIP-seq library quality using histone mark–specific enrichment metrics. `apex_qc()` reports two complementary measures per sample and mark:  
  
	•	Fragment number: total uniquely mapped fragments (proxy for library complexity and sequencing depth)  
	•	Enrichment score: signal-to-noise metric comparing normalized coverage over expected on-target versus off-target genomic regions  
	
```r
qc <- apex_qc(manifest)
```

**Recommended QC thresholds**
These thresholds were used during model training and benchmarking and serve as practical guidelines (not strict cutoffs):  
	•	H3K4me3: enrichment > 7 and > 1 million fragments  
	•	H3K27ac: enrichment > 2 and > 1 million fragments  
	•	H3K36me3: enrichment > 2 and > 2 million fragments  

Samples below these thresholds may still be informative but should be interpreted with caution.  

## Infer gene expression

Once samples pass basic QC, APEX extracts epigenomic and fragmentomic features and infers genome-wide gene expression using pretrained models.  

For individual analysis, use apex(), and for cohort-level analyses, use `apex_batch()`:  

```r
#Single APEX run
apex_single <- apex(frag_file_k4 = "/PATH/TO/H3K4me3/FRAGMENT/FILE", frag_file_k36 = "/PATH/TO/H3K36me3/FRAGMENT/FILE")

#Batch APEX run
apex_mat <- apex_batch(manifest = manifest)
```

The result is a genes × samples matrix analogous to bulk RNA-seq expression data.  

## Differential gene expression analysis

Because APEX outputs inferred expression in a familiar matrix format, results can be analyzed using standard transcriptomic workflows. apex_diff() performs a limma-based differential analysis to estimate log₂ fold changes and moderated statistics between groups.  

```r
de <- apex_diff(apex_mat, group = manifest$group)

apex_volcano_plot(
  de,
  highlight_genes = c("NECTIN4", "TACSTD2", "ERBB2"),
  title = "Responders vs non-responders"
)
```

This plot summarizes gene-level differential expression with the option to highlight genes of interest.  


## Geneset analysis

In addition to individual genes, APEX supports pathway- and program-level analyses. `apex_geneset_score()` computes gene set activity scores (e.g., using ssGSEA) from APEX-inferred expression.  

```r
gs_scores <- apex_geneset_score(
  apex_mat,
  method = "ssgsea",
  genesets = list(
    EMT     = c("VIM", "FN1", "ZEB1", "TWIST1"),
    LUMINAL = c("UPK1A", "GATA3", "PPARG"),
    BASAL   = c("KRT5", "KRT14", "TP63")
  )
)
```

Gene set scores can be analyzed analogously to gene-level data, including differential analysis (with `apex_geneset_diff()`) and volcano-style visualization (with `apex_geneset_volcano_plot()`).  

---

# Nominating expression-based cancer targets

The `apex_rank_targets()` function computes gene-wise expression percentiles for user-supplied APEX-inferred expression data using either a pan-cancer or cancer-specific reference. Higher percentiles indicate genes that are unusually highly expressed relative to comparable tumors profiled by plasma cfChIP-seq, supporting context-aware target nomination directly from plasma.


```{r}
ranked <- apex_rank_targets(
  apex_mat,
  top_n = 5,
  plot = TRUE
)
```

---

# Documentation

A detailed vignette with example data and workflow can be accessed here:  
```r
browseVignettes("apex")
```

Other references:  
	•	SNAP pipeline: https://github.com/prc992/SNAP  
	•	Manuscript: in preparation  
	
---

# Citation

If you use APEX in your work, please cite the accompanying manuscript (details pending).  

---

