# Installing APEX dependencies

APEX does not download or install third-party R packages. Users are responsible
for reviewing and installing the dependencies below before installing APEX.
Package names and minimum versions are declared in [`DESCRIPTION`](DESCRIPTION),
which is the authoritative dependency list for the current APEX release.

## Requirements

- R 4.4.0 or newer
- Bioconductor 3.20 or a version compatible with the installed R release
- Any operating-system libraries required by the selected CRAN and Bioconductor
  packages

## 1. Install the package-management tools

Install these tools only if they are not already available:

```r
install.packages(c("BiocManager", "remotes"))
```

## 2. Install required CRAN packages

```r
cran_packages <- c(
  "magrittr",
  "data.table",
  "reshape2",
  "dplyr",
  "xgboost",
  "ggrepel",
  "ggplot2",[/']'
  "msigdbr",
  "ggpubr",
  "plotly"
)

install.packages(cran_packages)
```

## 3. Install required Bioconductor packages

```r
bioc_packages <- c(
  "GenomicRanges",
  "IRanges",
  "GenomeInfoDb",
  "S4Vectors",
  "rtracklayer",
  "limma",
  "GSVA"
)

BiocManager::install(bioc_packages, ask = FALSE, update = FALSE)
```

## 4. Optional packages

The following packages support vignette building, annotation resources, or
development workflows but are not required for the core installed package:

```r
optional_packages <- c(
  "stringr",
  "BSgenome.Hsapiens.UCSC.hg19",
  "BiocIO",
  "Biostrings",
  "XVector",
  "BiocGenerics",
  "knitr",
  "rmarkdown"
)

BiocManager::install(optional_packages, ask = FALSE, update = FALSE)
```

## 5. Verify the required dependencies

This check reports missing packages but does not install anything:

```r
required_packages <- c(cran_packages, bioc_packages)
installed <- vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)

if (any(!installed)) {
  stop(
    "Install these missing APEX dependencies first: ",
    paste(required_packages[!installed], collapse = ", ")
  )
}

message("All required APEX dependencies are installed.")
```

## 6. Install APEX without dependency resolution

```r
remotes::install_github(
  "Baca-Lab/APEX",
  build_vignettes = FALSE,
  dependencies = FALSE,
  upgrade = "never"
)
```

If a required dependency is absent or too old, APEX installation will fail and
identify the dependency. APEX will not download or update it automatically.
