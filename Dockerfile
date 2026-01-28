FROM bioconductor/bioconductor_docker:RELEASE_3_20

## -----------------------------
## 1. System libraries
## -----------------------------
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    libhdf5-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    pandoc \
    && rm -rf /var/lib/apt/lists/*

## -----------------------------
## 2. CRAN packages (Imports)
## -----------------------------
RUN R -e "install.packages(c( \
    'magrittr', \
    'data.table', \
    'reshape2', \
    'dplyr', \
    'xgboost', \
    'ggplot2', \
    'ggrepel', \
    'ggpubr', \
    'plotly', \
    'msigdbr' \
), repos='https://cloud.r-project.org')"

## -----------------------------
## 3. Bioconductor packages (Imports)
## -----------------------------
RUN R -e "BiocManager::install(c( \
    'GenomicRanges', \
    'IRanges', \
    'GenomeInfoDb', \
    'S4Vectors', \
    'rtracklayer', \
    'limma', \
    'GSVA' \
), ask = FALSE)"

## -----------------------------
## 4. Suggested packages (optional but recommended)
## -----------------------------
RUN R -e "install.packages(c( \
    'stringr', \
    'knitr', \
    'rmarkdown', \
    'remotes' \
), repos='https://cloud.r-project.org')"

RUN R -e "BiocManager::install(c( \
    'BiocIO', \
    'Biostrings', \
    'XVector', \
    'BiocGenerics', \
    'BSgenome.Hsapiens.UCSC.hg19' \
), ask = FALSE)"

## -----------------------------
## 5. Install your package
## -----------------------------
WORKDIR /opt/pkg
COPY . /opt/pkg
RUN R CMD INSTALL /opt/pkg

## -----------------------------
## 6. HPC-safe defaults
## -----------------------------
ENV R_LIBS_USER=""
ENV OMP_NUM_THREADS=1
ENV MKL_NUM_THREADS=1

## -----------------------------
## 7. Default command
## -----------------------------
CMD ["R"]
