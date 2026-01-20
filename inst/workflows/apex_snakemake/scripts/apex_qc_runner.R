#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(apex)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: apex_qc_runner.R manifest.tsv qc_out.tsv")
}

manifest_file <- args[1]
qc_out <- args[2]

manifest <- read.delim(manifest_file, stringsAsFactors = FALSE)

qc <- apex_qc(manifest, plotQC = FALSE)

# Merge enrichment + fragment number into one table
qc_df <- merge(
  qc$enrichment,
  qc$frag_num,
  by = "SampleID",
  all = TRUE
)

write.table(
  qc_df,
  qc_out,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
