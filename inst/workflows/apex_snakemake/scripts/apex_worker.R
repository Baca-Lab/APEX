#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(apex)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: apex_worker.R sample_id frag_k4 frag_k27 frag_k36 out_tsv")
}

sample_id <- args[1]
frag_k4   <- ifelse(args[2] == "NA", NA, args[2])
frag_k27  <- ifelse(args[3] == "NA", NA, args[3])
frag_k36  <- ifelse(args[4] == "NA", NA, args[4])
out_tsv   <- args[5]

manifest <- data.frame(
  sample_id = sample_id,
  H3K4me3   = frag_k4,
  H3K27ac  = frag_k27,
  H3K36me3 = frag_k36,
  stringsAsFactors = FALSE
)

res <- apex_batch(manifest = manifest, fastMode = FALSE)

write.table(
  res,
  file = out_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = TRUE
)
