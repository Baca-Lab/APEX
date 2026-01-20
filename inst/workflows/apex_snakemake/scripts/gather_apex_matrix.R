#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
per_sample_dir <- args[1]
out_file <- args[2]

files <- list.files(per_sample_dir, pattern="\\.apex\\.tsv$", full.names=TRUE)
stopifnot(length(files) > 0)

read_one <- function(f) {
  x <- read.delim(f, check.names = FALSE)
  sid <- sub("\\.apex\\.tsv$", "", basename(f))
  out <- cbind.data.frame(rownames(x), unlist(x))
  colnames(out) <- c("gene", sid)
  out
}

mat <- Reduce(function(a, b) merge(a, b, by="gene", all=TRUE),
              lapply(files, read_one))

write.table(mat, out_file, sep="\t", quote=FALSE, row.names=FALSE)
