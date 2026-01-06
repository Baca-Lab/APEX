#' Quality control for histone mark fragment enrichment
#'
#' This function assesses cfChIP-seq or ChIP-seq library quality by quantifying
#' histone mark-specific enrichment of recovered fragments within predefined
#' on-target and off-target genomic regions. It accepts either (i) an APEX-compatible
#' BED file path or (ii) a manifest (as produced by \code{apex_batch}) containing
#' per-sample fragment file paths for multiple histone marks.
#'
#' @param frag_file A \code{character} string specifying the path to an
#'   APEX-compatible BED file containing fragment information.
#'
#'   The file should be tab-delimited with the following columns:
#'   \itemize{
#'     \item \code{chr} – chromosome name
#'     \item \code{start} – fragment start position (0-based)
#'     \item \code{end} – fragment end position
#'     \item \code{strand} – strand information ("+" or "-" or "*")
#'     \item \code{frag_length} – fragment length
#'     \item \code{gc_bias} – GC content bias estimate
#'     \item \code{End_motif_PosStrand} – 5' 4-bp end motif on positive strand
#'     \item \code{End_motif_NegStrand} – 5' 4-bp end motif on negative strand
#'   }
#'
#'   Only the genomic coordinates (\code{chr}, \code{start}, and \code{end})
#'   are used for enrichment calculations.
#'
#' @param histone_mark A \code{character} string specifying the histone modification
#'   of interest (default: \code{"H3K4me3"}). Supported options include
#'   \code{"H3K4me3"}, \code{"H3K27ac"}, and \code{"H3K36me3"}.
#'
#' @param manifest Optional \code{data.frame} / tibble with sample-level fragment paths.
#'   Manifest format:
#'   \itemize{
#'     \item col 1: sample name
#'     \item col 2: H3K4me3 fragment file path (required by default)
#'     \item col 3: H3K27ac fragment file path (optional; \code{NA} if not used)
#'     \item col 4: H3K36me3 fragment file path (optional; \code{NA} if not used)
#'   }
#'   If \code{manifest} is provided, \code{frag_file} is ignored and the function
#'   returns per-sample QC for all three marks when available.
#'
#' @details
#' This function quantifies two primary quality metrics:
#' \enumerate{
#'   \item \strong{Number of unique fragments} – the total number of fragments analyzed.
#'   \item \strong{Enrichment score} – a measure of histone mark specificity,
#'     defined as the ratio of normalized fragment densities between on-target
#'     and off-target regions:
#'     \deqn{Enrichment = (on\_reads / on\_bp) / (off\_reads / off\_bp)}
#' }
#'
#' On-target and off-target regions are predefined based on chromatin states from
#' the 18-state ChromHMM model of the Roadmap Epigenomics Project
#' (\url{https://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html#exp_18state}).
#' On/off target BED files are loaded from \code{inst/resources/histone_target_sites/<mark>/}.
#'
#' \strong{Recommended quality thresholds:}
#' \itemize{
#'   \item \code{H3K4me3}: enrichment score > 7 and > 1,000,000 fragments
#'   \item \code{H3K27ac}: enrichment score > 2 and > 1,000,000 fragments
#'   \item \code{H3K36me3}: enrichment score > 2 and > 2,000,000 fragments
#' }
#'
#' @return If \code{manifest} is \code{NULL}, returns a named \code{list} with:
#' \itemize{
#'   \item \code{histone_mark}
#'   \item \code{enrichment}
#'   \item \code{frag_num}
#' }
#'
#' If \code{manifest} is provided, returns a \code{data.frame} with one row per sample and
#' columns:
#' \itemize{
#'   \item \code{sample}
#'   \item \code{H3K4me3_enrichment}, \code{H3K4me3_frag_num}
#'   \item \code{H3K27ac_enrichment}, \code{H3K27ac_frag_num}
#'   \item \code{H3K36me3_enrichment}, \code{H3K36me3_frag_num}
#' }
#'
#' @importFrom rtracklayer import
#' @importFrom GenomicRanges GRanges countOverlaps start end
#' @importFrom IRanges IRanges
#' @importFrom utils read.table
#'
#' @export

apex_qc <- function(frag_file = NULL, histone_mark = "H3K4me3", manifest = NULL) {

  .qc_one <- function(frag_file, histone_mark) {

    if (is.null(frag_file) || is.na(frag_file) || !nzchar(frag_file)) {
      return(list(histone_mark = histone_mark, enrichment = NA_real_, frag_num = NA_integer_))
    }

    # Read in APEX-compatible BED file
    frags_df <- utils::read.table(
      frag_file,
      header = FALSE,
      sep = "\t",
      col.names = c(
        "chr", "start", "end", "strand",
        "frag_length", "gc_bias",
        "End_motif_PosStrand", "End_motif_NegStrand"
      )
    )

    # Convert to GRanges
    frags <- GenomicRanges::GRanges(
      seqnames = frags_df$chr,
      ranges = IRanges::IRanges(start = frags_df$start, end = frags_df$end),
      strand = frags_df$strand
    )

    # Construct file paths for target regions
    dir <- "resources/histone_target_sites"
    on_file  <- system.file(paste0(dir, "/", histone_mark, "/on.target.filt.bed"),  package = "apex")
    off_file <- system.file(paste0(dir, "/", histone_mark, "/off.target.filt.bed"), package = "apex")

    if (on_file == "" || off_file == "") {
      stop("Target BED files not found for histone_mark = ", histone_mark,
           ". Expected files under inst/resources/histone_target_sites/", histone_mark, "/")
    }

    # Import on/off target regions
    on  <- rtracklayer::import(on_file)
    off <- rtracklayer::import(off_file)

    # Compute base pair coverage
    on_bp  <- sum(GenomicRanges::end(on)  - GenomicRanges::start(on))
    off_bp <- sum(GenomicRanges::end(off) - GenomicRanges::start(off))

    # Compute read overlaps
    on_reads  <- sum(GenomicRanges::countOverlaps(frags, on,  ignore.strand = TRUE))
    off_reads <- sum(GenomicRanges::countOverlaps(frags, off, ignore.strand = TRUE))

    # Calculate enrichment
    enrichment <- (on_reads / on_bp) / (off_reads / off_bp)
    frag_num <- length(frags)

    list(
      histone_mark = histone_mark,
      enrichment = enrichment,
      frag_num = frag_num
    )
  }

  # -----------------------------
  # Manifest mode (multi-mark QC)
  # -----------------------------
  if (!is.null(manifest)) {

    if (!is.data.frame(manifest)) {
      stop("`manifest` must be a data.frame / tibble.")
    }
    if (ncol(manifest) < 4) {
      stop("`manifest` must have 4 columns: sample, H3K4me3, H3K27ac, H3K36me3.")
    }

    sample_names <- manifest[[1]]
    k4_files  <- manifest[[2]]
    k27_files <- manifest[[3]]
    k36_files <- manifest[[4]]

    out <- data.frame(
      sample = sample_names,
      H3K4me3_enrichment  = NA_real_,
      H3K4me3_frag_num    = NA_integer_,
      H3K27ac_enrichment  = NA_real_,
      H3K27ac_frag_num    = NA_integer_,
      H3K36me3_enrichment = NA_real_,
      H3K36me3_frag_num   = NA_integer_,
      stringsAsFactors = FALSE
    )

    for (i in seq_len(nrow(manifest))) {
      q4  <- .qc_one(k4_files[i],  "H3K4me3")
      q27 <- .qc_one(k27_files[i], "H3K27ac")
      q36 <- .qc_one(k36_files[i], "H3K36me3")

      out$H3K4me3_enrichment[i]  <- q4$enrichment
      out$H3K4me3_frag_num[i]    <- q4$frag_num
      out$H3K27ac_enrichment[i]  <- q27$enrichment
      out$H3K27ac_frag_num[i]    <- q27$frag_num
      out$H3K36me3_enrichment[i] <- q36$enrichment
      out$H3K36me3_frag_num[i]   <- q36$frag_num
    }

    return(out)
  }

  # -----------------------------
  # Single-file mode (backward compatible)
  # -----------------------------
  if (is.null(frag_file)) {
    stop("Provide `frag_file` (single-file mode) or `manifest` (multi-sample mode).")
  }

  histone_mark <- as.character(histone_mark)
  if (!histone_mark %in% c("H3K4me3", "H3K27ac", "H3K36me3")) {
    stop("Unsupported `histone_mark`: ", histone_mark,
         ". Supported: H3K4me3, H3K27ac, H3K36me3.")
  }

  return(.qc_one(frag_file, histone_mark))
}
