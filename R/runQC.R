#' Compute QC metrics for a manifest of cfChIP-seq samples
#'
#' Given a manifest containing fragment file paths for one or more histone marks,
#' computes enrichment scores and fragment counts for each sample and mark.
#'
#' @param manifest data.frame describing the cfChIP-seq samples. The first column
#'   should contain sample identifiers. Subsequent columns should be named
#'   \code{H3K4me3}, \code{H3K27ac}, and/or \code{H3K36me3}, and contain paths to
#'   fragment files for each histone mark (or \code{NA} if a mark was not assayed).
#'
#' @param plotQC logical. If TRUE (default), generates QC boxplots for enrichment
#'   scores and fragment counts using mark-specific thresholds.
#'
#' @return A list with two data.frames:
#'   \describe{
#'     \item{enrichment}{Enrichment scores (rows = samples, columns = histone marks).}
#'     \item{frag_num}{Fragment counts (rows = samples, columns = histone marks).}
#'   }
#'
#' @export

apex_qc <- function(manifest, plotQC = TRUE) {

    stopifnot(
    is.data.frame(manifest),
    ncol(manifest) >= 2
  )

  # Expected structure
  sample_col <- manifest[[1]]
  histone_cols <- c("H3K4me3", "H3K27ac", "H3K36me3")

  # Initialize result matrices
  enrichment_mat <- matrix(
    NA_real_,
    nrow = nrow(manifest),
    ncol = length(histone_cols),
    dimnames = list(sample_col, histone_cols)
  )

  fragnum_mat <- matrix(
    NA_real_,
    nrow = nrow(manifest),
    ncol = length(histone_cols),
    dimnames = list(sample_col, histone_cols)
  )

  # Iterate over samples and histone marks
  for (i in seq_len(nrow(manifest))) {
    for (h in 1:length(histone_cols)) {

      frag_file <- manifest[,h+1][i]

      if (is.na(frag_file) || is.null(frag_file)) {
        next
      }

      qc <- qualityControl(
        frags = frag_file,
        histone_mark = histone_cols[h]
      )

      enrichment_mat[i, h] <- qc$enrichment
      fragnum_mat[i, h]    <- qc$frag_num
    }
  }

  # Convert to data.frames
  enrichment_df <- data.frame(
    sample = sample_col,
    enrichment_mat,
    row.names = NULL,
    check.names = FALSE
  )
  colnames(enrichment_df) <- c("SampleID", paste0("EnrichmentScore_", histone_cols))

  fragnum_df <- data.frame(
    sample = sample_col,
    fragnum_mat,
    row.names = NULL,
    check.names = FALSE
  )
  colnames(fragnum_df) <- c("SampleID", paste0("FragNum_", histone_cols))

if(plotQC){
  p_enrichment <- plot_apex_qc(
    qc_df = enrichment_df,
    value_prefix = "EnrichmentScore",
    ylab = "Enrichment score",
    qc_thresholds = c(
      H3K36me3 = 2,
      H3K4me3  = 7,
      H3K27ac  = 7
    )
  )

  p_enrichment

  p_fragnum <- plot_apex_qc(
    qc_df = fragnum_df,
    value_prefix = "FragNum",
    ylab = "Fragment count",
    qc_thresholds = c(
      H3K4me3  = 1e6,
      H3K36me3 = 2e6,
      H3K27ac  = 1e6
    )
  ) +
    ggplot2::scale_y_log10()

  p_fragnum

  }

return(list(
  enrichment = enrichment_df,
  frag_num   = fragnum_df
))
}
