#' Compute QC metrics for a manifest of cfChIP-seq samples
#'
#' Given a manifest with fragment file paths per histone mark, computes
#' enrichment scores and fragment counts for each sample and mark.
#'
#' @param manifest data.frame with columns:
#'   sample, H3K4me3, H3K27ac, H3K36me3 (file paths or NA)
#'
#' @return A list with two data.frames:
#'   \describe{
#'     \item{enrichment}{Enrichment scores (rows = samples, cols = histone marks)}
#'     \item{frag_num}{Fragment counts (rows = samples, cols = histone marks)}
#'   }
#'
#' @export
apex_qc <- function(manifest) {

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

  return(list(
    enrichment = enrichment_df,
    frag_num   = fragnum_df
  ))
}
