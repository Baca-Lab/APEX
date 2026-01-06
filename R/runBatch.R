#' Run APEX across multiple samples from a manifest
#'
#' Applies \code{apex()} to each row of a sample manifest and combines
#' predictions into a single gene × sample matrix.
#'
#' Manifest format (data.frame / tibble):
#'   - col 1: sample name
#'   - col 2: H3K4me3 fragment file path (required)
#'   - col 3: H3K27ac fragment file path (optional; NA if not used)
#'   - col 4: H3K36me3 fragment file path (optional; NA if not used)
#'
#' @param manifest data.frame with sample metadata.
#' @param sample_col Column index or name for sample IDs. Default 1.
#' @param k4_col Column index or name for H3K4me3 fragment paths. Default 2.
#' @param k27_col Column index or name for H3K27ac fragment paths. Default 3.
#' @param k36_col Column index or name for H3K36me3 fragment paths. Default 4.
#' @param intersect_genes Logical; intersect genes across samples (default TRUE).
#' @param drop_failed Logical; skip failed samples (default TRUE).
#' @param ... Additional arguments passed to \code{apex()}.
#'
#' @return Numeric matrix (genes × samples) of APEX-inferred expression.
#' @export
apex_batch <- function(
    manifest,
    sample_col = 1,
    k4_col = 2,
    k27_col = 3,
    k36_col = 4,
    intersect_genes = TRUE,
    drop_failed = TRUE,
    ...
) {

  if (!is.data.frame(manifest)) {
    stop("`manifest` must be a data.frame or tibble.")
  }

  get_col <- function(df, col) {
    if (is.numeric(col)) return(df[[col]])
    if (is.character(col)) return(df[[col]])
    stop("Column specifiers must be numeric indices or column names.")
  }

  samples <- as.character(get_col(manifest, sample_col))
  k4 <- get_col(manifest, k4_col)
  k27 <- get_col(manifest, k27_col)
  k36 <- get_col(manifest, k36_col)

  clean_path <- function(x) {
    if (is.na(x) || is.null(x) || trimws(x) == "") return(NULL)
    as.character(x)
  }

  preds <- list()

  for (i in seq_len(nrow(manifest))) {

    s <- samples[i]

    res <- tryCatch(
      {
        apex(
          frag_file_k4  = clean_path(k4[i]),
          frag_file_k27 = clean_path(k27[i]),
          frag_file_k36 = clean_path(k36[i]),
          ...
        )
      },
      error = function(e) {
        if (drop_failed) {
          warning("APEX failed for sample '", s, "': ", conditionMessage(e),
                  call. = FALSE)
          return(NULL)
        }
        stop(e)
      }
    )

    if (is.null(res)) next
    if (!is.matrix(res)) next

    # normalize to named vector
    vec <- res[, 1]
    names(vec) <- rownames(res)

    preds[[s]] <- vec
  }

  if (length(preds) == 0) {
    stop("No samples were successfully processed.")
  }

  gene_sets <- lapply(preds, names)

  if (intersect_genes) {
    genes <- Reduce(intersect, gene_sets)
    if (length(genes) == 0) {
      stop("No overlapping genes across samples.")
    }
    genes <- sort(genes)
    mat <- sapply(preds, function(v) v[genes])
  } else {
    genes <- sort(unique(unlist(gene_sets)))
    mat <- sapply(preds, function(v) {
      out <- rep(NA_real_, length(genes))
      names(out) <- genes
      out[names(v)] <- v
      out
    })
  }

  rownames(mat) <- genes
  colnames(mat) <- names(preds)

  mat
}
