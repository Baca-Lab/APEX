#' Differential gene set analysis using APEX
#'
#' Computes gene set scores from APEX-inferred expression and performs
#' differential analysis between two groups using limma.
#'
#' @param apex_mat Numeric matrix of APEX-inferred expression (genes x samples).
#' @param group Factor or character vector defining two sample groups.
#' @param ... Passed to apex_geneset_score().
#'
#' @return A data.frame of differential gene set statistics.
#' @export
apex_geneset_diff <- function(apex_mat, group, ...) {

  # 1. score gene sets
  gs_scores <- apex_geneset_score(apex_mat, ...)

  # gs_scores: gene sets x samples
  # apex_diff expects genes x samples — perfect match
  res <- apex_diff(gs_scores, group)

  return(res)
}
