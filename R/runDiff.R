#' Differential expression analysis of APEX-inferred expression using limma
#'
#' Performs two-group differential expression analysis on a matrix of APEX-inferred
#' gene expression values using limma with an automatically defined contrast.
#'
#' @param apex_mat Numeric matrix of APEX-inferred expression values (genes x samples).
#'   Row names must correspond to gene identifiers.
#' @param group Factor or character vector defining sample groups
#'   (length must equal ncol(apex_mat)); must have exactly two levels.
#'
#' @return A data.frame with limma statistics for all genes. Attributes include:
#'   \describe{
#'     \item{contrast}{Human-readable contrast label.}
#'     \item{n_per_group}{Sample counts per group (original labels).}
#'   }
#'
#' @importFrom limma lmFit makeContrasts contrasts.fit eBayes topTable
#' @export
apex_diff <- function(apex_mat, group) {

  # ----------------------------
  # Input validation
  # ----------------------------
  if (!is.matrix(apex_mat) || !is.numeric(apex_mat)) {
    stop("`apex_mat` must be a numeric matrix (genes x samples).")
  }

  if (ncol(apex_mat) != length(group)) {
    stop("Length of `group` must equal ncol(apex_mat).")
  }

  if (is.null(rownames(apex_mat))) {
    stop("`apex_mat` must have rownames corresponding to gene identifiers.")
  }

  group <- factor(group)

  if (nlevels(group) != 2) {
    stop("`group` must have exactly two levels.")
  }

  # Preserve original labels for reporting
  group_levels <- levels(group)
  n_per_group  <- table(group)  # uses original labels

  # Make safe names for limma design matrix
  safe_levels <- make.names(group_levels)
  levels(group) <- safe_levels

  # ----------------------------
  # Design + contrast
  # ----------------------------
  design <- model.matrix(~ 0 + group)
  colnames(design) <- safe_levels

  contrast_str <- paste0(safe_levels[2], " - ", safe_levels[1])
  contrast <- limma::makeContrasts(contrasts = contrast_str, levels = design)

  # ----------------------------
  # Fit model
  # ----------------------------
  fit <- limma::lmFit(apex_mat, design)
  fit <- limma::contrasts.fit(fit, contrast)
  fit <- limma::eBayes(fit)

  # ----------------------------
  # Extract results
  # ----------------------------
  res <- limma::topTable(
    fit,
    coef    = 1,
    number  = Inf,
    sort.by = "none"
  )

  res$gene <- rownames(apex_mat)

  # limma topTable returns AveExpr (not Aveapex_mat)
  res <- res[, c(
    "gene",
    "logFC",
    "AveExpr",
    "t",
    "P.Value",
    "adj.P.Val",
    "B"
  )]
colnames(res)[which(colnames(res) %in% "logFC")] <- paste0("logFC_",safe_levels[2], "_", safe_levels[1])
  # Order by effect size
  res <- res[order(-res[,grep("logFC", colnames(res))]), ]
  return(res)
}
