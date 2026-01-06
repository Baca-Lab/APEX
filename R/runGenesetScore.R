#' Gene set scoring (ssGSEA or average log expression)
#'
#' Computes per-sample gene set scores from an APEX-inferred expression matrix using either
#' ssGSEA (via GSVA) or the average expression across genes in each set.
#'
#' Users can score (i) an MSigDB collection/subcollection via \pkg{msigdbr} or
#' (ii) custom gene sets supplied as a named list.
#'
#' @param apex_mat Numeric matrix of APEX-inferred expression values (genes x samples).
#'   Row names must be gene symbols.
#' @param method Character. Scoring method: \code{"ssgsea"} or \code{"mean"}.
#' @param geneset_source Character. \code{"msigdb"} or \code{"custom"}.
#' @param msigdb_collection Character. MSigDB collection abbreviation (e.g., \code{"H"}, \code{"C2"}).
#'   Used only when \code{geneset_source = "msigdb"}.
#' @param msigdb_subcollection Optional character. MSigDB subcollection.
#' @param msigdb_species Character. Species name for output genes.
#' @param msigdb_db_species Character. Database species abbreviation (\code{"HS"} or \code{"MM"}).
#' @param genesets Named list of character vectors (gene symbols). Required when \code{geneset_source = "custom"}.
#' @param geneset_names Optional character vector of gene set names to score.
#' @param min_genes Integer. Minimum number of overlapping genes required to compute a score.
#' @param verbose Logical. If \code{TRUE}, prints brief progress messages.
#' @param ... Additional arguments passed to the underlying method.
#'
#' @return A numeric matrix of gene set scores (gene sets x samples).
#'   Row names are gene set names; column names match \code{colnames(apex_mat)}.
#'
#' @importFrom msigdbr msigdbr
#' @importFrom GSVA ssgseaParam gsva
#' @export
#'
apex_geneset_score <- function(
    apex_mat,
    method = c("ssgsea", "mean"),
    geneset_source = c("msigdb", "custom"),
    msigdb_collection = "H",
    msigdb_subcollection = NULL,
    msigdb_species = "Homo sapiens",
    msigdb_db_species = "HS",
    genesets = NULL,
    geneset_names = NULL,
    min_genes = 5L,
    verbose = TRUE,
    ...
) {

  .msg <- function(...) if (isTRUE(verbose)) message(...)

  # -----------------------------
  # Validate input matrix
  # -----------------------------
  method <- match.arg(method)
  geneset_source <- match.arg(geneset_source)

  if (!is.matrix(apex_mat)) {
    stop("`apex_mat` must be a numeric matrix (genes x samples).")
  }
  if (!is.numeric(apex_mat)) {
    stop("`apex_mat` must be numeric.")
  }
  if (is.null(rownames(apex_mat))) {
    stop("`apex_mat` must have rownames containing gene symbols.")
  }
  if (is.null(colnames(apex_mat))) {
    stop("`apex_mat` must have colnames (sample IDs).")
  }

  # -----------------------------
  # Get gene sets
  # -----------------------------
  gs <- NULL

  if (geneset_source == "msigdb") {
    .msg("Loading MSigDB gene sets via msigdbr...")

    msig <- msigdbr::msigdbr(
      db_species = msigdb_db_species,
      species = msigdb_species,
      collection = msigdb_collection,
      subcollection = msigdb_subcollection
    )

    if (!all(c("gs_name", "gene_symbol") %in% colnames(msig))) {
      stop("Unexpected msigdbr output; expected columns `gs_name` and `gene_symbol`.")
    }

    gs <- split(msig$gene_symbol, msig$gs_name)

  } else {
    if (is.null(genesets)) {
      stop("When `geneset_source = 'custom'`, you must provide `genesets` as a named list.")
    }
    if (!is.list(genesets) || is.null(names(genesets)) || any(names(genesets) == "")) {
      stop("`genesets` must be a *named* list of gene vectors.")
    }
    gs <- genesets
  }

  # Optional subset
  if (!is.null(geneset_names)) {
    missing_names <- setdiff(geneset_names, names(gs))
    if (length(missing_names) > 0) {
      stop("Requested gene sets not found: ", paste(missing_names, collapse = ", "))
    }
    gs <- gs[geneset_names]
  }

  # Apply minimum overlap filter
  overlap_n <- vapply(gs, function(x) sum(rownames(apex_mat) %in% x), integer(1))
  keep <- overlap_n >= as.integer(min_genes)
  dropped <- names(gs)[!keep]
  gs <- gs[keep]

  if (length(gs) == 0) {
    stop("No gene sets passed the `min_genes` filter.")
  }

  if (length(dropped) > 0) {
    .msg("Dropping ", length(dropped), " gene sets with < ", min_genes, " overlapping genes.")
  }

  .msg("Scoring ", length(gs), " gene sets across ", ncol(apex_mat),
       " samples using method = ", method, ".")

  # -----------------------------
  # Score
  # -----------------------------
  if (method == "ssgsea") {
    y <- GSVA::ssgseaParam(data.matrix(apex_mat), gs, ...)
    out <- GSVA::gsva(y)
    return(as.matrix(data.frame(out)))
  }

  # mean method
  out_list <- lapply(names(gs), function(nm) {
    genes <- intersect(rownames(apex_mat), gs[[nm]])
    colMeans(apex_mat[genes, , drop = FALSE], na.rm = TRUE)
  })

  out <- do.call(rbind, out_list)
  rownames(out) <- names(gs)
  colnames(out) <- colnames(apex_mat)

  return(as.matrix(data.frame(out)))
}
