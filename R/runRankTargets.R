#' Rank cancer targets by APEX-inferred expression percentile
#'
#' Computes gene-wise expression percentiles for therapeutic cancer targets
#' in user-supplied APEX-inferred expression data relative to a curated
#' pan-cancer or cancer-specific APEX reference cohort.
#'
#' The reference cohort contains 573 plasma samples with summarized cancer
#' annotations stored as a character vector. Supported
#' cancer categories include:
#' \itemize{
#'   \item \code{"Breast"}
#'   \item \code{"CNS"}
#'   \item \code{"Colorectal"}
#'   \item \code{"Gynecologic"}
#'   \item \code{"Kidney"}
#'   \item \code{"Liver"}
#'   \item \code{"Melanoma"}
#'   \item \code{"NSCLC"}
#'   \item \code{"Pancreatobiliary"}
#'   \item \code{"Prostate / NEPC"}
#'   \item \code{"SCLC / Neuroendocrine"}
#'   \item \code{"Thymoma / Mesothelioma"}
#'   \item \code{"Upper GI"}
#'   \item \code{"Urothelial"}
#'   \item \code{"Other"}
#' }
#'
#' Sample counts per cancer type can be queried programmatically using
#' \code{table(ref$pheno)}.
#'
#' @param apex_mat Numeric matrix (genes x samples) of APEX-inferred expression.
#'        Row names must be gene symbols.
#' @param top_n Integer; number of top-ranked cancer targets to return per sample.
#' @param ref_cancer Optional character; cancer type used to subset the reference
#'        cohort (must be one of the categories listed above). If \code{NULL},
#'        percentiles are computed relative to the full pan-cancer reference.
#' @param cancer_col Character; retained for backward compatibility but ignored.
#'        Cancer labels are stored as a character vector in \code{ref$pheno}.
#' @param plot Logical; whether to generate a simple bar plot of top-ranked targets.
#'
#' @return A list containing:
#' \describe{
#'   \item{percentiles}{Numeric matrix of gene-wise percentiles (genes x samples).}
#'   \item{top_targets}{List of data.frames with top-ranked targets per sample.}
#'   \item{reference_n}{Table of sample counts per cancer type used as reference.}
#'   \item{reference_set}{Character string indicating pan-cancer or cancer-specific reference.}
#'   \item{plot}{ggplot object (if \code{plot = TRUE}).}
#' }
#'
#' @importFrom utils read.table
#' @importFrom ggplot2 scale_color_manual element_text element_blank element_line
#' @export
apex_rank_targets <- function(
    apex_mat,
    top_n = 3,
    ref_cancer = NULL,
    cancer_col = "cancer_type",
    plot = TRUE
) {

  # -----------------------------
  # Load packaged resources
  # -----------------------------
  targets <- read.table(
    system.file("resources/cancer_target_genes.txt", package = "apex"),
    stringsAsFactors = FALSE, header = T
  )[, 1]

  ref <- readRDS(
    system.file("resources/APEX_inferred_GEXP_573samples.rds", package = "apex")
  )

  ref_mat       <- ref$mat
  cancer_labels <- ref$pheno

  # -----------------------------
  # Optional cancer-specific reference
  # -----------------------------
  if (!is.null(ref_cancer)) {
    keep <- cancer_labels == ref_cancer
    if (sum(keep, na.rm = TRUE) == 0) {
      stop("No reference samples found for cancer type: ", ref_cancer)
    }
    ref_mat       <- ref_mat[, keep, drop = FALSE]
    cancer_labels <- cancer_labels[keep]
  }

  # reference sample counts (for the reference actually used)
  ref_n <- table(cancer_labels)

  # -----------------------------
  # Restrict to shared target genes
  # -----------------------------
  genes <- intersect(
    targets,
    intersect(rownames(apex_mat), rownames(ref_mat))
  )

  if (length(genes) == 0) {
    stop("No overlapping cancer target genes found.")
  }

  apex_mat <- apex_mat[genes, , drop = FALSE]
  ref_mat  <- ref_mat[genes, , drop = FALSE]

  # -----------------------------
  # Percentile transform (gene-wise)
  # -----------------------------
  pct <- sapply(seq_len(ncol(apex_mat)), function(j) {
    sapply(seq_len(nrow(apex_mat)), function(i) {
      mean(ref_mat[i, ] <= apex_mat[i, j], na.rm = TRUE) * 100
    })
  })

  rownames(pct) <- genes
  colnames(pct) <- colnames(apex_mat)

  # -----------------------------
  # Top targets per sample
  # -----------------------------
  top_targets <- lapply(seq_len(ncol(pct)), function(j) {
    ord <- order(pct[, j], decreasing = TRUE)
    k   <- min(top_n, length(ord))
    data.frame(
      gene       = rownames(pct)[ord][seq_len(k)],
      percentile = pct[ord, j][seq_len(k)],
      sample     = colnames(pct)[j],
      rank       = seq_len(k),
      row.names  = NULL
    )
  })
  names(top_targets) <- colnames(pct)

  # -----------------------------
  # Publication-quality plot
  # -----------------------------
  p <- NULL
  if (isTRUE(plot)) {

    df <- do.call(rbind, top_targets)

    n_samples   <- length(unique(df$sample))
    ncol_facets <- ceiling(sqrt(n_samples))

    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(
        x = percentile,
        y = reorder(gene, percentile),
        fill = rank == 1
      )
    ) +
      ggplot2::geom_col(
        width = 0.75,
        color = "black",
        linewidth = 0.25
      ) +
      ggplot2::scale_fill_manual(
        values = c("TRUE" = "#e85c47", "FALSE" = "grey80"),
        guide = "none"
      ) +
      ggplot2::facet_wrap(
        ~sample,
        ncol = ncol_facets,
        scales = "free_y"
      ) +
      ggplot2::scale_x_continuous(
        limits = c(0, 100),
        expand = c(0, 0)
      ) +
      ggplot2::labs(
        x = "APEX expression percentile (reference cohort)",
        y = NULL,
        title = "Top cancer targets by APEX-inferred expression"
      ) +
      ggplot2::theme_classic(base_size = 14) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          hjust = 0.5,
          face = "bold",
          size = 16
        ),
        strip.background = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(
          face = "bold",
          size = 13
        ),
        panel.border = ggplot2::element_rect(
          color = "black",
          fill = NA,
          linewidth = 0.6
        ),
        axis.line = ggplot2::element_line(color = "black"),
        axis.ticks = ggplot2::element_line(color = "black"),
        panel.spacing = grid::unit(0.6, "lines")
      )
  }

  # -----------------------------
  # Return
  # -----------------------------
  list(
    percentiles   = pct,
    top_targets   = top_targets,
    reference_n   = ref_n,
    reference_set = ifelse(is.null(ref_cancer), "Pan-cancer", ref_cancer),
    plot          = p
  )
}
