#' Volcano plot for differential expression results
#'
#' Generates a volcano plot from limma-style differential expression output,
#' with optional interactive (plotly) support.
#'
#' @param df Data.frame returned by \code{apex_diff()}.
#' @param highlight_genes Optional character vector of gene names to highlight.
#' @param highlight_col Color for points to highlight.
#' @param signif_col Color for significant points.
#' @param nonsignif_col Color for non-significant points.
#' @param p_cutoff P-value threshold.
#' @param logFC_cutoff Absolute log2 fold-change threshold.
#' @param title Plot title.
#' @param interactive Logical; return an interactive plotly object. Default FALSE.
#'
#' @return A ggplot object (static) or plotly object (interactive).
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_vline labs theme coord_cartesian
#' @importFrom ggrepel geom_text_repel
#' @export
apex_volcano_plot <- function(
    df,
    highlight_genes = NULL,
    highlight_col = "#e85c47",
    signif_col = "#f4a683",
    nonsignif_col = "grey80",
    p_cutoff = 0.05,
    logFC_cutoff = 1,
    title = "",
    interactive = FALSE
) {

  var1 <- unlist(lapply(strsplit(colnames(df)[2], "_"), function(x) x[[2]]))
  var2 <- unlist(lapply(strsplit(colnames(df)[2], "_"), function(x) x[[3]]))
  colnames(df)[2] <- "logFC"

  df <- df %>%
    dplyr::mutate(
      negLogP = -log10(P.Value),
      signif  = ifelse(
        P.Value < p_cutoff & abs(logFC) > logFC_cutoff,
        "Significant",
        "Not Significant"
      ),
      highlight = if (!is.null(highlight_genes)) {
        ifelse(gene %in% highlight_genes, gene, NA)
      } else {
        NA
      }
    )

  x_pad <- diff(range(df$logFC,   na.rm = TRUE)) * 0.12
  y_pad <- diff(range(df$negLogP, na.rm = TRUE)) * 0.12

  # ---- Axis labels (ggplot vs plotly safe) ----
  if (interactive) {
    x_lab <- paste0("log2 fold change ", var1, " - ", var2)
    y_lab <- "-log10(P)"
  } else {
    x_lab <- bquote(log[2] * " fold change " * .(var1) * " - " * .(var2))
    y_lab <- bquote(-log[10] * "(P)")
  }

  # ---- Base ggplot ----
  p <- ggplot(
    df,
    aes(
      x = logFC,
      y = negLogP,
      text = paste0(
        "Gene: ", gene,
        "<br>log2FC: ", round(logFC, 3),
        "<br>P-value: ", signif(P.Value, 3)
      )
    )
  ) +
    geom_point(aes(color = signif), alpha = 0.65, size = 1.6) +

    geom_point(
      data = dplyr::filter(df, !is.na(highlight)),
      shape = 21,
      fill  = highlight_col,
      color = "black",
      size  = 3
    ) +

    { if (!interactive)
      ggrepel::geom_text_repel(
        data = dplyr::filter(df, !is.na(highlight)),
        aes(label = highlight),
        size = 4,
        fontface = "bold",
        box.padding   = 0.8,
        point.padding = 0.3,
        min.segment.length = 0,
        segment.color = "black",
        segment.size  = 0.3
      )
    } +

    geom_hline(
      yintercept = -log10(p_cutoff),
      linetype = "dashed",
      color = "grey50"
    ) +
    geom_vline(
      xintercept = c(-logFC_cutoff, logFC_cutoff),
      linetype = "dashed",
      color = "grey50"
    ) +

    scale_color_manual(values = c(
      "Not Significant" = nonsignif_col,
      "Significant"     = signif_col
    )) +

    labs(
      title = title,
      x = x_lab,
      y = y_lab,
      color = NULL
    ) +

    ggpubr::theme_pubr() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90")
    ) +

    coord_cartesian(
      xlim = c(min(df$logFC,   na.rm = TRUE) - x_pad,
               max(df$logFC,   na.rm = TRUE) + x_pad),
      ylim = c(0, max(df$negLogP, na.rm = TRUE) + y_pad),
      clip = "off"
    )

  # ---- Return interactive or static ----
  if (interactive) {
    return(
      plotly::ggplotly(p, tooltip = "text") %>%
        plotly::layout(
          hoverlabel = list(bgcolor = "white"),
          xaxis = list(title = list(text = x_lab)),
          yaxis = list(title = list(text = y_lab))
        )
    )
  }

  return(p)
}
