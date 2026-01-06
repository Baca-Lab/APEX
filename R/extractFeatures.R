#' Extract Features from Fragment Files
#'
#' This function extracts features from fragment files and generates a dataframe
#' with several epigenetic features for each gene.
#'
#' @param frag_file_k4 Path to the H3K4me3 fragment file (required).
#' @param frag_file_k36 Path to the H3K36me3 fragment file (optional).
#' @param frag_file_k27 Path to the H3K27ac fragment file (optional). Only used when \code{fastMode = FALSE}.
#' @param upstream Number of base pairs upstream to include in promoter regions. Default: 3000.
#' @param downstream Number of base pairs downstream to include in promoter regions. Default: 3000.
#' @param num_of_promoter_tiles Number of promoter tiles to generate. Default: 20.
#' @param num_of_genebody_tiles Number of gene body tiles to generate. Default: 20.
#' @param num_of_intergenic_tiles Number of tiles for intergenic enhancer regions. Only used when \code{fastMode = FALSE}. Default: 1.
#' @param num_of_intragenic_tiles Number of tiles for intragenic enhancer regions. Only used when \code{fastMode = FALSE}. Default: 1.
#' @param fastMode Logical; if \code{TRUE} (default), runs the faster model using promoter/gene-body features.
#'        If \code{FALSE}, includes H3K27ac and intergenic/intragenic enhancer features.
#' @param useAltPromoter Logical; whether to use alternate promoters based on highest H3K4me3 coverage. Default: TRUE.
#' @param verbose Logical. If TRUE (default), progress messages are printed during execution. If FALSE, all messages are suppressed.
#'
#' @return A dataframe containing the extracted features.
#' @export
#'
extractFeatures <- function(
    frag_file_k4 = NULL,
    frag_file_k36 = NULL,
    frag_file_k27 = NULL,
    upstream = 3000,
    downstream = 3000,
    num_of_promoter_tiles = 20L,
    num_of_genebody_tiles = 20L,
    num_of_intergenic_tiles = 1L,
    num_of_intragenic_tiles = 1L,
    fastMode = TRUE,
    useAltPromoter = TRUE,
    verbose = TRUE
) {


  ## Welcome banner
  .msg("
----------------------------------------------
       ___      .______    __________   ___
      /   \\     |   _  \\  |   ____\\  \\ /  /
     /  ^  \\    |  |_)  | |  |__   \\  V  /
    /  /_\\  \\   |   ___/  |   __|   >   <
   /  _____  \\  |  |      |  |____ /  .  \\
  /__/     \\__\\ | _|      |_______/__/ \\__\\
----------------------------------------------
Integrating plasma cell-free epigenomic features to predict gene expression...
", verbose = verbose)

  # Load reference sites
  sites <- readRDS(system.file(
    "resources/promoter_and_genebody_enhancer_sites.rds",
    package = "apex"
  ))
  promoter_sites <- sites$promoter_sites
  genebody_sites <- sites$genebody_sites

  if (!fastMode) {
    .msg(
      "fastMode = FALSE: Using full model (H3K27ac + intergenic/intragenic enhancers).",
      verbose = verbose
    )
    intergenic_sites <- sites$intergenic_sites
    intragenic_sites <- sites$intragenic_sites
  }

  quality_k4 <- quality_k36 <- quality_k27 <- NULL

  ## Load fragment files
  frag_orig_k4 <- if (!is.null(frag_file_k4)) {
    .msg("Loading fragment file: ", frag_file_k4, verbose = verbose)
    frag <- readFragBed(frag_file_k4)
    frag <- keepSeqlevels(
      frag,
      intersect(seqlevels(frag), seqlevels(genebody_sites)),
      pruning.mode = "coarse"
    )
    quality_k4 <- qualityControl(frag, histone_mark = "H3K4me3")
    frag
  }

  frag_orig_k36 <- if (!is.null(frag_file_k36)) {
    .msg("Loading fragment file: ", frag_file_k36, verbose = verbose)
    frag <- readFragBed(frag_file_k36)
    frag <- keepSeqlevels(
      frag,
      intersect(seqlevels(frag), seqlevels(genebody_sites)),
      pruning.mode = "coarse"
    )
    quality_k36 <- qualityControl(frag, histone_mark = "H3K36me3")
    frag
  }

  if (!fastMode && !is.null(frag_file_k27)) {
    .msg("Loading fragment file: ", frag_file_k27, verbose = verbose)
    frag <- readFragBed(frag_file_k27)
    frag <- keepSeqlevels(
      frag,
      intersect(seqlevels(frag), seqlevels(genebody_sites)),
      pruning.mode = "coarse"
    )
    quality_k27 <- qualityControl(frag, histone_mark = "H3K27ac")
    frag_orig_k27 <- frag
  } else {
    frag_orig_k27 <- NULL
  }

  ## QC aggregation
  quality_results <- list(
    H3K4me3 = quality_k4,
    H3K36me3 = quality_k36,
    H3K27ac = quality_k27
  )
  quality_results <- qc(quality_results)

  if (is.null(frag_file_k4) &&
      is.null(frag_file_k36) &&
      is.null(frag_file_k27)) {
    .msg(
      "No valid fragment files provided. At least H3K4me3 or H3K36me3 is required.",
      verbose = verbose
    )
  }

  if (is.null(frag_file_k4)) {
    .msg(
      "No H3K4me3 provided; alternate promoter selection disabled.",
      verbose = verbose
    )
    useAltPromoter <- FALSE
  }

  ## Alternate promoter selection
  if (useAltPromoter) {
    .msg(
      "Determining optimal TSS per gene based on H3K4me3 coverage...",
      verbose = verbose
    )

    coverage_by_promoter <- countOverlaps(
      resize(promoter_sites, width = 1001, fix = "center"),
      frag_orig_k4, ignore.strand = TRUE
    )

    max_tss <- unlist(tapply(
      seq_along(coverage_by_promoter),
      promoter_sites$mcols.gene_name,
      function(i) i[which.max(coverage_by_promoter[i])]
    ))

    promoter_sites <- promoter_sites[max_tss]
    tss <- mid(promoter_sites)
    tts <- start(getTTS(genebody_sites))

    promoter_sites <- promoters(
      resize(promoter_sites, width = 1, fix = "center"),
      upstream = upstream,
      downstream = downstream
    )

    genebody_sites <- GRanges(
      seqnames = seqnames(genebody_sites),
      ranges = IRanges(start = pmin(tss, tts), end = pmax(tss, tts)),
      strand = strand(genebody_sites),
      mcols.gene_name = genebody_sites$mcols.gene_name
    )

    keep <- width(genebody_sites) >= 500
    promoter_sites <- promoter_sites[keep]
    genebody_sites <- genebody_sites[keep]
  } else {
    promoter_sites <- promoters(
      getTSS(genebody_sites),
      upstream = upstream,
      downstream = downstream
    )
  }

  .msg(
    paste0(
      "Final reference contains ",
      length(unique(promoter_sites$mcols.gene_name)),
      " genes."
    ),
    verbose = verbose
  )

  ## Tiling
  promoter_sites <- tile_site(promoter_sites, num_of_tiles = num_of_promoter_tiles)
  genebody_sites <- tile_site(genebody_sites, num_of_tiles = num_of_genebody_tiles)

  if (!fastMode) {
    intergenic_sites <- tile_site(
      intergenic_sites[width(intergenic_sites) > num_of_intergenic_tiles],
      num_of_tiles = num_of_intergenic_tiles
    )
    intragenic_sites <- tile_site(
      intragenic_sites[width(intragenic_sites) > num_of_intragenic_tiles],
      num_of_tiles = num_of_intragenic_tiles
    )
  }

  ## Assign unique IDs
  promoter_sites$mcols.uniqueid <- seq_along(promoter_sites)
  genebody_sites$mcols.uniqueid <- seq_along(genebody_sites)

  if (!fastMode) {
    intergenic_sites$mcols.uniqueid <- seq_along(intergenic_sites)
    intragenic_sites$mcols.uniqueid <- seq_along(intragenic_sites)
  }

  ## Feature generation (unchanged logic)
  feature_list <- list()

  if (!is.null(frag_orig_k4)) {
    feature_list <- c(
      feature_list,
      list(
        generateFeatures(promoter_sites, frag_orig_k4, "K4_promoter"),
        generateFeatures(genebody_sites, frag_orig_k4, "K4_genebody")
      )
    )
    if (!fastMode) {
      feature_list <- c(
        feature_list,
        list(
          generateFeatures(intergenic_sites, frag_orig_k4, "K4_intergenic"),
          generateFeatures(intragenic_sites, frag_orig_k4, "K4_intragenic")
        )
      )
    }
  }

  if (!is.null(frag_orig_k36)) {
    feature_list <- c(
      feature_list,
      list(
        generateFeatures(promoter_sites, frag_orig_k36, "K36_promoter"),
        generateFeatures(genebody_sites, frag_orig_k36, "K36_genebody")
      )
    )
    if (!fastMode) {
      feature_list <- c(
        feature_list,
        list(
          generateFeatures(intergenic_sites, frag_orig_k36, "K36_intergenic"),
          generateFeatures(intragenic_sites, frag_orig_k36, "K36_intragenic")
        )
      )
    }
  }

  if (!fastMode && !is.null(frag_orig_k27)) {
    feature_list <- c(
      feature_list,
      list(
        generateFeatures(promoter_sites, frag_orig_k27, "K27_promoter"),
        generateFeatures(genebody_sites, frag_orig_k27, "K27_genebody"),
        generateFeatures(intergenic_sites, frag_orig_k27, "K27_intergenic"),
        generateFeatures(intragenic_sites, frag_orig_k27, "K27_intragenic")
      )
    )
  }

  feature_list <- lapply(
    feature_list,
    function(x) setnames(x, 1, "Gene")
  )

  final_df <- Reduce(
    function(x, y) merge(x, y, by = "Gene", all = TRUE),
    feature_list
  )

  return(final_df)
}
