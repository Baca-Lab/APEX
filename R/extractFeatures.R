#' Extract Features from Fragment Files
#'
#' This function extracts epigenomic features from fragment
#' files and returns a dataframe of features for each gene.
#'
#' @param frag_file_k4 Path to the H3K4me3 fragment file (optional). Used for
#'        promoter features and, when \code{useAltPromoter = TRUE}, to select
#'        the transcription start site (TSS) with the highest H3K4me3 coverage.
#' @param frag_file_k36 Path to the H3K36me3 fragment file (optional).
#' @param frag_file_k27 Path to the H3K27ac fragment file (optional). Used for
#'        enhancer features when \code{fastMode = FALSE}, and as the fallback
#'        mark for alternate TSS selection when \code{useAltPromoter = TRUE}
#'        and H3K4me3 is not provided.
#' @param custom_coordinates Optional path to a BED6 file containing one
#’        user-defined genomic interval per gene. BED column 4 must contain the
#’        gene symbol, and BED column 6 must contain ‘+’ or ‘-’ strand
#’        information. These intervals replace the default gene-body
#’        coordinates, and TSSs are derived from their strand-aware boundaries.
#’        When \code{fastMode = FALSE}, intergenic and intragenic enhancer
#’        regions remain those from the package’s default reference annotation
#’        and are not redefined from \code{custom_coordinates}.
#' @param upstream Number of base pairs upstream of the TSS to include in
#'        promoter regions. Default: 3000.
#' @param downstream Number of base pairs downstream of the TSS to include in
#'        promoter regions. Default: 3000.
#' @param num_of_promoter_tiles Number of promoter tiles. Default: 20.
#' @param num_of_genebody_tiles Number of gene-body tiles. Default: 20.
#' @param num_of_intergenic_tiles Number of intergenic enhancer tiles. Only
#'        used when \code{fastMode = FALSE}. Default: 1.
#' @param num_of_intragenic_tiles Number of intragenic enhancer tiles. Only
#'        used when \code{fastMode = FALSE}. Default: 1.
#' @param fastMode Logical; if \code{TRUE} (default), uses promoter and
#'        gene-body features. If \code{FALSE}, additionally includes H3K27ac
#'        and intergenic/intragenic enhancer features.
#' @param useAltPromoter Logical; if \code{TRUE} (default), selects the TSS
#'        with the highest H3K4me3 coverage for each gene, or H3K27ac coverage
#'        if H3K4me3 is unavailable. If \code{FALSE}, user-defined coordinates
#'        from \code{custom_coordinates} are used when provided; otherwise,
#'        default annotated gene coordinates and TSSs are used.
#' @param verbose Logical; if \code{TRUE} (default), progress messages are
#'        printed. If \code{FALSE}, messages are suppressed.
#'
#' @return A dataframe containing the extracted features.
#' @export
#'
extractFeatures <- function(
    frag_file_k4 = NULL, frag_file_k36 = NULL, frag_file_k27 = NULL,
    custom_coordinates = NULL, upstream = 3000, downstream = 3000,
    num_of_promoter_tiles = 20L, num_of_genebody_tiles = 20L,
    num_of_intergenic_tiles = 1L, num_of_intragenic_tiles = 1L,
    fastMode = TRUE, useAltPromoter = TRUE, verbose = TRUE
) {
  
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
  
  if (all(vapply(list(frag_file_k4, frag_file_k36, frag_file_k27), is.null, logical(1))))
    stop("At least one fragment file must be provided.")
  
  ## Load reference sites
  sites <- readRDS(system.file(
    "resources/promoter_and_genebody_enhancer_sites.rds", package = "apex"
  ))
  promoter_sites <- sites$promoter_sites
  genebody_sites <- sites$genebody_sites
  
  if (!fastMode) {
    .msg("fastMode = FALSE: Using full model (H3K27ac + intergenic/intragenic enhancers).",
         verbose = verbose)
    intergenic_sites <- sites$intergenic_sites
    intragenic_sites <- sites$intragenic_sites
  }
  
  ## Custom gene coordinates
  if (!is.null(custom_coordinates)) {
    .msg("Using user-provided gene coordinates...", verbose = verbose)
    genebody_sites <- rtracklayer::import(custom_coordinates)
    
    if (!length(genebody_sites))
      stop("custom_coordinates contains no genomic intervals.")
    
    ## BED column 4 is imported as 'name'; match APEX reference convention
    if ("name" %in% names(mcols(genebody_sites)))
      names(mcols(genebody_sites))[names(mcols(genebody_sites)) == "name"] <- "mcols.gene_name"
    if (!"mcols.gene_name" %in% names(mcols(genebody_sites)))
      stop("BED column 4 must contain a gene identifier.")
    if (any(is.na(genebody_sites$mcols.gene_name) | genebody_sites$mcols.gene_name == ""))
      stop("All custom coordinates must have a non-empty gene identifier.")
    if (any(strand(genebody_sites) == "*"))
      stop("BED column 6 must contain '+' or '-' strand information.")
    if (anyDuplicated(genebody_sites$mcols.gene_name))
      stop("custom_coordinates must contain one interval per gene.")
    
    promoter_sites <- promoters(
      getTSS(genebody_sites), upstream = upstream, downstream = downstream
    )

  
    if (!fastMode)
      .msg("Custom coordinates replace promoter/gene-body regions; default enhancer annotations are retained.",
           verbose = verbose)
  }
  
  ## Disable alternate promoter selection if no suitable promoter mark exists
  if (useAltPromoter && is.null(frag_file_k4) && is.null(frag_file_k27)) {
    .msg("Neither H3K4me3 nor H3K27ac was provided; useAltPromoter has been set to FALSE.",
         verbose = verbose)
    useAltPromoter <- FALSE
  }
  
  quality_k4 <- quality_k36 <- quality_k27 <- NULL
  
  ## Helper for loading fragment files
  load_frag <- function(file, mark) {
    if (is.null(file)) return(NULL)
    .msg("Loading fragment file: ", file, verbose = verbose)
    frag <- readFragBed(file)
    frag <- keepSeqlevels(
      frag, intersect(seqlevels(frag), seqlevels(genebody_sites)),
      pruning.mode = "coarse"
    )
    assign(paste0("quality_", switch(mark,
                                     H3K4me3 = "k4", H3K36me3 = "k36", H3K27ac = "k27")),
           qualityControl(frag, histone_mark = mark), envir = parent.frame())
    frag
  }
  
  frag_orig_k4 <- load_frag(frag_file_k4, "H3K4me3")
  frag_orig_k36 <- load_frag(frag_file_k36, "H3K36me3")
  frag_orig_k27 <- if (!is.null(frag_file_k27) && (!fastMode || useAltPromoter))
    load_frag(frag_file_k27, "H3K27ac") else NULL
  
  quality_results <- qc(list(
    H3K4me3 = quality_k4, H3K36me3 = quality_k36, H3K27ac = quality_k27
  ), verbose = verbose)
  
  ## Promoter selection
  if (useAltPromoter) {
    if (!is.null(frag_orig_k4)) {
      .msg("Determining optimal TSS per gene based on H3K4me3 coverage...",
           verbose = verbose)
      frag_promoter <- frag_orig_k4
    } else {
      .msg("H3K4me3 was not provided. Determining optimal TSS per gene based on H3K27ac coverage...",
           verbose = verbose)
      frag_promoter <- frag_orig_k27
    }
    
    coverage_by_promoter <- countOverlaps(
      resize(promoter_sites, width = 1001, fix = "center"),
      frag_promoter, ignore.strand = TRUE
    )
    
    max_tss <- unlist(tapply(
      seq_along(coverage_by_promoter), promoter_sites$mcols.gene_name,
      function(i) i[which.max(coverage_by_promoter[i])]
    ))
    
    promoter_sites <- promoter_sites[max_tss]
    tss <- mid(promoter_sites)
    tts <- start(getTTS(genebody_sites))
    
    promoter_sites <- promoters(
      resize(promoter_sites, width = 1, fix = "center"),
      upstream = upstream, downstream = downstream
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
    
  } else if (is.null(custom_coordinates)) {
    .msg("Using default annotated transcription start sites...", verbose = verbose)
    promoter_sites <- promoters(
      getTSS(genebody_sites), upstream = upstream, downstream = downstream
    )
  }
  
  .msg(paste0(
    "Final reference contains ",
    length(unique(promoter_sites$mcols.gene_name)), " genes."
  ), verbose = verbose)
  
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
  
  promoter_sites$mcols.uniqueid <- seq_along(promoter_sites)
  genebody_sites$mcols.uniqueid <- seq_along(genebody_sites)
  
  if (!fastMode) {
    intergenic_sites$mcols.uniqueid <- seq_along(intergenic_sites)
    intragenic_sites$mcols.uniqueid <- seq_along(intragenic_sites)
  }
  
  ## Feature generation
  feature_list <- list()
  
  add_features <- function(frag, prefix) {
    if (is.null(frag)) return(list())
    out <- list(
      generateFeatures(promoter_sites, frag, paste0(prefix, "_promoter")),
      generateFeatures(genebody_sites, frag, paste0(prefix, "_genebody"))
    )
    if (!fastMode)
      out <- c(out, list(
        generateFeatures(intergenic_sites, frag, paste0(prefix, "_intergenic")),
        generateFeatures(intragenic_sites, frag, paste0(prefix, "_intragenic"))
      ))
    out
  }
  
  feature_list <- c(
    feature_list,
    add_features(frag_orig_k4, "K4"),
    add_features(frag_orig_k36, "K36")
  )
  
  if (!fastMode && !is.null(frag_orig_k27))
    feature_list <- c(feature_list, add_features(frag_orig_k27, "K27"))
  
  if (!length(feature_list))
    stop("No features could be generated from the provided fragment files.")
  
  feature_list <- lapply(feature_list, function(x) setnames(x, 1, "Gene"))
  final_df <- Reduce(function(x, y) merge(x, y, by = "Gene", all = TRUE),
                     feature_list)
  
  return(final_df)
}