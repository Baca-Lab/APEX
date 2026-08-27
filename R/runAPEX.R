#' APEX (Associating Plasma Epigenomics with eXpression), framework to predict gene expression from plasma epigenomics
#'
#' This function extracts features from fragment files and generates (1) a dataframe
#' with several epigenetic features for each gene and (2) a numeric vector of predicted gene expression.
#'
#' @importFrom dplyr left_join group_by mutate select rename case_when ungroup if_else matches
#' @importFrom magrittr %>%
#' @importFrom data.table as.data.table data.table .SD .N := setnames rbindlist
#' @importFrom reshape2 dcast
#' @importFrom GenomicRanges resize promoters GRanges seqnames strand width findOverlaps tile
#' @importFrom IRanges IRanges mid countOverlaps start end
#' @importFrom GenomeInfoDb seqnames keepSeqlevels seqlevels
#' @importFrom S4Vectors mcols mcols<- subjectHits queryHits
#' @importFrom stats predict
#' @importFrom xgboost xgb.DMatrix xgb.load
#' @importFrom utils globalVariables
#'
#' @param frag_file_k4 Optional path to an H3K4me3 fragment file.
#'   H3K4me3 is recommended and, when available, is used for alternate-promoter
#'   selection. At least one supported histone-mark fragment file must be supplied.
#' @param frag_file_k36 Path to the H3K36me3 fragment file (optional).
#' @param frag_file_k27 Path to the H3K27ac fragment file (optional).
#' @param upstream Integer. Number of base pairs upstream to include in promoter regions. Default: 3000.
#' @param downstream Integer. Number of base pairs downstream to include in promoter regions. Default: 3000.
#' @param custom_coordinates Optional path to a BED file containing
#'        user-defined genomic coordinates for each gene. The file must contain
#'        a column named \code{gene_name} and strand information. These
#'        coordinates are used to define gene-body regions, and transcription
#'        start sites are derived from the coordinates according to strand.
#' @param num_of_promoter_tiles Integer. Number of promoter tiles to generate. Default: 20.
#' @param num_of_genebody_tiles Integer. Number of gene body tiles to generate. Default: 20.
#' @param num_of_intergenic_tiles Integer. Number of intergenic enhancer tiles to generate. Default: 1.
#' @param num_of_intragenic_tiles Integer. Number of intragenic enhancer tiles to generate. Default: 1.
#' @param fastMode Logical. If TRUE (default), runs an efficient model that excludes H3K27ac, intergenic enhancer, and intragenic enhancer sites.
#'   If FALSE, these features are included, which may improve performance but increase runtime.
#' @param useAltPromoter Logical. Whether to use alternate promoters based on coverage. Default: TRUE.
#' @param verbose Logical. If TRUE, prints progress messages. Default: TRUE.
#'
#' @return A named numeric vector of predicted gene expression values.
#' @export
apex <- function(frag_file_k4,
                 frag_file_k36 = NULL,
                 frag_file_k27 = NULL,
                 upstream = 3000,
                 downstream = 3000,
                 custom_coordinates = NULL,
                 num_of_promoter_tiles = 20L,
                 num_of_genebody_tiles = 20L,
                 num_of_intergenic_tiles = 1L,
                 num_of_intragenic_tiles = 1L,
                 fastMode = TRUE,
                 useAltPromoter = TRUE,
                 minNorm = TRUE,
                 verbose = TRUE) {

#Suppress data table outut
options(
    datatable.verbose = FALSE,
    datatable.showProgress = FALSE
  )

  .msg("Running APEX...", verbose = verbose)

  features <- extractFeatures(
    frag_file_k4,
    frag_file_k36 = frag_file_k36,
    frag_file_k27 = frag_file_k27,
    upstream = upstream,
    downstream = downstream,
    num_of_promoter_tiles = num_of_promoter_tiles,
    num_of_genebody_tiles = num_of_genebody_tiles,
    num_of_intergenic_tiles = num_of_intergenic_tiles,
    num_of_intragenic_tiles = num_of_intragenic_tiles,
    fastMode = fastMode,
    custom_coordinates = custom_coordinates,
    useAltPromoter = useAltPromoter,
    verbose = verbose  # <-- key: pass through
  )

  gene_names <- features[, 1]
  features <- apply(features[, -1], 2, empirical_percentile)
  rownames(features) <- gene_names

  .msg("Loading model...", verbose = verbose)

  if (fastMode == FALSE && !is.null(frag_file_k4) && !is.null(frag_file_k27) && !is.null(frag_file_k36)) {
    xgboost_model <- xgboost::xgb.load(system.file("resources/APEX_full_model.ubj", package = "apex"))
  } else if (fastMode == TRUE && !is.null(frag_file_k4) && !is.null(frag_file_k36)) {
    xgboost_model <- xgboost::xgb.load(system.file("resources/APEX_K4_K36_model.ubj", package = "apex"))
  } else if (fastMode == TRUE && !is.null(frag_file_k4) && is.null(frag_file_k36)) {
    xgboost_model <- xgboost::xgb.load(system.file("resources/APEX_K4_model.ubj", package = "apex"))
  } else if (fastMode == TRUE && is.null(frag_file_k4) && !is.null(frag_file_k36)) {
    xgboost_model <- xgboost::xgb.load(system.file("resources/APEX_K36_model.ubj", package = "apex"))
  } else if (fastMode == FALSE && is.null(frag_file_k4) && is.null(frag_file_k36) && !is.null(frag_file_k27)){
    xgboost_model <- xgboost::xgb.load(system.file("resources/APEX_K27_custom_model.ubj", package = "apex"))
  } else if (fastMode == FALSE && !is.null(frag_file_k4) && is.null(frag_file_k36) && !is.null(frag_file_k27)){
    xgboost_model <- xgboost::xgb.load(system.file("resources/APEX_K4_K27_custom_model.ubj", package = "apex"))
  } else {
    stop("Invalid combination of inputs: require frag_file_k4 and/or frag_file_k36, and fastMode settings consistent with provided files.")
  }

  .msg("Predicting expression...", verbose = verbose)

  dtable <- xgboost::xgb.DMatrix(data = as.matrix(features))
  apex_prediction <- stats::predict(xgboost_model, dtable)
  if(minNorm == TRUE){
  apex_prediction <- apex_prediction - min(apex_prediction)
  }
  names(apex_prediction) <- rownames(features)

  .msg("Done.", verbose = verbose)
  return(data.matrix(apex_prediction))
}

