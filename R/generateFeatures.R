#' Generate epigenomic features for a set of target regions from a fragment GRanges
#'
#' Given a set of tiled target regions (e.g., promoter or gene-body tiles) and a
#' fragment-level \code{GRanges} (cfChIP-seq fragments with per-fragment metadata),
#' this function aggregates multiple feature classes per gene and per bin,
#' including fragment length statistics, coverage, fragment-length composition,
#' entropy metrics, end-motif frequencies/entropy, and GC bias summaries.
#'
#' @param targets A \code{GRanges} of target intervals (e.g., tiled promoter/gene-body sites).
#'        Must contain metadata columns \code{mcols.gene_name}, \code{mcols.bin},
#'        \code{mcols.binsize}, and \code{mcols.uniqueid}. If \code{targets$score}
#'        is present, it is used as a weight.
#' @param frag A \code{GRanges} of fragments overlapping the reference genome.
#'        The fragment metadata (via \code{mcols(frag)}) must include:
#'        \code{frag_length}, \code{gc_bias}, \code{End_motif_PosStrand},
#'        and \code{End_motif_NegStrand}.
#' @param histone Character string used as a prefix for output column names
#'        (e.g., \code{"K4_promoter"}).
#' @param verbose Logical; if \code{TRUE} (default), progress messages are printed
#'        via \code{.msg}. If \code{FALSE}, messages are suppressed.
#'
#' @return A data.frame-like object (from \pkg{data.table}) with one row per gene
#'         and feature columns prefixed by \code{histone}.
#'
#' @importFrom dplyr left_join case_when
#' @importFrom data.table as.data.table data.table .SD .N setnames rbindlist setDT
#' @importFrom reshape2 dcast
#' @importFrom GenomicRanges findOverlaps seqnames
#' @importFrom IRanges start end
#' @importFrom S4Vectors mcols subjectHits queryHits
#' @importFrom stats model.matrix ecdf reorder weighted.mean
#' @importFrom graphics hist
#' @export
#'
#Main function to generate features from a fragment file given a target region (promoter or gene body)
generateFeatures <- function(targets, frag, histone, verbose = TRUE) {

.msg(paste0("Generating features for ", histone, "..."), verbose = verbose)

  # Map individual fragments to target sites
  ol <- GenomicRanges::findOverlaps(targets, frag, ignore.strand = TRUE)
  df <- as.data.table(
    cbind.data.frame(
      ol,
      mcols(frag[subjectHits(ol)]),
      Coverage = 1,
      Gene = targets[queryHits(ol)]$mcols.gene_name,
      TargetCoordinates = paste(seqnames(targets[queryHits(ol)]), start(targets[queryHits(ol)]), end(targets[queryHits(ol)]), sep = "_"),
      Bin = targets[queryHits(ol)]$mcols.bin,
      BinSize = targets[queryHits(ol)]$mcols.binsize,
      UniqueID = targets[queryHits(ol)]$mcols.uniqueid,
      Weight = if(!is.null(targets$score)) {targets$score[queryHits(ol)]} else 1
    )
  )

  # 1 -- Mean and variance of fragment length per bin
.msg("Calculating mean and variance of fragment lengths per bin...", verbose = verbose)
  fl_stats <- df[, .(Mean_FL = weighted.mean(frag_length, Weight, na.rm = TRUE), Var_FL  = weighted.var(frag_length, Weight, na.rm = TRUE)), by = .(Bin, Gene)]
  fl_stats <- unique(mergeDF(fl_stats, targets))

  fl_mean <- dcast(as.data.table(fl_stats), Gene ~ Bin, value.var = "Mean_FL")
  setnames(fl_mean, c("Gene", setdiff(names(fl_mean), "Gene")))

  fl_var <- dcast(as.data.table(fl_stats), Gene ~ Bin, value.var = "Var_FL")
  setnames(fl_var, c("Gene", setdiff(names(fl_var), "Gene")))

  # 2 -- Frequency of differently sized fragments per bin
.msg("Calculating frequency of differently sized fragments per bin...", verbose = verbose)
  frag_length <- df$frag_length
  fb <- case_when(
      frag_length >= 20 & frag_length <= 80 ~ "FragFreq_20_80",
      frag_length >= 80 & frag_length <= 120 ~ "FragFreq_80_120",
      frag_length >= 160 & frag_length <= 200 ~ "FragFreq_160_200",
      frag_length >= 280 & frag_length <= 320 ~ "FragFreq_280_320",
      frag_length >= 400 & frag_length <= 440 ~ "FragFreq_400_440",
      TRUE ~ "FragFreq_remaining"
    )

  fb <- df[, .N, by = .(UniqueID, Bin, Weight, BinSize, TargetCoordinates, Gene, fb)]
  fb <- dcast(
    as.data.table(fb),
    UniqueID + Bin + BinSize + Weight + TargetCoordinates + Gene ~ fb,
    value.var = "N",
    fill = 0
  )
  fb <- data.table(fb)

  # Normalize coverage
  total_coverage <- fb[, lapply(.SD, sum), .SDcols = grep("FragFreq_", names(fb)), by = .(Gene, Bin, BinSize, Weight)]
  total_coverage[, RPK := rowSums(.SD) / (BinSize / 1e3), .SDcols = grep("^FragFreq_", names(total_coverage), value = TRUE)]
  total_rpk_sum <- sum(total_coverage$RPK, na.rm = TRUE)
  total_coverage[, TotalCoverage := RPK / total_rpk_sum * 1e6]
  total_coverage <- as.data.table(total_coverage)[, .(TotalCoverage = weighted.mean(TotalCoverage, Weight, na.rm = TRUE)), by = .(Gene, Bin)]
  total_coverage <- unique(mergeDF(total_coverage, targets))
  total_coverage <- dcast(
    as.data.table(total_coverage),
    Gene ~ Bin,
    value.var = "TotalCoverage",
    fun.aggregate = sum,
    na.rm = TRUE
  )

  #Maximum Coverage
  max_coverage <- data.table(total_coverage)[, MaxCoverage := apply(.SD, 1, max, na.rm = TRUE), .SDcols = -1]
  max_coverage <- max_coverage[, .SD, .SDcols = c(1, length(max_coverage))]

  #Fragment Length Frequency
  frag_cols <- grep("^FragFreq_", colnames(fb), value = TRUE)
  binned_coverage <- fb[, c("Gene", "Weight", frag_cols), with = FALSE][, lapply(.SD, function(col) weighted.mean(col, Weight, na.rm = TRUE)),  by = .(Gene), .SDcols = frag_cols]
  binned_coverage[, (frag_cols) := lapply(.SD, function(x) x / rowSums(.SD, na.rm = TRUE)), .SDcols = frag_cols]
  binned_coverage <- left_join(data.table(Gene = unique(targets$mcols.gene_name)), binned_coverage, by = "Gene")

  # Set missing coverage values to 0
  total_coverage[is.na(total_coverage)] <- 0
  binned_coverage[is.na(binned_coverage)] <- 0
  max_coverage[is.na(max_coverage)] <- 0

  # 3 -- Shannon entropy of fragment length per bin
.msg("Calculating Shannon entropy of fragment length sizes per bin...", verbose = verbose)
  bin_range <- 20:500
  fent <- df[, .(entropy = shannon_entropy_fragLen(frag_length, bin_range = bin_range) / log2(length(bin_range))), by = .(Bin, Gene, Weight)]
  fent <- fent[, .(entropy = weighted.mean(entropy, Weight, na.rm = TRUE)), by = .(Bin, Gene)]
  fent <- unique(mergeDF(fent, targets))
  fent <- dcast(as.data.table(fent), Gene ~ Bin, value.var = "entropy")
  setnames(fent, c("Gene", setdiff(names(fent), "Gene")))

  # 4 -- Frequency of different 4-mer motifs per region
.msg("Calculating 4-mer motif frequency across region...", verbose = verbose)
  nucleotides <- c("A", "C", "T", "G")
  combinations <- expand.grid(nucleotides, nucleotides, nucleotides, nucleotides)
  all_combinations <- sort(apply(combinations, 1, paste0, collapse = ""))

  # Combine both strands' motifs
  mf <- rbindlist(list(df[, .(Gene, Bin, Weight, end_motif = End_motif_PosStrand)], df[, .(Gene, Bin, Weight, end_motif = End_motif_NegStrand)]))
  mf[, end_motif := factor(end_motif, levels = all_combinations)]
  mf <- mf[!is.na(end_motif)]

  # Calculate frequency table
  mf_freq <- mf[, .N, by = .(Gene, Weight, end_motif)]
  mf_freq[, freq := N / sum(N), by = .(Gene, Weight)]
  mf_freq <- unique(mf_freq[, .(freq = weighted.mean(freq, Weight, na.rm = TRUE)), by = .(Gene, end_motif)])
  mf_freq <- dcast(as.data.table(mf_freq), Gene ~ end_motif, value.var = "freq", fill = 0)
  mf_freq <- left_join(data.table(Gene = unique(targets$mcols.gene_name)), mf_freq, by = "Gene")

  # 5 -- Shannon entropy of 4-mer motifs per bin
.msg("Calculating 4-mer motif entropy across bins...", verbose = verbose)

  mf_ent <- mf[,.(entropy = shannon_entropy_endmotif(end_motif) / log2(length(levels(end_motif)))),by = .(Bin, Gene, Weight)]
  mf_ent <- mf_ent[,.(entropy = weighted.mean(entropy, Weight, na.rm = TRUE)), by = .(Bin, Gene)]
  mf_ent <- unique(mergeDF(mf_ent, targets))
  mf_ent <- dcast(as.data.table(mf_ent), Gene ~ Bin, value.var = "entropy")
  setnames(mf_ent, c("Gene", setdiff(names(mf_ent), "Gene")))

  # 6 -- Average GC % of fragments
.msg("Calculating average GC% of fragments per bin...", verbose = verbose)
  gcbias <- df[, .(GCperc = weighted.mean(gc_bias, Weight, na.rm = TRUE)), by = .(Bin, Gene)]
  gcbias <- unique(mergeDF(gcbias, targets))
  gcbias <- dcast(as.data.table(gcbias), Gene ~ Bin, value.var = "GCperc")
  setnames(gcbias, c("Gene", setdiff(names(gcbias), "Gene")))

  # 7 -- Final feature compilation
.msg("Compiling all extracted features into final table...", verbose = verbose)
  combined <- cbind(
    Gene = fl_mean$Gene,
    FragLen_Mean = fl_mean[, -1],
    FragLen_Var = fl_var[, -1],
    Total_Coverage = log(total_coverage[,-1] + 0.01, 2),  # Log-transform with 0.01 pseudocount
    Max_Coverage = unlist(log(max_coverage[, -1] + 0.01, 2)),
    binned_coverage[, -1],
    Frag_Entropy = fent[,-1],
    Endmotif_Freq = mf_freq[,-1],
    Endmotif_Entropy = mf_ent[, -1],
    GCperc = gcbias[, -1]
  )
  colnames(combined) <- paste0(histone, "_", colnames(combined))

.msg(paste0("Feature extraction complete for ", histone, "!\n"), verbose = verbose)
  return(combined)
}
