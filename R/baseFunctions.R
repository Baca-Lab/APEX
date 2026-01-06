#Basic functions for analyzing cfChIP-seq fragment files
readFragBed <- function(file){
  df <- suppressMessages(data.frame(data.table::fread(file, col.names = c("chr","start","end","strand","frag_length", "gc_bias", "End_motif_PosStrand", "End_motif_NegStrand"))))
  df <- df[which(df$frag_length >= 20 & df$frag_length <= 500),]
  return(GenomicRanges::makeGRangesFromDataFrame(df=df,keep.extra.columns = TRUE,starts.in.df.are.0based = TRUE))
}

mergeDF <- function(df, targets) {
  if ("UniqueID" %in% colnames(df)) {
    all_rows <- data.frame(UniqueID = targets$mcols.uniqueid,
                           Bin = targets$mcols.bin,
                           Gene = targets$mcols.gene_name)
    df <- left_join(all_rows, df, by = "UniqueID")
    df <- df %>%
      select(-matches("Bin.y|Gene.y")) %>%
      rename(Bin = Bin.x, Gene = Gene.x)
  } else {
    all_rows <- data.frame(Gene = targets$mcols.gene_name, Bin = targets$mcols.bin)
    df <- left_join(all_rows, df, by = c("Gene", "Bin"))
  }
  return(df)
}


tile_site <- function(gr, num_of_tiles = 50L) {
  tiles <- unlist(tile(gr, n = num_of_tiles))
  mcols(tiles) <- rep(mcols(gr), each = num_of_tiles)
  tiles$mcols.bin <- rep(1:num_of_tiles, length(gr))
  tiles$mcols.binsize <- rep(round(width(gr)/num_of_tiles), each = num_of_tiles)
  tiles <- rev_strand(tiles)
  return(tiles)
}

rev_strand <- function(gr){
  rev_gr <- data.frame(gr) %>%
    group_by(mcols.gene_name) %>%
    mutate(mcols.bin = if_else(strand == "-", rev(mcols.bin), mcols.bin)) %>%
    ungroup()
  return(GRanges(rev_gr))
}

getTSS <- function(gr) {
  gr2 <- GRanges(
    seqnames = seqnames(gr),        # Chromosome
    ranges = IRanges(
      start = ifelse(strand(gr) == "+", start(gr), end(gr)),  # Start for +, End for -
      width = 1),
    strand = strand(gr)             # Maintain the strand information
  )
  mcols(gr2) <- data.frame(mcols(gr))
  return(gr2)
}

getTTS <- function(gr) {
  gr <- GRanges(
    seqnames = seqnames(gr),        # Chromosome
    ranges = IRanges(
      start = ifelse(strand(gr) == "+", end(gr), start(gr)),  # Start for +, End for -
      width = 1,
      mcols =data.frame(mcols(gr))# TTS is a single base pair position
    ),
    strand = strand(gr)             # Maintain the strand information
  )
  return(gr)
}

shannon_entropy_endmotif <- function(endmotif) {
  bin_motif <- table(endmotif)
  # Compute probabilities
  total_endmotif <- sum(bin_motif)
  if (total_endmotif == 0) return(NA)  # Avoid division by zero
  p1 <- bin_motif/total_endmotif

  # Remove zero probabilities to avoid log(0)
  p1 <- p1[p1 > 0]

  # Compute Shannon entropy
  shannon_entropy <- -sum(p1 * log2(p1))
  return(shannon_entropy)
}

shannon_entropy_fragLen <- function(frag_lengths, bin_range = 20:500) {
  # Define bin edges
  bin_edges <- seq(min(bin_range), max(bin_range), length.out = max(bin_range) - min(bin_range) + 1)

  # Count fragments in each bin
  bin_counts <- hist(frag_lengths, breaks = bin_edges, plot = FALSE)$counts

  # Compute probabilities
  total_fragments <- sum(bin_counts)
  if (total_fragments == 0) return(NA)  # Avoid division by zero
  p1 <- bin_counts / total_fragments

  # Remove zero probabilities to avoid log(0)
  p1 <- p1[p1 > 0]

  # Compute Shannon entropy
  shannon_entropy <- -sum(p1 * log2(p1))
  return(shannon_entropy)
}

empirical_percentile <- function(x, treat_zeros_as_missing = TRUE) {
  valid_x <- x[!is.na(x)]
  if (treat_zeros_as_missing) valid_x <- valid_x[valid_x != log(0.01, 2)]

  if (length(valid_x) == 0) {
    return(rep(NA_real_, length(x)))  # nothing valid to rank
  }

  ecdf_fn <- ecdf(valid_x)
  percentiles <- ifelse(is.na(x) | (treat_zeros_as_missing & x == 0),
                        NA_real_,
                        ecdf_fn(x))

  return(percentiles)
}


weighted_means_df <- function(df, tgts) {
  tmp <- setDT(cbind.data.frame(df, genes = tgts$mcols.gene_name, bins = tgts$mcols.bin, weights = tgts$score, binsize = tgts$mcols.binsize))
  tmp <- tmp[, lapply(.SD, function(col) weighted.mean(col, weights, na.rm = TRUE)), by = .(genes, bins), .SDcols = c(colnames(df), "binsize")]
  return(tmp)
}


weighted.var <- function(x, w, na.rm = TRUE) {
  if (na.rm) {
    keep <- !is.na(x) & !is.na(w)
    x <- x[keep]
    w <- w[keep]
  }

  # If fewer than 2 values, variance is undefined
  if (length(x) < 2) {
    return(NA_real_)
  }
  w <- w/sum(w)
  w_mean <- weighted.mean(x, w)
  denom <- (sum(w)^2 - sum(w^2)) / sum(w)
  var <- sum(w * (x - w_mean)^2) / denom
  return(var)
}

qualityControl <- function(frags, histone_mark = "H3K4me3") {
  if(unlist(class(frags)) != "GRanges"){
    frags <- readFragBed(frags)
  }

  # Construct file paths
  dir <- "resources/histone_target_sites"
  on_file <- system.file(paste0(dir, "/", histone_mark, "/", "on.target.filt.bed"), package = "apex")
  off_file <- system.file(paste0(dir, "/", histone_mark, "/", "off.target.filt.bed"), package = "apex")

  # Import on/off target regions
  on <- rtracklayer::import(on_file)
  off <- rtracklayer::import(off_file)

  # Compute base pair coverage
  on_bp <- sum(GenomicRanges::end(on) - GenomicRanges::start(on))
  off_bp <- sum(GenomicRanges::end(off) - GenomicRanges::start(off))

  # Compute read overlaps
  on_reads <- sum(GenomicRanges::countOverlaps(frags, on, ignore.strand = TRUE))
  off_reads <- sum(GenomicRanges::countOverlaps(frags, off, ignore.strand = TRUE))

  # Calculate enrichment
  enrichment <- (on_reads / on_bp) / (off_reads / off_bp)
  frag_num <- length(frags)

  # Return as named list
  return(list(
    histone_mark = histone_mark,
    enrichment = enrichment,
    frag_num = frag_num
  ))
}

qc <- function(quality_results) {
  mark <- lapply(quality_results, function(x) x[[1]])
  enrichment_scores <- lapply(quality_results, function(x) x[[2]])
  frag_num <- lapply(quality_results, function(x) x[[3]])

  cat("\n")
  cat("===========================================\n")
  cat("             QC SUMMARY REPORT            \n")
  cat("===========================================\n")
  cat(" Histone Mark   Enrichment   Fragments\n")
  cat("-------------------------------------------------\n")

  qc <- logical(length(enrichment_scores))
  status <- character(length(enrichment_scores))

  for (i in seq_along(enrichment_scores)) {
    enrich <- enrichment_scores[[i]]
    fragments <- frag_num[[i]]
    cat(sprintf(" %-12s   %-11.2f  %-12s\n", mark[[i]], enrich, format(fragments, big.mark = ",")))
}

  cat("=================================================\n\n")
  return(qc)
}

.msg <- function(..., verbose = TRUE) {
  if (isTRUE(verbose)) message(...)
}
