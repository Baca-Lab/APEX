.onLoad <- function(libname, pkgname) {
  utils::globalVariables(c(

    # data.table symbols
    ".", ":=", ".N", ".SD",
    "Weight", "N", "entropy", "freq", "weights",
    "genes", "bins",


    # ggplot2 / dplyr aesthetics
    "gene", "P.Value", "logFC", "negLogP",
    "highlight", "percentile",

    # fragment / feature metadata
    "Bin", "Bin.x", "BinSize", "Gene", "Gene.x", "UniqueID",
    "TargetCoordinates", "RPK",
    "TotalCoverage", "MaxCoverage",
    "frag_length", "gc_bias",
    "end_motif",
    "End_motif_PosStrand", "End_motif_NegStrand",

    # GenomicRanges metadata
    "mcols.gene_name", "mcols.bin", "strand"
  ))
}
