#' @title \code{EpiTxDb} - Storing and accessing epitranscriptomic information
#' using the AnnotationDbi interface
#'
#' @author Felix G M Ernst [aut]
#'
#' @description
#' title
#'
#' @references
#' Jia-Jia Xuan, Wen-Ju Sun, Ke-Ren Zhou, Shun Liu, Peng-Hui Lin, Ling-Ling
#' Zheng, Liang-Hu Qu, Jian-Hua Yang (2017): "RMBase v2.0: Deciphering the Map of
#' RNA Modifications from Epitranscriptome Sequencing Data." Nucleic Acids
#' Research, Volume 46, Issue D1, 4 January 2018, Pages D327–D334.
#' doi: 10.1093/nar/gkx934
#'
#' Jühling, Frank; Mörl, Mario; Hartmann, Roland K.; Sprinzl, Mathias; Stadler,
#' Peter F.; Pütz, Joern (2009): "TRNAdb 2009: Compilation of tRNA Sequences and
#' tRNA Genes." Nucleic Acids Research 37 (suppl_1): D159–D162.
#' doi: 10.1093/nar/gkn772
#'
#' Sprinzl, Mathias; Vassilenko, Konstantin S. (2005): "Compilation of tRNA
#' Sequences and Sequences of tRNA Genes." Nucleic Acids Research 33 (suppl_1):
#' D139–D140. doi: 10.1093/nar/gki012
#'
#' @name EpiTxDb-package#'
NULL

#' @keywords internal
"_PACKAGE"


#' @import methods
#' @import GenomicFeatures
#' @import txdbmaker
#' @import BiocGenerics
#' @import S4Vectors
#' @import IRanges
#' @import AnnotationDbi
#' @import RSQLite
#' @import GenomeInfoDb
#' @import Modstrings
NULL

#' @name EpiTxDb-data
#'
#' @title EpiTxDb internal data
#'
#' @usage data(rmbase_data)
#' @format data.frame
"rmbase_data"