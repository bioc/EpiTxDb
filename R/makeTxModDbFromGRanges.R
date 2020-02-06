#' @include TxModDb-class.R
#' @include makeTxModDb.R
NULL

#' @name makeTxModDbfromtRNAdb
#'
#' @title makeTxModDbfromtRNAdb
#'
#' @description
#' title
#'
#'
NULL

# makeTxModDbfromtRNAdb --------------------------------------------------------

#' @rdname makeTxModDbfromtRNAdb
#' @export
makeTxModDbfromGRanges <- function(gr, metadata = NULL, reassign.ids = FALSE){
    modifications <- data.frame(mod_id = mcols(gr)$mod_id,
                                mod_type = mcols(gr)$mod_type,
                                mod_name = mcols(gr)$mod_name,
                                mod_start = start(gr),
                                mod_end = end(gr),
                                transcript_id = mcols(gr)$transcript_id,
                                transcript_name= mcols(gr)$transcript_name,
                                ensembltrans = mcols(gr)$ensembltrans,
                                check.names = FALSE, stringsAsFactors = FALSE)
    reactions <- list(mod_id = mcols(gr)$mod_id,
                      mod_rank = mcols(gr)$mod_rank,
                      reaction_genename = mcols(gr)$reaction_genename,
                      reaction_ensembl = mcols(gr)$reaction_ensembl,
                      reaction_ensembltrans = mcols(gr)$reaction_ensembltrans,
                      reaction_entrezid = mcols(gr)$reaction_entrezid,
                      reaction_enzyme= mcols(gr)$reaction_enzyme)
    reactions <- reactions[!vapply(reactions, is.null, logical(1))]
    if(length(reactions) > 1L){
        reactions <- data.frame(reactions, check.names = FALSE,
                                stringsAsFactors = FALSE)
    } else {
        reactions <- NULL
    }
    specifier <- list(mod_id = mcols(gr)$mod_id,
                      specifier_type = mcols(gr)$specifier_type,
                      specifier_genename = mcols(gr)$specifier_genename,
                      specifier_entrezid = mcols(gr)$reaction_entrezid,
                      specifier_ensembl = mcols(gr)$specifier_ensembl)
    specifier <- specifier[!vapply(specifier, is.null, logical(1))]
    if(length(specifier) > 1L){
        specifier <- data.frame(specifier, check.names = FALSE,
                                stringsAsFactors = FALSE)
    } else {
        specifier <- NULL
    }
    makeTxModDb(modifications, reactions = reactions, specifier = specifier,
                metadata = metadata, reassign.ids = reassign.ids)
}