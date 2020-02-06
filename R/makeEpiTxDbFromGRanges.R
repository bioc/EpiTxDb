#' @include EpiTxDb-class.R
#' @include makeEpiTxDb.R
NULL

#' @name makeEpiTxDbfromGRanges
#'
#' @title makeEpiTxDbfromGRanges
#'
#' @description
#' title
#'
#'
NULL

# makeEpiTxDbfromGRanges --------------------------------------------------------

#' @rdname makeEpiTxDbfromGRanges
#' @export
makeEpiTxDbfromGRanges <- function(gr, metadata = NULL, reassign.ids = FALSE){
    # extract modification table
    modifications <- data.frame(mod_id = mcols(gr)$mod_id,
                                mod_type = mcols(gr)$mod_type,
                                mod_name = mcols(gr)$mod_name,
                                mod_start = start(gr),
                                mod_end = end(gr),
                                transcript_id = mcols(gr)$transcript_id,
                                transcript_name= mcols(gr)$transcript_name,
                                ensembltrans = mcols(gr)$ensembltrans,
                                check.names = FALSE, stringsAsFactors = FALSE)
    # extract reactions table
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
    # extract specifiers table
    specifiers <- list(mod_id = mcols(gr)$mod_id,
                      specifier_type = mcols(gr)$specifier_type,
                      specifier_genename = mcols(gr)$specifier_genename,
                      specifier_entrezid = mcols(gr)$reaction_entrezid,
                      specifier_ensembl = mcols(gr)$specifier_ensembl)
    specifiers <- specifiers[!vapply(specifiers, is.null, logical(1))]
    if(length(specifiers) > 1L){
        specifiers <- data.frame(specifiers, check.names = FALSE,
                                stringsAsFactors = FALSE)
    } else {
        specifiers <- NULL
    }
    # extract references table
    ref_lengths <- lengths(mcols(gr)$reference_type)
    references <- list(mod_id = unlist(Map(rep,mcols(gr)$mod_id,ref_lengths)),
                       reference_type = unlist(mcols(gr)$reference_type),
                       reference = unlist(mcols(gr)$reference))
    references <- references[!vapply(references, is.null, logical(1))]
    if(length(references) > 1L){
        references <- data.frame(references, check.names = FALSE,
                                 stringsAsFactors = FALSE)
    } else {
        references <- NULL
    }
    makeEpiTxDb(modifications, reactions = reactions, specifiers = specifiers,
                references = references, metadata = metadata,
                reassign.ids = reassign.ids)
}