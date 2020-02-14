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

# makeEpiTxDbfromGRanges -------------------------------------------------------

# extract modification table
.extract_modifications <- function(gr){
    modifications <- list(mod_id = mcols(gr)$mod_id,
                          mod_type = mcols(gr)$mod_type,
                          mod_name = mcols(gr)$mod_name,
                          mod_start = start(gr),
                          mod_end = end(gr),
                          transcript_id = mcols(gr)$transcript_id,
                          transcript_name= mcols(gr)$transcript_name,
                          ensembltrans = mcols(gr)$ensembltrans)
    modifications <- modifications[!vapply(modifications, is.null, logical(1))]
    if(length(modifications) > 1L){
        modifications <- data.frame(modifications, check.names = FALSE,
                                    stringsAsFactors = FALSE)
    } else {
        stop("Couldn't extranct modification information from 'GRanges' object.",
             call. = FALSE)
    }
    modifications
}

# extract reactions table
.extract_reactions <- function(gr){
    reactions <- list(mod_id = mcols(gr)$mod_id,
                      mod_rank = mcols(gr)$mod_rank,
                      reaction_genename = mcols(gr)$reaction_genename,
                      reaction_ensembl = mcols(gr)$reaction_ensembl,
                      reaction_ensembltrans = mcols(gr)$reaction_ensembltrans,
                      reaction_entrezid = mcols(gr)$reaction_entrezid)
    reactions <- reactions[!vapply(reactions, is.null, logical(1))]
    if(length(reactions) > 1L){
        reactions <- data.frame(reactions, check.names = FALSE,
                                stringsAsFactors = FALSE)
    } else {
        reactions <- NULL
    }
    reactions
}

# extract specifiers table
.extract_specifiers <- function(gr){
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
    specifiers
}

# extract references table
.extract_references <- function(gr){
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
    references
}

#' @rdname makeEpiTxDbfromGRanges
#' @export
makeEpiTxDbfromGRanges <- function(gr, metadata = NULL, reassign.ids = FALSE){
    modifications <- .extract_modifications(gr)
    reactions <- .extract_reactions(gr)
    specifiers <- .extract_specifiers(gr)
    references <- .extract_references(gr)
    makeEpiTxDb(modifications = modifications, reactions = reactions,
                specifiers = specifiers, references = references,
                metadata = metadata, reassign.ids = reassign.ids)
}