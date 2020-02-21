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
.get_gr_modifications <- function(gr){
    modifications <- list(mod_id = mcols(gr)$mod_id,
                          mod_type = mcols(gr)$mod_type,
                          mod_name = mcols(gr)$mod_name,
                          mod_start = start(gr),
                          mod_end = end(gr),
                          tx_id = mcols(gr)$tx_id,
                          tx_name= mcols(gr)$tx_name,
                          tx_ensembl = mcols(gr)$tx_ensembl)
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
.get_gr_reactions <- function(gr){
    reactions <- list(mod_id = mcols(gr)$mod_id,
                      rx_genename = mcols(gr)$rx_genename,
                      rx_rank = mcols(gr)$rx_rank,
                      rx_ensembl = mcols(gr)$rx_ensembl,
                      rx_ensembltrans = mcols(gr)$rx_ensembltrans,
                      rx_entrezid = mcols(gr)$rx_entrezid)
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
.get_gr_specifiers <- function(gr){
    specifiers <- list(mod_id = mcols(gr)$mod_id,
                       spec_type = mcols(gr)$spec_type,
                       spec_genename = mcols(gr)$spec_genename,
                       spec_ensembl = mcols(gr)$spec_ensembl,
                       spec_ensembltrans = mcols(gr)$spec_ensembltrans,
                       spec_entrezid = mcols(gr)$spec_entrezid)
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
.get_gr_references <- function(gr){
    ref_lengths <- lengths(mcols(gr)$ref_type)
    references <- list(mod_id = unlist(Map(rep,mcols(gr)$mod_id,ref_lengths)),
                       ref_type = unlist(mcols(gr)$ref_type),
                       ref = unlist(mcols(gr)$ref))
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
    modifications <- .get_gr_modifications(gr)
    reactions <- .get_gr_reactions(gr)
    specifiers <- .get_gr_specifiers(gr)
    references <- .get_gr_references(gr)
    makeEpiTxDb(modifications = modifications, reactions = reactions,
                specifiers = specifiers, references = references,
                metadata = metadata, reassign.ids = reassign.ids)
}