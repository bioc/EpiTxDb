#' @include EpiTxDb-class.R
#' @include makeEpiTxDb.R
NULL

#' @name makeEpiTxDbFromGRanges
#'
#' @title Create a \code{EpiTxDb} object from a \code{GRanges} object
#'
#' @description
#' \code{makeEpiTxDbFromGRanges} extracts informations from a
#' \code{\link[GenomicRanges:GRanges-class]{GRanges}} object. The following
#' metadata columns can be used:
#' \itemize{
#' \item{\code{mod_id}, \code{mod_type}, \code{mod_name} and \code{tx_ensembl}.
#' The first three are mandatory, whereas \code{tx_ensembl} is optional.}
#' \item{\code{rx_genename}, \code{rx_rank}, \code{rx_ensembl},
#' \code{rx_ensembltrans} and \code{rx_entrezid}}
#' \item{\code{spec_type}, \code{spec_genename}, \code{spec_ensembl},
#' \code{spec_ensembltrans} and \code{spec_entrezid}}
#' \item{\code{ref_type} and \code{ref}}
#' }
#' ... and passed on the \code{\link[=makeEpiTxDb]{makeEpiTxDb}}.
#'
#' @param gr A \code{\link[GenomicRanges:GRanges-class]{GRanges}} object, which
#' contains at least the mandatory columns.
#' @param metadata A 2-column \code{data.frame} containing meta information to
#' be included in the \code{EpiTxDb} object. This \code{data.frame} is just
#' passed to \code{\link[=makeEpiTxDb]{makeEpiTxDb}}. See
#' \code{\link[=makeEpiTxDb]{makeEpiTxDb}} for more information about the format
#' of metadata. (default: \code{metadata = NULL})
#' @param reassign.ids = FALSE
#'
#' @return a \code{EpiTxDb} object.
#'
#' @export
#'
#' @examples
#' library(GenomicRanges)
#' gr <- GRanges(seqnames = "test",
#'               ranges = IRanges::IRanges(1,1),
#'               strand = "+",
#'               DataFrame(mod_id = 1L,
#'                         mod_type = "Am",
#'                         mod_name = "Am_1"))
#' etdb <- makeEpiTxDbFromGRanges(gr)
NULL

# makeEpiTxDbfromGRanges -------------------------------------------------------

# extract modification table
.get_gr_modifications <- function(gr){
    modifications <- list(mod_id = mcols(gr)$mod_id,
                          mod_type = mcols(gr)$mod_type,
                          mod_name = mcols(gr)$mod_name,
                          mod_start = start(gr),
                          mod_end = end(gr),
                          sn_id = as.integer(seqnames(gr)),
                          sn_name= seqnames(gr))
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
    if(length(ref_lengths) == 0L){
        return(NULL)
    }
    references <- list(mod_id = unlist(Map(rep,mcols(gr)$mod_id,ref_lengths)),
                       ref_type = unlist(mcols(gr)$ref_type),
                       ref = unlist(mcols(gr)$ref))
    references <- references[!vapply(references, is.null, logical(1))]
    if(length(references) == 0L){
        return(NULL)
    }
    data.frame(references, check.names = FALSE,
               stringsAsFactors = FALSE)
}

#' @rdname makeEpiTxDbFromGRanges
#' @export
makeEpiTxDbFromGRanges <- function(gr, metadata = NULL, reassign.ids = FALSE){
    modifications <- .get_gr_modifications(gr)
    reactions <- .get_gr_reactions(gr)
    specifiers <- .get_gr_specifiers(gr)
    references <- .get_gr_references(gr)
    makeEpiTxDb(modifications = modifications, reactions = reactions,
                specifiers = specifiers, references = references,
                metadata = metadata, reassign.ids = reassign.ids)
}