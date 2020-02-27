
#' @rdname positionSequence
setGeneric("positionSequence",
           function(x, order = FALSE, decreasing = FALSE)
               standardGeneric("positionSequence"))

#' @rdname rescale
setGeneric("rescale",
           function(x, to = 1L, from = 1L)
               standardGeneric("rescale"))

#' @rdname shiftGenomicToTranscript
setGeneric("shiftTranscriptToGenomic",
           function(subject, tx)
               standardGeneric("shiftTranscriptToGenomic"))

#' @rdname shiftGenomicToTranscript
setGeneric("shiftGenomicToTranscript",
           function(subject, tx)
               standardGeneric("shiftGenomicToTranscript"))


#' @rdname modifications
setGeneric("modifications",
           function(x, columns  = c("mod_id","mod_type","mod_name"),
                    filter = NULL, use.names = FALSE, ...)
               standardGeneric("modifications"))

#' @rdname modifications
setGeneric("modificationsBy",
           function(x, by = c("seqnames","mod_type","reaction","specifier",
                              "specifier_type"),
                    ...)
               standardGeneric("modificationsBy"))

#' @rdname modifications
setGeneric("modifiedSeqsByTranscript",
           function(x, sequences, ...)
               standardGeneric("modifiedSeqsByTranscript"))

