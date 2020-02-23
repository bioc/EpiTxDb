
#' @rdname positionSequence
setGeneric("positionSequence",
           function(x, order, decreasing)
               standardGeneric("positionSequence"))

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
           function(x, columns, filter, use.names, ...)
               standardGeneric("modifications"))

#' @rdname modifications
setGeneric("modificationsBy",
           function(x, by, ...) standardGeneric("modificationsBy"))

#' @rdname modifications
setGeneric("modifiedSeqsByTranscript",
           function(x, sequences, ...)
               standardGeneric("modifiedSeqsByTranscript"))

