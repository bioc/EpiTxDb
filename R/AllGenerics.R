
setGeneric("positionSequence",
           function(x, ...)
               standardGeneric("positionSequence"))

setGeneric("shiftTranscriptToGenomic",
           function(subject, tx)
               standardGeneric("shiftTranscriptToGenomic"))

setGeneric("shiftGenomicToTranscript",
           function(subject, tx)
               standardGeneric("shiftGenomicToTranscript"))


setGeneric("modifications",
           function(x, ...) standardGeneric("modifications"))

setGeneric("modificationsBy",
           function(x, ...) standardGeneric("modificationsBy"))

setGeneric("modifiedSeqsByTranscript",
           function(x, sequences, ...)
               standardGeneric("modifiedSeqsByTranscript"))

