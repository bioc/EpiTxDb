
setGeneric("positionSequence",
           function(x, ...)
               standardGeneric("positionSequence"))

setGeneric("shiftToGenomicCoordinates",
           function(subject, tx)
               standardGeneric("shiftToGenomicCoordinates"))

setGeneric("shiftToTranscriptCoordinates",
           function(subject, tx)
               standardGeneric("shiftToTranscriptCoordinates"))


setGeneric("modifications",
           function(x, ...) standardGeneric("modifications"))

setGeneric("modificationsBy",
           function(x, ...) standardGeneric("modificationsBy"))

setGeneric("modifiedSeqsByTranscript",
           function(x, sequences, ...)
               standardGeneric("modifiedSeqsByTranscript"))

