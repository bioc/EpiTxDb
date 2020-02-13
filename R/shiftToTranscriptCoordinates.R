#' @include EpiTxDb-class.R
NULL

#' @name shiftToTranscriptCoordinates
#'
#' @title Shift \code{GRanges} coordinates based on another \code{GRanges}
#' object
#'
#' @description
#' \code{shiftToTranscriptCoordinates} shifts positions of a
#' \code{\link[GenomicRanges:GRanges-class]{GRanges}}
#' object based on coordinates of another \code{GRanges} object. The most common
#' application is to shift genomic coordinates to transcript coordinates, which
#' is reflected in the name. \code{shiftToGenomicCoordinates} implements the
#' reverse operation.
#'
#' @param subject a \code{\link[GenomicRanges:GRanges-class]{GRanges}} or
#' \code{\link[GenomicRanges:GRangesList-class]{GRangesList}} object
#' @param tx a \code{\link[GenomicRanges:GRangesList-class]{GRangesList}} object
#' with a \code{type} metadata column. In addition the metadata \code{mcols(tx)}
#' must contain at least two columns \code{gene_name} and \code{transcript_name}
#' defining the relation of gene and transcript names.
#'
#' @examples
#' # Construct some example data
#' gr1 <- GRanges("chr2", IRanges(3, 6),
#'                strand="+", gene_name=c("c_g"))
#' gr2 <- GRanges(c("chr1", "chr1"), IRanges(c(7,13), width=3),
#'                strand=c("+", "-"), gene_name=c("a_g","b_g"))
#' gr3 <- GRanges(c("chr1", "chr2"), IRanges(c(1, 4), c(3, 9)),
#'                strand=c("-", "-"), gene_name=c("a_g","c_g"))
#' grl <- GRangesList(gr1=gr1, gr2=gr2, gr3=gr3)
#' tx1 <- GRanges("chr1", IRanges(1, 10),
#'                strand="+", type="exon")
#' tx2 <- GRanges("chr1", IRanges(10, 20),
#'                strand="+", type="exon")
#' tx3 <- GRanges("chr2", IRanges(1, 10),
#'                strand="-", type="exon")
#' tx <- GRangesList(a=tx1, b=tx2, c=tx3)
#' mcols(tx) <- DataFrame(gene_name = c("a_g","b_g","c_g"),
#'                        transcript_name = c("a","b","c"))
#'
#' # shift to transcript coordinates
#' shifted_grl <- shiftToTranscriptCoordinates(grl,tx)
#' # ... and back
#' shifted_grl2 <- shiftToGenomicCoordinates(shifted_grl,tx)
#' # comparison of ranges work. However the seqlevels differ
#' ranges(shifted_grl2) == ranges(grl)
NULL

# helper functions -------------------------------------------------------------

.REQUIRED_SUBJECT_COLNAMES_TO_TRANS <- c("gene_name")
.REQUIRED_SUBJECT_COLNAMES_TO_GENOMIC <- c("transcript_name")
.REQUIRED_TX_COLNAMES <- c("type")
.REQUIRED_TX_MCOLS_COLNAMES <- c("gene_name","transcript_name")

.norm_subject <- function(subject, .required_columns){
    if(is(subject,"subjectangesList")){
        tmp_subject <- unlist(subject)
    } else {
        tmp_subject <- subject
    }
    if(!all(.required_columns %in% colnames(mcols(tmp_subject)))){
        stop("'subject' must contain the following columns: '",
             paste(.required_columns, collapse = "','"),
             "'",
             call. = FALSE)
    }
    subject
}

.norm_tx <- function(tx){
    if(missing(tx) || !is(tx,"GRangesList")){
        stop("'tx' must be a named 'GRangesList'.", call. = FALSE)
    }
    if(is.null(names(tx))){
        stop("'tx' must be a named 'GRangesList'.", call. = FALSE)
    }
    if(anyDuplicated(names(tx))){
        stop("names of 'tx' must be unique.", call. = FALSE)
    }
    if(!all(.REQUIRED_TX_COLNAMES %in% colnames(mcols(unlist(tx))))){
        stop("'tx' must contain the following columns: '",
             paste(.REQUIRED_TX_COLNAMES, collapse = "','"),
             "'",
             call. = FALSE)
    }
    unique_types <- unique(mcols(tx, level="within")[,"type"])
    if(any(lengths(unique_types) != 1L)){
        stop("'type' must be a unique value within each element",
             call. = FALSE)
    }
    mc <- mcols(tx)
    if(!all(.REQUIRED_TX_MCOLS_COLNAMES %in% colnames(mc))){
        stop("'tx' must contain the following columns in its metadata columns",
             ": '",
             paste(.REQUIRED_TX_MCOLS_COLNAMES, collapse = "','"),
             "'",
             call. = FALSE)
    }
    if(!all(mcols(tx)$transcript_name == names(tx))){
        if(!all(names(tx) %in% mcols(tx)$transcript_name)){
            stop("All names(tx) must be in mcols(tx)$transcript_name.",
                 call. = FALSE)
        }
        if(anyDuplicated(mcols(tx)$transcript_name)){
            stop("mcols(tx)$transcript_name must be unique and match names(tx).",
                 call. = FALSE)
        }
        mcols(tx) <- mcols(tx)[match(names(tx),mcols(tx)$transcript_name),]
    }
    tx
}

# shiftToGenomicCoordinates ----------------------------------------------------

.shiftToGenomicCoordinates <- function(subject, tx){
    # expand by gene names
    transcript_name <- IRanges::CharacterList(strsplit(as.character(mcols(subject)$transcript_name),","))
    subject <- subject[unlist(Map(rep,seq_along(transcript_name),lengths(transcript_name)))]
    mcols(subject)$transcript_name <- unlist(transcript_name)
    # match seqnames of subject and names of tx
    hits_names <- suppressWarnings(findMatches(mcols(tx)$transcript_name,
                                               as.character(mcols(subject)$transcript_name)))
    n <- length(hits_names)
    f_non_hit <- seq_along(subject)
    f_non_hit <- f_non_hit[!(f_non_hit %in% subjectHits(hits_names))]
    # match the remaining elements be coordinates
    hits_coordinates <- new("SortedByQueryHits")
    if(length(f_non_hit) > 0L){
        hits_coordinates <- suppressWarnings(findOverlaps(tx,subject[f_non_hit]))
        n <- n + length(hits_coordinates)
    }
    # assemble results vectors
    if(n == 0L){
        stop("No 'subject'/'tx' matches found.",
             call. = FALSE)
    }
    qh <- c(queryHits(hits_names),queryHits(hits_coordinates))
    sh <- c(subjectHits(hits_names), f_non_hit[subjectHits(hits_coordinates)])
    f <- seq_along(subject) %in% sh
    if(!all(f)){
        warning("Coordinates for some ranges of 'subject' not found: '",
                paste(unique(mcols(subject)$transcript_name[!f]),collapse = "','"),
                "'.",
                call. = FALSE)
    }
    tx_hits <- tx[qh]
    subject_hits <- subject[sh]
    # reposition modification on transcript
    chunk_size <- 10^5
    chunks_n <- ceiling(length(qh) / chunk_size)
    f <- factor(unlist(Map(rep.int,seq_len(chunks_n),chunk_size)))
    f <- f[seq_along(qh)]
    seqs <- positionSequence(tx)
    starts <- Map(
        function(y, z){
            seqs[y][IRanges::IntegerList(as.list(z))]
        },
        split(qh, f),
        split(start(subject_hits), f))
    starts <- do.call(c,unname(starts))
    # drop problematic modifications a bit later
    invalid_pos <- lengths(starts) < 1L
    if(any(invalid_pos)){
        warning("Dropping elements of 'subject' which could not be ",
                "repositioned on the genome.",
                call. = FALSE)
    }
    tx_hits <- tx_hits[!invalid_pos]
    subject_hits <- subject_hits[!invalid_pos]
    starts <- starts[!invalid_pos]
    # transfer gene information
    mcols(subject_hits)$gene_name <- ""
    mcols(subject_hits)$gene_name <- mcols(tx_hits)$gene_name
    # prepare reformat of GRanges results
    seqnames <- unlist(unique(seqnames(tx_hits)),use.names = FALSE)
    strand <- unlist(unique(strand(tx_hits)),use.names = FALSE)
    mcols <- mcols(subject_hits)[,!(colnames(mcols(subject)) %in% "transcript_name"),
                            drop=FALSE]
    ranges <- IRanges::IRanges(start = unlist(starts),
                               width = width(subject_hits))
    # reformat GRanges results
    ans <- GenomicRanges::GRanges(seqnames = seqnames,
                                  ranges = ranges,
                                  strand = strand,
                                  mcols)
    ans
}

#' @rdname shiftToTranscriptCoordinates
#' @export
setMethod("shiftToGenomicCoordinates", c("GRanges","GRangesList"),
    function(subject, tx){
        subject <- .norm_subject(subject, .REQUIRED_SUBJECT_COLNAMES_TO_GENOMIC)
        tx <- .norm_tx(tx)
        .shiftToGenomicCoordinates(subject, tx)
    })

#' @rdname shiftToTranscriptCoordinates
#' @export
setMethod("shiftToGenomicCoordinates", c("GRangesList","GRangesList"),
    function(subject, tx){
        subject <- .norm_subject(subject, .REQUIRED_SUBJECT_COLNAMES_TO_GENOMIC)
        tx <- .norm_tx(tx)
        ans <- lapply(subject,.shiftToGenomicCoordinates,tx)
        GenomicRanges::GRangesList(ans)
    })


# shiftToTranscriptCoordinates ----------------------------------------------------

.shiftToTranscriptCoordinates <- function(subject, tx){
    # expand by gene names
    gene_name <- IRanges::CharacterList(strsplit(as.character(mcols(subject)$gene_name),","))
    subject <- subject[unlist(Map(rep,seq_along(gene_name),lengths(gene_name)))]
    mcols(subject)$gene_name <- unlist(gene_name)
    # match seqnames of subject and names of tx
    hits_names <- suppressWarnings(findMatches(mcols(tx)$gene_name,
                                               as.character(mcols(subject)$gene_name)))
    n <- length(hits_names)
    f_non_hit <- seq_along(subject)
    f_non_hit <- f_non_hit[!(f_non_hit %in% subjectHits(hits_names))]
    # match the remaining elements be coordinates
    hits_coordinates <- new("SortedByQueryHits")
    if(length(f_non_hit) > 0L){
        hits_coordinates <- suppressWarnings(findOverlaps(tx,subject[f_non_hit]))
        n <- n + length(hits_coordinates)
    }
    # assemble results vectors
    if(n == 0L){
        stop("No 'subject'/'tx' matches found.",
             call. = FALSE)
    }
    qh <- c(queryHits(hits_names),queryHits(hits_coordinates))
    sh <- c(subjectHits(hits_names), f_non_hit[subjectHits(hits_coordinates)])
    f <- seq_along(subject) %in% sh
    if(!all(f)){
        warning("Coordinates for some ranges of 'subject' not found: '",
                paste(unique(mcols(subject)$gene_name[!f]),collapse = "','"),
                "'.",
                call. = FALSE)
    }
    tx_hits <- tx[qh]
    subject_hits <- subject[sh]
    # reposition modification on transcript because of size constraint split it
    # up into batches
    chunk_size <- 10^5
    chunks_n <- ceiling(length(qh) / chunk_size)
    f <- factor(unlist(Map(rep.int,seq_len(chunks_n),chunk_size)))
    f <- f[seq_along(qh)]
    seqs <- positionSequence(tx)
    starts <- Map(
        function(y, z){
            which(seqs[y] == z)
        },
        split(qh, f),
        split(start(subject_hits), f))
    starts <- do.call(c,unname(starts))
    # drop problematic modifications
    invalid_pos <- lengths(starts) < 1L
    if(any(invalid_pos)){
        warning("Dropping elements of 'subject' which could not be ",
                "repositioned on the transcript.",
                call. = FALSE)
    }
    tx_hits <- tx_hits[!invalid_pos]
    subject_hits <- subject_hits[!invalid_pos]
    starts <- starts[!invalid_pos]
    # transfer transcript information
    mcols(subject_hits)$transcript_name <- ""
    mcols(subject_hits)$transcript_name <- mcols(tx_hits)$transcript_name
    # reformat GRanges results
    seqnames <- mcols(subject_hits)$transcript_name
    if(anyNA(seqnames)){
        stop("Some 'transcript_name' matching the ranges of 'subject' are NA.",
             call. = FALSE)
    }
    mcols <- mcols(subject_hits)[,!(colnames(mcols(subject_hits)) %in% "gene_name"),
                            drop=FALSE]
    seqnames <- seqnames
    ranges <- IRanges::IRanges(start = unlist(starts),
                               width = width(subject_hits))
    ans <- GenomicRanges::GRanges(seqnames = seqnames,
                                  ranges = ranges,
                                  strand = "+",
                                  mcols)
    ans
}

#' @rdname shiftToTranscriptCoordinates
#' @export
setMethod("shiftToTranscriptCoordinates", c("GRanges","GRangesList"),
    function(subject, tx){
        subject <- .norm_subject(subject, .REQUIRED_SUBJECT_COLNAMES_TO_TRANS)
        tx <- .norm_tx(tx)
        .shiftToTranscriptCoordinates(subject, tx)
    })

#' @rdname shiftToTranscriptCoordinates
#' @export
setMethod("shiftToTranscriptCoordinates", c("GRangesList","GRangesList"),
    function(subject, tx){
        subject <- .norm_subject(subject, .REQUIRED_SUBJECT_COLNAMES_TO_TRANS)
        tx <- .norm_tx(tx)
        ans <- lapply(subject,.shiftToTranscriptCoordinates,tx)
        GenomicRanges::GRangesList(ans)
    })
