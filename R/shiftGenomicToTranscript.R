#' @include EpiTxDb-class.R
NULL

#' @name shiftGenomicToTranscript
#'
#' @title Shift \code{GRanges} coordinates based on another \code{GRanges}
#' object
#'
#' @description
#' \code{shiftGenomicToTranscript} shifts positions of a
#' \code{\link[GenomicRanges:GRanges-class]{GRanges}}
#' object based on coordinates of another \code{GRanges} object. The most common
#' application is to shift genomic coordinates to transcript coordinates, which
#' is reflected in the name. \code{shiftTranscriptToGenomic} implements the
#' reverse operation.
#'
#' Matches are determined by
#' \code{\link[GenomicRanges:findOverlaps-methods]{findOverlaps}} for
#' \code{shiftGenomicToTranscript} and by
#' \code{\link[S4Vectors:Vector-comparison]{findMatches}} for
#' \code{shiftTranscriptToGenomic} using the \code{seqnames} of the
#' \code{subject} and the \code{names} of \code{tx}.
#'
#' @param subject a \code{\link[GenomicRanges:GRanges-class]{GRanges}} or
#'   \code{\link[GenomicRanges:GRangesList-class]{GRangesList}} object
#' @param tx a named \code{\link[GenomicRanges:GRangesList-class]{GRangesList}}
#'   object.
#'
#' @examples
#' library(GenomicRanges)
#' # Construct some example data
#' subject1 <- GRanges("chr1", IRanges(3, 6),
#'                     strand = "+")
#' subject2 <- GRanges("chr1", IRanges(c(17,23), width=3),
#'                     strand = c("+","-"))
#' subject3 <- GRanges("chr2", IRanges(c(51, 54), c(53, 59)),
#'                     strand = "-")
#' subject <- GRangesList(a=subject1, b=subject2, c=subject3)
#' tx1 <- GRanges("chr1", IRanges(1, 40),
#'                strand="+")
#' tx2 <- GRanges("chr1", IRanges(10, 30),
#'                strand="+")
#' tx3 <- GRanges("chr2", IRanges(50, 60),
#'                strand="-")
#' tx <- GRangesList(a=tx1, b=tx2, c=tx3)
#'
#' # shift to transcript coordinates. Since the third subject does not have
#' # a match in tx it is dropped with a warning
#' shifted_grl <- shiftGenomicToTranscript(subject,tx)
#'
#' # ... and back
#' shifted_grl2 <- shiftTranscriptToGenomic(shifted_grl,tx)
#'
#' # comparison of ranges work. However the seqlevels differ
#' ranges(shifted_grl2) == ranges(subject[list(1,c(1,1),c(1,2))])
NULL

# helper functions -------------------------------------------------------------

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
    tx
}

# shiftTranscriptToGenomic ----------------------------------------------------

#' @importFrom S4Vectors findMatches
.hitsFUN_TX_TO_GENOME <- function(subject,tx){
    hits <- suppressWarnings(S4Vectors::findMatches(names(tx),
                                                    seqnames(subject)))
    hits
}

#' @importFrom GenomicRanges findOverlaps
.hitsFUN_GENOME_TO_TX <- function(subject,tx){
    hits <- suppressWarnings(GenomicRanges::findOverlaps(tx,subject))
    hits
}

.compareFUN_TX_TO_GENOME <- function(y, z, seqs){
    seqs[y][IRanges::IntegerList(as.list(z))]
}

.compareFUN_GENOME_TO_TX <- function(y, z, seqs){
    which(seqs[y] == z)
}

.resultFUN_TX_TO_GENOME <- function(tx_hits, subject_hits, starts, ends){
    # transfer gene information
    mcols(subject_hits)$seq_start <- start(subject_hits)
    mcols(subject_hits)$seq_end <- end(subject_hits)
    mcols(subject_hits)$seq_strand <- strand(subject_hits)
    mcols(subject_hits)$seq_name <- seqnames(subject_hits)
    # prepare reformat of GRanges results
    seqnames <- unlist(unique(seqnames(tx_hits)),use.names = FALSE)
    strand <- unlist(unique(strand(tx_hits)),use.names = FALSE)
    mcols <- mcols(subject_hits)
    tmp <- starts
    starts[strand == "-"] <- ends[strand == "-"]
    ends[strand == "-"] <- tmp[strand == "-"]
    ranges <- IRanges::IRanges(start = unlist(starts), end = unlist(ends))
    # reformat GRanges results
    ans <- GenomicRanges::GRanges(seqnames = seqnames,
                                  ranges = ranges,
                                  strand = strand,
                                  mcols)
    ans
}

.resultFUN_GENOME_TO_TX <- function(tx_hits, subject_hits, starts, ends){
    # transfer transcript information
    mcols(subject_hits)$seq_start <- start(subject_hits)
    mcols(subject_hits)$seq_end <- end(subject_hits)
    mcols(subject_hits)$seq_strand <- strand(subject_hits)
    mcols(subject_hits)$seq_name <- seqnames(subject_hits)
    # reformat GRanges results
    seqnames <- names(tx_hits)
    mcols <- mcols(subject_hits)
    tmp <- starts
    starts[strand(subject_hits) == "-"] <- ends[strand(subject_hits) == "-"]
    ends[strand(subject_hits) == "-"] <- tmp[strand(subject_hits) == "-"]
    ranges <- IRanges::IRanges(start = unlist(starts), end = unlist(ends))
    #
    ans <- GenomicRanges::GRanges(seqnames = seqnames,
                                  ranges = ranges,
                                  strand = "*",
                                  mcols)
    ans
}

.shift <- function(subject, tx, seqs, .hitsFUN, .compareFUN, .resultFUN){
    # match seqnames of subject and names of tx
    hits <- .hitsFUN(subject, tx)
    if(length(hits) == 0L){
        stop("No 'subject'/'tx' matches found.",
             call. = FALSE)
    }
    # assemble results vectors
    qh <- queryHits(hits)
    sh <- subjectHits(hits)
    f <- seq_along(subject) %in% sh
    if(!all(f)){
        warning("Coordinates for ",sum(!f)," ranges of 'subject' not found: '",
                paste(unique(as.character(subject[!f])),collapse = "','"),
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
    starts <- mapply(.compareFUN,
                     split(qh, f),
                     split(start(subject_hits), f),
                     MoreArgs = list(seqs = seqs),
                     SIMPLIFY = FALSE)
    starts <- do.call(c,unname(starts))
    ends <- mapply(.compareFUN,
                   split(qh, f),
                   split(end(subject_hits), f),
                   MoreArgs = list(seqs = seqs),
                   SIMPLIFY = FALSE)
    ends <- do.call(c,unname(ends))
    # drop problematic modifications a bit later
    invalid_pos <- lengths(starts) < 1L | lengths(ends) < 1L
    if(any(invalid_pos)){
        warning("Dropping elements of 'subject' which could not be ",
                "repositioned on the genome.",
                call. = FALSE)
    }
    ans <- .resultFUN(tx_hits[!invalid_pos], subject_hits[!invalid_pos],
                      starts[!invalid_pos], ends[!invalid_pos])
    ans
}


#' @rdname shiftGenomicToTranscript
#' @export
setMethod("shiftTranscriptToGenomic", c("GRanges","GRangesList"),
    function(subject, tx){
        tx <- .norm_tx(tx)
        seqs <- positionSequence(tx)
        .shift(subject, tx, seqs, .hitsFUN_TX_TO_GENOME,
               .compareFUN_TX_TO_GENOME, .resultFUN_TX_TO_GENOME)
    })

#' @rdname shiftGenomicToTranscript
#' @export
setMethod("shiftTranscriptToGenomic", c("GRangesList","GRangesList"),
    function(subject, tx){
        tx <- .norm_tx(tx)
        seqs <- positionSequence(tx)
        ans <- lapply(subject, .shift, tx, seqs, .hitsFUN_TX_TO_GENOME,
                      .compareFUN_TX_TO_GENOME, .resultFUN_TX_TO_GENOME)
        GenomicRanges::GRangesList(ans)
    })


# shiftGenomicToTranscript ----------------------------------------------------

#' @rdname shiftGenomicToTranscript
#' @export
setMethod("shiftGenomicToTranscript", c("GRanges","GRangesList"),
    function(subject, tx){
        tx <- .norm_tx(tx)
        seqs <- positionSequence(tx)
        .shift(subject, tx, seqs, .hitsFUN_GENOME_TO_TX, .compareFUN_GENOME_TO_TX,
               .resultFUN_GENOME_TO_TX)
    })

#' @rdname shiftGenomicToTranscript
#' @export
setMethod("shiftGenomicToTranscript", c("GRangesList","GRangesList"),
    function(subject, tx){
        tx <- .norm_tx(tx)
        seqs <- positionSequence(tx)
        ans <- lapply(subject, .shift, tx, seqs, .hitsFUN_GENOME_TO_TX,
                      .compareFUN_GENOME_TO_TX, .resultFUN_GENOME_TO_TX)
        GenomicRanges::GRangesList(ans)
    })
