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
#' @param subject a \code{\link[GenomicRanges:GRanges-class]{GRanges}} or
#' \code{\link[GenomicRanges:GRangesList-class]{GRangesList}} object
#' @param tx a \code{\link[GenomicRanges:GRangesList-class]{GRangesList}} object
#' with a \code{type} metadata column. In addition the metadata \code{mcols(tx)}
#' must contain at least two columns \code{gene_name} and \code{transcript_name}
#' defining the relation of gene and transcript names.
#'
#' @examples
#' \dontrun{
#' library(GenomicRanges)
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
#'                        tx_name = c("a","b","c"))
#'
#' # shift to transcript coordinates
#' shifted_grl <- shiftGenomicToTranscript(grl,tx)
#' # ... and back
#' shifted_grl2 <- shiftTranscriptToGenomic(shifted_grl,tx)
#' # comparison of ranges work. However the seqlevels differ
#' ranges(shifted_grl2) == ranges(grl)
#' }
NULL

# helper functions -------------------------------------------------------------

.REQUIRED_TX_COLNAMES <- c("tx_id")

.norm_subject <- function(subject){
    if(is(subject,"GRangesList")){
        tmp_subject <- unlist(subject)
    } else {
        tmp_subject <- subject
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
             "'.",
             call. = FALSE)
    }
    if(anyNA(unlist(mcols(tx, level="within")[,"tx_id"]))){
        stop("'tx_id' must be a character vector with no NAs.",
             call. = FALSE)
    }
    tx <- tx[lengths(tx) != 0L]
    unique_tx_id <- unique(mcols(tx, level="within")[,"tx_id"])
    if(any(lengths(unique_tx_id) != 1L)){
        stop("'tx_id' must be a unique value within each element.",
             call. = FALSE)
    }
    tx
}

# shiftTranscriptToGenomic ----------------------------------------------------

.compareFUN_TX_TO_GENOME <- function(y, z, seqs){
    seqs[y][IRanges::IntegerList(as.list(z))]
}

.compareFUN_GENOME_TO_TX <- function(y, z, seqs){
    which(seqs[y] == z)
}

.resultFUN_TX_TO_GENOME <- function(tx_hits, subject_hits, starts, ends){
    # transfer gene information
    mcols(subject_hits)$sn_id <- ""
    mcols(subject_hits)$sn_id <-
        unlist(unique(mcols(tx_hits,level="within")[,"sn_id"]))
    mcols(subject_hits)$seq_start <- start(subject_hits)
    mcols(subject_hits)$seq_end <- end(subject_hits)
    mcols(subject_hits)$seq_strand <- strand(subject_hits)
    mcols(subject_hits)$seq_name <- seqnames(subject_hits)
    # prepare reformat of GRanges results
    seqnames <- unlist(unique(seqnames(tx_hits)),use.names = FALSE)
    strand <- unlist(unique(strand(tx_hits)),use.names = FALSE)
    mcols <- mcols(subject_hits)
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
    mcols(subject_hits)$sn_id <- ""
    mcols(subject_hits)$sn_id <-
        unlist(unique(mcols(tx_hits,level="within")[,"sn_id"]))
    mcols(subject_hits)$seq_start <- start(subject_hits)
    mcols(subject_hits)$seq_end <- end(subject_hits)
    mcols(subject_hits)$seq_strand <- strand(subject_hits)
    mcols(subject_hits)$seq_name <- seqnames(subject_hits)
    # reformat GRanges results
    seqnames <- mcols(subject_hits)$sn_id
    mcols <- mcols(subject_hits)
    ranges <- IRanges::IRanges(start = unlist(starts), end = unlist(ends))
    #
    ans <- GenomicRanges::GRanges(seqnames = seqnames,
                                  ranges = ranges,
                                  strand = "+",
                                  mcols)
    ans
}

.shift <- function(subject, tx, .compareFUN, .resultFUN){
    # match seqnames of subject and names of tx
    hits_coordinates <-
        suppressWarnings(GenomicRanges::findOverlaps(tx,subject))
    if(length(hits_coordinates) == 0L){
        stop("No 'subject'/'tx' matches found.",
             call. = FALSE)
    }
    # assemble results vectors
    qh <- queryHits(hits_coordinates)
    sh <- subjectHits(hits_coordinates)
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
    seqs <- positionSequence(tx)
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
        subject <- .norm_subject(subject)
        tx <- .norm_tx(tx)
        .shift(subject, tx, .compareFUN_TX_TO_GENOME, .resultFUN_TX_TO_GENOME)
    })

#' @rdname shiftGenomicToTranscript
#' @export
setMethod("shiftTranscriptToGenomic", c("GRangesList","GRangesList"),
    function(subject, tx){
        subject <- .norm_subject(subject)
        tx <- .norm_tx(tx)
        ans <- lapply(subject, .shift, tx, .compareFUN_TX_TO_GENOME,
                      .resultFUN_TX_TO_GENOME)
        GenomicRanges::GRangesList(ans)
    })


# shiftGenomicToTranscript ----------------------------------------------------

#' @rdname shiftGenomicToTranscript
#' @export
setMethod("shiftGenomicToTranscript", c("GRanges","GRangesList"),
    function(subject, tx){
        subject <- .norm_subject(subject)
        tx <- .norm_tx(tx)
        .shift(subject, tx, .compareFUN_GENOME_TO_TX, .resultFUN_GENOME_TO_TX)
    })

#' @rdname shiftGenomicToTranscript
#' @export
setMethod("shiftGenomicToTranscript", c("GRangesList","GRangesList"),
    function(subject, tx){
        subject <- .norm_subject(subject)
        tx <- .norm_tx(tx)
        ans <- lapply(subject, .shift, tx, .compareFUN_GENOME_TO_TX,
                      .resultFUN_GENOME_TO_TX)
        GenomicRanges::GRangesList(ans)
    })
