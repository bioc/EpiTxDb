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
#' @return a \code{\link[GenomicRanges:GRanges-class]{GRanges}} or
#'   \code{\link[GenomicRanges:GRangesList-class]{GRangesList}} object depending
#'   on the type of \code{subject}
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
    FUN <- S4Vectors::findMatches
    seqnames <- seqnames(subject)
    if(is(seqnames,"RleList")){
        seqnames <- as(seqnames,"CharacterList")
        hits <- suppressWarnings(lapply(seqnames,function(s){FUN(names(tx),s)}))
    } else {
        seqnames <- as.character(seqnames)
        hits <- suppressWarnings(FUN(names(tx), seqnames))
    }
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

.compareFUN_TX_TO_GENOME_List <- function(y, z, seqs){
    seqs[y][z]
}

.compareFUN_GENOME_TO_TX_List <- function(y, z, seqs){
    ans <- Map(
        function(a,b){
            match(b,seqs[[a]],nomatch = 0L)
        },
        y,
        z)
    IRanges::IntegerList(ans)
}

.resultFUN_TX_TO_GENOME <- function(tx_hits, subject_hits, starts, ends){
    # drop problematic subjects
    invalid_pos <- lengths(starts) < 1L | lengths(ends) < 1L
    if(any(invalid_pos)){
        warning("Dropping elements of 'subject' which could not be ",
                "repositioned.",
                call. = FALSE)
    }
    tx_hits <- tx_hits[!invalid_pos]
    subject_hits <- subject_hits[!invalid_pos]
    starts <- starts[!invalid_pos]
    ends <- ends[!invalid_pos]
    # transfer subject information
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
    non_dup <- !duplicated(paste0(as.character(ans),"_",mcols(ans)$mod_type,
                                  "_",mcols(ans)$seq_name))
    ans <- ans[non_dup]
    ans
}

.resultFUN_GENOME_TO_TX <- function(tx_hits, subject_hits, starts, ends){
    # drop problematic subjects
    invalid_pos <- lengths(starts) < 1L | lengths(ends) < 1L|
        strand(subject_hits) != unlist(unique(strand(tx_hits)))
    if(any(invalid_pos)){
        warning("Dropping elements of 'subject' which could not be ",
                "repositioned.",
                call. = FALSE)
    }
    tx_hits <- tx_hits[!invalid_pos]
    subject_hits <- subject_hits[!invalid_pos]
    starts <- starts[!invalid_pos]
    ends <- ends[!invalid_pos]
    # transfer subject information
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

.resultFUN_TX_TO_GENOME_List <- function(tx_hits, subject_hits, starts, ends){
    # drop problematic subjects
    invalid_pos <- starts == 0L | ends == 0L
    if(any(any(invalid_pos))){
        warning("Dropping elements of 'subject' which could not be ",
                "repositioned.",
                call. = FALSE)
    }
    tx_hits <- tx_hits[any(!invalid_pos)]
    subject_hits <- subject_hits[!invalid_pos]
    starts <- starts[!invalid_pos]
    ends <- ends[!invalid_pos]
    tx_hits <- tx_hits[lengths(tx_hits) > 0L]
    subject_hits <- subject_hits[lengths(subject_hits) > 0L]
    starts <- starts[lengths(starts) > 0L]
    ends <- ends[lengths(ends) > 0L]
    # transfer subject information
    mcols(subject_hits,level="within")[,"seq_start"] <- start(subject_hits)
    mcols(subject_hits,level="within")[,"seq_end"] <- end(subject_hits)
    mcols(subject_hits,level="within")[,"seq_strand"] <- strand(subject_hits)
    mcols(subject_hits,level="within")[,"seq_name"] <- seqnames(subject_hits)
    # prepare reformat of GRanges results
    seqnames <- rep(unlist(unique(seqnames(tx_hits)),use.names = FALSE),
                    lengths(subject_hits))
    strand <- rep(unlist(unique(strand(tx_hits)),use.names = FALSE),
                  lengths(subject_hits))
    mcols <- unlist(mcols(subject_hits,level="within"), use.names = FALSE)
    u_starts <- unlist(starts, use.names = FALSE)
    u_ends <- unlist(ends, use.names = FALSE)
    u_tmp <- u_starts
    u_starts[strand == "-"] <- u_ends[strand == "-"]
    u_ends[strand == "-"] <- u_tmp[strand == "-"]
    ranges <- IRanges::IRanges(start = u_starts, end = u_ends)
    # reformat GRanges results
    ans <- GenomicRanges::GRanges(seqnames = seqnames,
                                  ranges = ranges,
                                  strand = strand,
                                  mcols)
    ans <- relist(ans,subject_hits)
    names(ans) <- names(tx_hits)
    ans
}

.resultFUN_GENOME_TO_TX_List <- function(tx_hits, subject_hits, starts, ends){
    # drop problematic subjects
    invalid_pos <- starts == 0L | ends == 0L |
        unname(as(strand(subject_hits) != unique(strand(tx_hits)),"LogicalList"))
    if(any(any(invalid_pos))){
        warning("Dropping elements of 'subject' which could not be ",
                "repositioned.",
                call. = FALSE)
    }
    tx_hits <- tx_hits[any(!invalid_pos)]
    subject_hits <- subject_hits[!invalid_pos]
    starts <- starts[!invalid_pos]
    ends <- ends[!invalid_pos]
    tx_hits <- tx_hits[lengths(tx_hits) > 0L]
    subject_hits <- subject_hits[lengths(subject_hits) > 0L]
    starts <- starts[lengths(starts) > 0L]
    ends <- ends[lengths(ends) > 0L]
    # transfer subject information
    mcols(subject_hits,level="within")[,"seq_start"] <- start(subject_hits)
    mcols(subject_hits,level="within")[,"seq_end"] <- end(subject_hits)
    mcols(subject_hits,level="within")[,"seq_strand"] <- strand(subject_hits)
    mcols(subject_hits,level="within")[,"seq_name"] <- seqnames(subject_hits)
    # reformat GRanges results
    seqnames <- rep(names(tx_hits),lengths(subject_hits))
    mcols <- unlist(mcols(subject_hits,level="within"), use.names = FALSE)
    # switch starts ends depending on strand
    strand <- unlist(starts > ends, use.names = FALSE)
    u_starts <- unlist(starts, use.names = FALSE)
    u_ends <- unlist(ends, use.names = FALSE)
    u_tmp <- u_starts
    u_starts[strand] <- u_ends[strand]
    u_ends[strand] <- u_tmp[strand]
    ranges <- IRanges::IRanges(start = u_starts, end = u_ends)
    #
    ans <- GenomicRanges::GRanges(seqnames = seqnames,
                                  ranges = ranges,
                                  strand = "*",
                                  mcols)
    ans <- relist(ans,subject_hits)
    ans
}

.fixed_queryHits <- function(hits){
    if(is(hits,"list")){
        qh <- as(lapply(hits,queryHits),"IntegerList")
        qh <- unlist(unique(qh),use.names = FALSE)
    } else {
        qh <- queryHits(hits)
    }
    qh
}
.fixed_subjectHits <- function(hits){
    if(is(hits,"list")){
        sh <- as(lapply(hits,subjectHits),"IntegerList")

    } else {
        sh <- subjectHits(hits)
    }
    sh
}

.non_hits <- function(subject, sh){
    if((is(subject,"GRanges") | is(subject,"GRangesList")) & is.numeric(sh)){
        return(seq_along(subject) %in% sh)
    } else if(is(subject,"GRangesList") & is(sh,"IntegerList")) {
        return(names(subject) %in% names(sh))
    }
    stop(".")
}

.shift <- function(subject, tx, seqs, .hitsFUN, .compareFUN, .resultFUN){
    # match seqnames of subject and names of tx
    hits <- .hitsFUN(subject, tx)
    if(length(hits) == 0L){
        stop("No 'subject'/'tx' matches found.",
             call. = FALSE)
    }
    # assemble results vectors
    qh <- .fixed_queryHits(hits)
    sh <- .fixed_subjectHits(hits)
    f <- .non_hits(subject, sh)
    if(!all(f)){
        warning("Coordinates for ",sum(!f)," ranges of 'subject' not found: '",
                paste(unique(as.character(head(subject[!f], 10L))),
                      collapse = "','"),
                "'",ifelse(sum(!f) > 10L, " and more ..",""),".",
                call. = FALSE)
    }
    tx_hits <- tx[qh]
    subject_hits <- subject[sh]
    # reposition modification on transcript
    chunk_size <- 10^5
    chunks_n <- ceiling(length(qh) / chunk_size)
    f <- factor(unlist(Map(rep.int,seq_len(chunks_n),chunk_size)))
    f <- f[seq_along(qh)]
    # reposition starts and ends on sequence of positions
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
    ans <- .resultFUN(tx_hits, subject_hits, starts, ends)
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
    }
)

#' @rdname shiftGenomicToTranscript
#' @export
setMethod("shiftTranscriptToGenomic", c("GRangesList","GRangesList"),
    function(subject, tx){
        tx <- .norm_tx(tx)
        seqs <- positionSequence(tx)
        ans <- .shift(subject, tx, seqs, .hitsFUN_TX_TO_GENOME,
                      .compareFUN_TX_TO_GENOME_List,
                      .resultFUN_TX_TO_GENOME_List)
        ans
    }
)


# shiftGenomicToTranscript ----------------------------------------------------

#' @rdname shiftGenomicToTranscript
#' @export
setMethod("shiftGenomicToTranscript", c("GRanges","GRangesList"),
    function(subject, tx){
        tx <- .norm_tx(tx)
        seqs <- positionSequence(tx)
        .shift(subject, tx, seqs, .hitsFUN_GENOME_TO_TX,
               .compareFUN_GENOME_TO_TX,
               .resultFUN_GENOME_TO_TX)
    }
)

#' @rdname shiftGenomicToTranscript
#' @export
setMethod("shiftGenomicToTranscript", c("GRangesList","GRangesList"),
    function(subject, tx){
        tx <- .norm_tx(tx)
        seqs <- positionSequence(tx)
        ans <- .shift(subject, tx, seqs, .hitsFUN_GENOME_TO_TX,
                      .compareFUN_GENOME_TO_TX_List,
                      .resultFUN_GENOME_TO_TX_List)
        ans
    }
)
