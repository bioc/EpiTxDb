#' @include TxModDb-class.R
#' @include makeTxModDb.R
NULL

#' @name makeTxModDbfromtRNAdb
#'
#' @title makeTxModDbfromtRNAdb
#'
#' @description
#' title
#'
#'
NULL

# makeTxModDbfromtRNAdb --------------------------------------------------------

#' @rdname makeTxModDbfromtRNAdb
#' @importFrom tRNAdbImport import.tRNAdb
#' @export
makeTxModDbfromtRNAdb <- function(organism, txdb, genome, metadata = NULL){
    if(!assertive::is_a_non_empty_string(organism) ||
       !(organism %in% listAvailableOrganismsfromtRNAdb())){
        stop("'organism' must be a single character value and match an entry ",
             "from listAvailableOrganismsfromtRNAdb()")
    }

    # get tRNAdb information
    gr <- lapply(c("allothers","plastid","mitochondrial"),
                 function(origin){
                     try(tRNAdbImport::import.tRNAdb(organism = organism,
                                                     database = "RNA",
                                                     origin = origin),
                         silent = TRUE)
                 })
    gr <- gr[!vapply(gr,is,logical(1),"try-error")]
    if(length(gr) == 0L){
        stop(".")
    }
    gr <- suppressWarnings(do.call(c,gr))
    mod <- separate(seq_gr)

    # assemble metadata columns
    mcols <- mcols(mod)
    colnames(mcols) <- gsub("^mod$","mod_type",colnames(mcols))
    mcols$mod_name <- paste0(mcols$mod_type,"_",start(mod))
    mcols$mod_id <- seq_len(nrow(mcols))
    mcols$transcript_id <- as.integer(seqnames(mod))
    mcols(mod) <- mcols
    # try to create transcript name and ensembltrans
    tx <- transcripts(txdb, filter = list(TXCHROM = seqnames(seqinfo(genome))))
    tx <- tx[width(tx) %in% width(gr)]
    seq <- getSeq(genome,tx)
    seqtype(seq) <- "RNA"
    seq_gr <- as(gr$tRNA_seq,"RNAStringSet")
    width <- width(seq_gr)
    width[gr$tRNA_CCA.end] <- width[gr$tRNA_CCA.end] - 3L
    seq_gr <- subseq(seq_gr, 1L, width)
    m_seq <- match(seq_gr,seq)
    f <- !is.na(m_seq)
    tx_ensembltrans <- tx$tx_name[m_seq[f]]
    m <- as.integer(match(seqnames(mod), seqnames(gr)[f]))
    mod$transcript_name <- ""
    mod$ensembltrans <- ""
    mcols(mod[which(seqnames(mod) %in% seqnames(gr)[f])])$ensembltrans <-
        tx_ensembltrans[m[!is.na(m)]]
    #
    makeTxModDbfromGRanges(mod, txdb, genome, metadata = NULL)
}


#' @rdname makeTxModDbfromtRNAdb
#' @importFrom httr modify_url POST content
#' @importFrom xml2 xml_find_all
#' @export
listAvailableOrganismsfromtRNAdb <- function(){
    res <- httr::POST(paste0(httr::modify_url(tRNAdbImport:::TRNA_DB_URL),
                             "DataOutput/Organisms"))
    html <- httr::content(res)
    organsisms <- as.character(xml2::xml_find_all(html,'.//div[@id="inhalt"]//span[count(.//span)=0]//a//b/text()'))
    rna_seq <- as.character(xml2::xml_find_all(html,'.//div[@id="inhalt"]//span[count(.//span)=0]//a//small/text()'))
    rna_n <- vapply(regmatches(rna_seq,regexec("genes,[ ]{1}([0-9]+)[ ]{1}RNA",rna_seq)),"[",character(1),2)
    organsisms[rna_n != 0]
}