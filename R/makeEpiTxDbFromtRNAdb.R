#' @include EpiTxDb-class.R
#' @include makeEpiTxDb.R
NULL

#' @name makeEpiTxDbFromtRNAdb
#'
#' @title makeEpiTxDbFromtRNAdb
#'
#' @description
#' title
#'
#'
NULL

# makeEpiTxDbFromtRNAdb --------------------------------------------------------

#' @rdname makeEpiTxDbFromtRNAdb
#' @importFrom tRNAdbImport import.tRNAdb
#' @export
makeEpiTxDbFromtRNAdb <- function(organism, txdb = NULL, sequences = NULL,
                                  metadata = NULL){
    if(!assertive::is_a_non_empty_string(organism) ||
       !(organism %in% listAvailableOrganismsFromtRNAdb())){
        stop("'organism' must be a single character value and match an entry ",
             "from listAvailableOrganismsFromtRNAdb()")
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
    seq_gr <- gr$tRNA_seq
    mod <- separate(seq_gr)

    # assemble metadata columns
    mcols <- mcols(mod)
    colnames(mcols) <- gsub("^mod$","mod_type",colnames(mcols))
    mcols$mod_name <- paste0(mcols$mod_type,"_",start(mod))
    mcols$mod_id <- seq_len(nrow(mcols))
    mcols$transcript_id <- as.integer(seqnames(mod))
    mcols(mod) <- mcols
    # try to create transcript name and ensembltrans
    if(!is.null(txdb) && !is.null(genome)){
        tx <- transcripts(txdb,
                          filter = list(TXCHROM = seqnames(seqinfo(sequences))))
        tx <- tx[width(tx) %in% width(gr)]
        seq <- getSeq(sequences,tx)
        seqtype(seq) <- "RNA"
        seq_gr <- as(gr$tRNA_seq,"RNAStringSet")
        width <- width(seq_gr)
        width[gr$tRNA_CCA.end] <- width[gr$tRNA_CCA.end] - 3L
        seq_gr <- subseq(seq_gr, 1L, width)
        m_seq <- S4Vectors::match(seq_gr,seq)
        f <- !is.na(m_seq)
        tx_ensembltrans <- tx$tx_name[m_seq[f]]
        m <- as.integer(S4Vectors::match(seqnames(mod), seqnames(gr)[f]))
        mod$transcript_name <- ""
        mod$ensembltrans <- ""
        mcols(mod[which(seqnames(mod) %in% seqnames(gr)[f])])$ensembltrans <-
            tx_ensembltrans[m[!is.na(m)]]
    }
    # collect references
    m <- as.integer(S4Vectors::match(seqnames(mod), seqnames(gr)))
    reftype <- pc(IRanges::CharacterList(as.list(rep("tRNAdb_REF",length(gr))))[mcols(gr)$tRNAdb_reference != ""],
                  IRanges::CharacterList(as.list(rep("PMID",length(gr))))[mcols(gr)$tRNAdb_pmid != ""])
    references <- pc(mcols(gr)$tRNAdb_reference[mcols(gr)$tRNAdb_reference != ""],
                     mcols(gr)$tRNAdb_pmid[mcols(gr)$tRNAdb_pmid != ""])
    mcols(mod)$reference_type <- reftype[m]
    mcols(mod)$reference <- references[m]
    #
    makeEpiTxDbfromGRanges(mod, metadata = NULL)
}


#' @rdname makeEpiTxDbFromtRNAdb
#' @importFrom httr modify_url POST content
#' @importFrom xml2 xml_find_all
#' @export
listAvailableOrganismsFromtRNAdb <- function(){
    res <- httr::POST(paste0(httr::modify_url(tRNAdbImport:::TRNA_DB_URL),
                             "DataOutput/Organisms"))
    html <- httr::content(res)
    organsisms <- as.character(xml2::xml_find_all(html,'.//div[@id="inhalt"]//span[count(.//span)=0]//a//b/text()'))
    rna_seq <- as.character(xml2::xml_find_all(html,'.//div[@id="inhalt"]//span[count(.//span)=0]//a//small/text()'))
    rna_n <- vapply(regmatches(rna_seq,regexec("genes,[ ]{1}([0-9]+)[ ]{1}RNA",rna_seq)),"[",character(1),2)
    organsisms[rna_n != 0]
}