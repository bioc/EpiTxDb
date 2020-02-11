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
makeEpiTxDbFromtRNAdb <- function(organism, tx = NULL, sequences = NULL,
                                  metadata = NULL){
    if(!assertive::is_a_non_empty_string(organism) ||
       !(organism %in% listAvailableOrganismsFromtRNAdb())){
        stop("'organism' must be a single character value and match an entry ",
             "from listAvailableOrganismsFromtRNAdb()")
    }

    # get tRNAdb information
    message("Downloading data from tRNAdb ...")
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
    message("Assembling data ...")
    seq_gr <- gr$tRNA_seq
    mod <- separate(seq_gr)
    # collect references
    m <- as.integer(S4Vectors::match(seqnames(mod), seqnames(gr)))
    reftype <- pc(IRanges::CharacterList(as.list(rep("tRNAdb_REF",length(gr))))[mcols(gr)$tRNAdb_reference != ""],
                  IRanges::CharacterList(as.list(rep("PMID",length(gr))))[mcols(gr)$tRNAdb_pmid != ""])
    references <- pc(mcols(gr)$tRNAdb_reference[mcols(gr)$tRNAdb_reference != ""],
                     mcols(gr)$tRNAdb_pmid[mcols(gr)$tRNAdb_pmid != ""])
    mcols(mod)$reference_type <- reftype[m]
    mcols(mod)$reference <- references[m]
    # assemble metadata columns
    mcols <- mcols(mod)
    colnames(mcols) <- gsub("^mod$","mod_type",colnames(mcols))
    mcols$mod_name <- paste0(mcols$mod_type,"_",start(mod))
    mcols$mod_id <- seq_len(nrow(mcols))
    mcols$transcript_id <- as.integer(seqnames(mod))
    mcols$transcript_name <- as.character(seqnames(mod))
    mcols(mod) <- mcols
    #
    if(!is.null(tx) && !is.null(sequences)){
        message("Trying to associate tRNAdb entries with transcripts ...")
        tx <- .norm_tx(tx)
        .check_tx_sequences(tx, sequences)
        seq_rna_gr <- as(seq_gr,"RNAStringSet")
        width <- nchar(seq_rna_gr)
        width[gr$tRNA_CCA.end] <- width[gr$tRNA_CCA.end] - 3L
        seq_rna_gr <- subseq(seq_rna_gr, 1L, width)
        hits <- findMatches(sequences, seq_rna_gr)
        hits_mod <- findMatches(seqnames(mod),names(seq_rna_gr))
        # transfer genenames to modifications
        transcript_name <- mcols(tx)[queryHits(hits),"transcript_name"]
        transcript_name <- IRanges::CharacterList(split(transcript_name,subjectHits(hits)))
        transcript_name <- vapply(transcript_name,paste,character(1),collapse=",")
        hits_mod <- hits_mod[subjectHits(hits_mod) %in% unique(subjectHits(hits))]
        mod[queryHits(hits_mod)]$transcript_name <- transcript_name[as.character(subjectHits(hits_mod))]
        transcript_name <- IRanges::CharacterList(strsplit(as.character(mcols(mod)$transcript_name),","))
        mod <- mod[unlist(Map(rep,seq_along(transcript_name),lengths(transcript_name)))]
        mcols(mod)$transcript_name <- unlist(transcript_name)
        mcols(mod)$mod_id <- seq_along(mod)
        mod <- GenomicRanges::GRanges(seqnames = mcols(mod)$transcript_name,
                                      ranges = ranges(mod),
                                      strand = strand(mod),
                                      mcols(mod))
    }
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