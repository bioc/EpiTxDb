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

.add_tRNAdb_metadata_columns <- function(gr){
    # assemble metadata columns
    mcols <- mcols(gr)
    colnames(mcols) <- gsub("^mod$","mod_type",colnames(mcols))
    mcols$mod_name <- paste0(mcols$mod_type,"_",start(gr))
    mcols$mod_id <- seq_len(nrow(mcols))
    mcols$transcript_id <- as.integer(seqnames(gr))
    mcols$transcript_name <- as.character(seqnames(gr))
    mcols(gr) <- mcols
    gr
}

#' @rdname makeEpiTxDbFromtRNAdb
#' @importFrom tRNAdbImport import.tRNAdb
#' @export
gettRNAdbDataAsGRanges <- function(organism, tx = NULL, sequences = NULL){
    if(!assertive::is_a_non_empty_string(organism) ||
       !(organism %in% listAvailableOrganismsFromtRNAdb())){
        stop("'organism' must be a single character value and match an entry ",
             "from listAvailableOrganismsFromtRNAdb()")
    }

    # get tRNAdb information
    message("Loading data from tRNAdb ...")
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
    f_ref <- mcols(gr)$tRNAdb_reference != ""
    f_pmid <- mcols(gr)$tRNAdb_pmid != ""

    ref_type <- list(IRanges::CharacterList(as.list(rep("tRNAdb_ID",
                                                        length(gr)))),
                     IRanges::CharacterList(as.list(rep("tRNAdb_REF",
                                                        length(gr)))),
                     IRanges::CharacterList(as.list(rep("PMID",
                                                        length(gr)))))
    ref_type[[2]][!f_ref] <- ""
    ref_type[[3]][!f_pmid] <- ""
    ref_type <- do.call(pc,ref_type)
    ref_type <- ref_type[ref_type != ""]

    references <- list(mcols(gr)$tRNAdb_ID,
                      mcols(gr)$tRNAdb_reference,
                      mcols(gr)$tRNAdb_pmid)
    references[[2]][!f_ref] <- ""
    references[[3]][!f_pmid] <- ""
    references <- do.call(pc,references)
    references <- references[references != ""]

    m <- as.integer(S4Vectors::match(seqnames(mod), seqnames(gr)))
    mcols(mod)$reference_type <- ref_type[m]
    mcols(mod)$reference <- references[m]
    # assemble metadata columns
    mod <- .add_tRNAdb_metadata_columns(mod)
    #
    if(!is.null(tx) && !is.null(sequences)){
        message("Trying to associate tRNAdb entries with transcripts ...")
        tx <- .norm_tx(tx)
        .check_tx_sequences(tx, sequences)
        seq_rna_gr <- as(seq_gr,"RNAStringSet")
        width <- nchar(seq_rna_gr)
        width[gr$tRNA_CCA.end] <- width[gr$tRNA_CCA.end] - 3L
        seq_rna_gr <- subseq(seq_rna_gr, 1L, width)
        # remove sequences which are definitly to long
        max_length <- max(lengths(seq_rna_gr))
        min_length <- min(lengths(seq_rna_gr))
        sequences[lengths(sequences) > (max_length + 50)] <-
            do.call(class(sequences),
                    list(paste0(rep("A",(max_length + 50)),collapse = "")))
        sequences[lengths(sequences) < min_length] <-
            do.call(class(sequences),
                    list(paste0(rep("A",(max_length + 50)),collapse = "")))
        hits <- vwhichPDict(sequences, seq_rna_gr, with.indels = TRUE,
                            max.mismatch = 6)
        hits <- Hits(unlist(hits),as.integer(unlist(Map(rep.int,seq_along(hits),lengths(hits)))),
                     length(sequences), length(seq_rna_gr))
        hits_mod <- findMatches(seqnames(mod),names(seq_rna_gr))
        # transfer genenames to modifications
        transcript_name <-
            unlist(unique(mcols(tx,level="within")[queryHits(hits),"tx_id"]))
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
        # assemble metadata columns
        mod <- .add_tRNAdb_metadata_columns(mod)
    }
    mod
}

#' @rdname makeEpiTxDbFromtRNAdb
#' @export
makeEpiTxDbFromtRNAdb <- function(organism, tx = NULL, sequences = NULL,
                                  metadata = NULL){
    gr <- gettRNAdbDataAsGRanges(organism, tx = tx, sequences = sequences)
    if(!is.null(sequences)){
        gr <- gr[!duplicated(paste0(as.character(gr),"-",gr$mod_type))]
        colnames(mcols(gr)) <- gsub("mod_type","mod",colnames(mcols(gr)))
        gr <- Modstrings::removeIncompatibleModifications(gr, sequences)
        colnames(mcols(gr)) <- gsub("^mod$","mod_type",colnames(mcols(gr)))
    }
    makeEpiTxDbfromGRanges(gr, metadata = NULL)
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