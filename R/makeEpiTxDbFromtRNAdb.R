#' @include EpiTxDb-class.R
#' @include makeEpiTxDb.R
NULL

#' @name makeEpiTxDbFromtRNAdb
#'
#' @title Create a \code{EpiTxDb} object from tRNAdb resources
#'
#' @description
#' \code{makeEpiTxDbFromtRNAdb} will make use of the tRNAdb online
#' resources.
#'
#' \code{makeEpiTxDbFromtRNAdb} uses the functions provided by the
#' \code{\link[tRNAdbImport:tRNAdbImport]{tRNAdbImport}} package.
#' \code{\link[tRNAdbImport:import.tRNAdb]{import.tRNAdb}} will be used with
#' \code{database = "RNA"} and the three different values for \code{origin}.
#'
#' @param organism A \code{character} value for an organism available on the
#'   tRNAdb website.
#' @param tx A \code{\link[GenomicRanges:GRangesList-class]{GRangesList}} object
#'   which will be used to shift the genomic coordinates to transcript
#'   coordinates. This is optional, but highly recommended (default: \code{tx =
#'   NULL}).
#' @param sequences A named \code{DNAStringSet} or \code{RNAStringSet}, which
#'   will be used to associate a tRNAdb result with a specific transcript.
#' @param dbURL The URL to the tRNA db website.
#' @param metadata See \code{\link[=makeEpiTxDb]{makeEpiTxDb}}
#'
#' @references
#' Jühling F, Mörl M, Hartmann RK, Sprinzl M, Stadler PF, Pütz J. 2009. "tRNAdb
#' 2009: compilation of tRNA sequences and tRNA genes." Nucleic Acids Research,
#' Volume 37 (suppl_1): D159–162. doi:10.1093/nar/gkn772.
#'
#' @export
#'
#' @examples
#'
NULL

# makeEpiTxDbFromtRNAdb --------------------------------------------------------

.add_tRNAdb_metadata_columns <- function(gr){
    # assemble metadata columns
    mcols <- mcols(gr)
    colnames(mcols) <- gsub("^mod$","mod_type",colnames(mcols))
    mcols$mod_name <- paste0(mcols$mod_type,"_",start(gr))
    mcols$mod_id <- seq_len(nrow(mcols))
    mcols$tx_id <- as.integer(seqnames(gr))
    mcols$tx_name <- as.character(seqnames(gr))
    mcols(gr) <- mcols
    gr
}

#' @rdname makeEpiTxDbFromtRNAdb
#' @importFrom tRNAdbImport import.tRNAdb
#' @export
gettRNAdbDataAsGRanges <- function(organism, tx = NULL, sequences = NULL,
                                   dbURL = tRNAdbImport::TRNA_DB_URL){
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
                                                     origin = origin,
                                                     dbURL = dbURL),
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
    mcols(mod)$ref_type <- ref_type[m]
    mcols(mod)$ref <- references[m]
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
        seq_rna_gr <- Biostrings::subseq(seq_rna_gr, 1L, width)
        # remove sequences which are definitly to long
        max_length <- max(lengths(seq_rna_gr))
        min_length <- min(lengths(seq_rna_gr))
        sequences[lengths(sequences) > (max_length + 50)] <-
            do.call(class(sequences),
                    list(paste0(rep("A",(max_length + 50)),collapse = "")))
        sequences[lengths(sequences) < min_length] <-
            do.call(class(sequences),
                    list(paste0(rep("A",(max_length + 50)),collapse = "")))
        hits <- Biostrings::vwhichPDict(sequences, seq_rna_gr,
                                        with.indels = TRUE, max.mismatch = 6)
        hits <- Hits(unlist(hits),as.integer(unlist(Map(rep.int,seq_along(hits),lengths(hits)))),
                     length(sequences), length(seq_rna_gr))
        hits_mod <- findMatches(seqnames(mod),names(seq_rna_gr))
        # transfer genenames to modifications
        tx_name <-
            unlist(unique(mcols(tx,level="within")[queryHits(hits),"tx_id"]))
        tx_name <- IRanges::CharacterList(split(tx_name,subjectHits(hits)))
        tx_name <- vapply(tx_name,paste,character(1),collapse=",")
        hits_mod <- hits_mod[subjectHits(hits_mod) %in% unique(subjectHits(hits))]
        mod[queryHits(hits_mod)]$tx_name <- tx_name[as.character(subjectHits(hits_mod))]
        tx_name <- IRanges::CharacterList(strsplit(as.character(mcols(mod)$tx_name),","))
        mod <- mod[unlist(Map(rep,seq_along(tx_name),lengths(tx_name)))]
        mcols(mod)$mod_id <- seq_along(mod)
        mod <- GenomicRanges::GRanges(seqnames = unlist(tx_name),
                                      ranges = IRanges::ranges(mod),
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
                                  metadata = NULL,
                                  dbURL = tRNAdbImport::TRNA_DB_UR){
    gr <- gettRNAdbDataAsGRanges(organism, tx = tx, sequences = sequences,
                                 dbURL = dbURL)
    if(!is.null(sequences)){
        gr <- gr[!duplicated(paste0(as.character(gr),"-",gr$mod_type))]
        colnames(mcols(gr)) <- gsub("mod_type","mod",colnames(mcols(gr)))
        gr <- Modstrings::removeIncompatibleModifications(gr, sequences)
        colnames(mcols(gr)) <- gsub("^mod$","mod_type",colnames(mcols(gr)))
    }
    makeEpiTxDbFromGRanges(gr, metadata = NULL, reassign.ids = FALSE)
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