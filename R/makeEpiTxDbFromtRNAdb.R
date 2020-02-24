#' @include EpiTxDb-class.R
#' @include makeEpiTxDb.R
NULL

#' @name makeEpiTxDbFromtRNAdb
#'
#' @title Create a \code{EpiTxDb} object from tRNAdb resources
#'
#' @description
#' \code{makeEpiTxDbFromtRNAdb} will make use of the tRNAdb online
#' resources and extract the modification information from the RNA database.
#'
#' When a
#'
#'
#'
#' \code{makeEpiTxDbFromtRNAdb} uses the functions provided by the
#' \code{\link[tRNAdbImport:tRNAdbImport]{tRNAdbImport}} package.
#' \code{\link[tRNAdbImport:import.tRNAdb]{import.tRNAdb}} will be used with
#' \code{database = "RNA"} and the three different values for \code{origin}.
#'
#' @param organism A \code{character} value for an organism available on the
#'   tRNAdb website.
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

.add_reference_information_tRNAdb <- function(mod, gr){
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
    mod
}

.map_modifications_to_sequences <- function(mod, gr, sequences){
    seq <- as(gr$tRNA_seq,"RNAStringSet")
    width <- nchar(seq)
    width[gr$tRNA_CCA.end] <- width[gr$tRNA_CCA.end] - 3L
    seq <- Biostrings::subseq(seq, 1L, width)
    # remove sequences which are definitly to long
    max_length <- max(lengths(seq))
    min_length <- min(lengths(seq))
    sequences[lengths(sequences) > (max_length + 50)] <-
        do.call(class(sequences),
                list(paste0(rep("A",(max_length + 50)),collapse = "")))
    sequences[lengths(sequences) < min_length] <-
        do.call(class(sequences),
                list(paste0(rep("A",(max_length + 50)),collapse = "")))
    hits <- Biostrings::vwhichPDict(sequences, seq,
                                    with.indels = TRUE, max.mismatch = 5L)
    hits <- Hits(unlist(hits),as.integer(unlist(Map(rep.int,seq_along(hits),lengths(hits)))),
                 length(sequences), length(seq))
    hits_mod <- findMatches(seqnames(mod),names(seq))
    # transfer names of sequences to seqnames for modifications and expand the
    # result with multiple matches
    sn_name <- names(sequences)[queryHits(hits)]
    sn_name <- IRanges::CharacterList(split(sn_name,subjectHits(hits)))
    sn_name <- vapply(sn_name,paste,character(1),collapse=",")
    hits_mod <- hits_mod[subjectHits(hits_mod) %in% unique(subjectHits(hits))]
    mcols(mod)[queryHits(hits_mod),"sn_name"] <- sn_name[as.character(subjectHits(hits_mod))]
    sn_name <- IRanges::CharacterList(strsplit(as.character(mcols(mod)$sn_name),","))
    # expand results
    mod <- mod[unlist(Map(rep,seq_along(sn_name),lengths(sn_name)))]
    # assemble new result
    mcols <- mcols(mod)[,colnames(mcols(mod)) != "sn_name"]
    mod <- GenomicRanges::GRanges(seqnames = unlist(sn_name),
                                  ranges = IRanges::ranges(mod),
                                  strand = strand(mod),
                                  mcols)
    mcols(mod)$mod_id <- seq_along(mod)
    mod
}

.add_tRNAdb_metadata_columns <- function(gr){
    # assemble metadata columns
    mcols <- mcols(gr)
    colnames(mcols) <- gsub("^mod$","mod_type",colnames(mcols))
    mcols$mod_name <- paste0(mcols$mod_type,"_",start(gr),"_",seqnames(gr))
    mcols$mod_id <- seq_len(nrow(mcols))
    mcols(gr) <- mcols
    gr
}

#' @rdname makeEpiTxDbFromtRNAdb
#' @importFrom tRNAdbImport import.tRNAdb
#' @export
gettRNAdbDataAsGRanges <- function(organism, sequences = NULL,
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
    mod <- separate(gr$tRNA_seq)
    mod <- .add_reference_information_tRNAdb(mod, gr)
    #
    if(!is.null(sequences)){
        message("Trying to associate tRNAdb entries with sequences ...")
        if(is.null(names(sequences)) ||
           anyDuplicated(names(sequences)) ||
           any(names(sequences) == "")){
            stop("names() of 'sequences' must be set, unique and not empty.",
                 call. = FALSE)
        }
        mod <- .map_modifications_to_sequences(mod, gr, sequences)
    }
    # assemble metadata columns
    mod <- .add_tRNAdb_metadata_columns(mod)
    mod
}

#' @rdname makeEpiTxDbFromtRNAdb
#' @export
makeEpiTxDbFromtRNAdb <- function(organism, sequences = NULL, metadata = NULL,
                                  dbURL = tRNAdbImport::TRNA_DB_URL){
    gr <- gettRNAdbDataAsGRanges(organism, sequences = sequences, dbURL = dbURL)
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