#' @include EpiTxDb-class.R
#' @include makeEpiTxDb.R
NULL

EPITXDB_RMBASE_URL <- "http://rna.sysu.edu.cn/rmbase/download/"

#' @name makeEpiTxDbfromRMBase
#'
#' @title makeEpiTxDbfromRMBase
#'
#' @description
#' title
#'
#'
NULL

# makeEpiTxDbfromRMBase --------------------------------------------------------

.norm_tx <- function(tx){
    if(missing(tx) || !is(tx,"GRanges")){
        stop("'tx' must be an object of type 'TxDb' or 'GRanges'.")
    }
    if(!any("gene_name" %in% colnames(mcols(tx)))){
        stop("'tx' must contain a 'gene_name' column.")
    }
    tx
}

#' @importFrom BiocFileCache bfcquery bfcadd BiocFileCache
.download_RMBase_files <- function(organism, genome, type, files){
    urls <- paste0(EPITXDB_RMBASE_URL,organism,"/zip/",files)
    bfc <- BiocFileCache()
    rnames <- paste0("RMBase_",genome,"_",type)
    rnames_available <- rnames %in% bfcinfo(bfc)[,"rname"]
    files <- as.list(rep("",length(rnames)))
    files[rnames_available] <-
        Map(function(rname){
            bfcquery(bfc,rname)[,"rpath"]
        },
        rnames[rnames_available])
    files[!rnames_available] <-
        Map(function(url, rname){
            bfcadd(bfc,rname,url)
        },
        urls[!rnames_available],
        rnames[!rnames_available])
    unname(unlist(files))
}

#' @rdname makeEpiTxDbfromRMBase
#' @export
makeEpiTxDbfromRMBase <- function(organism, genome, type, tx,
                                  metadata = NULL, reassign.ids = FALSE){
    tx <- .norm_tx(tx)
    types <- listAvailableModFromRMBase(organism, genome)
    if(!(type %in% types)){
        stop("'type' must be a valid modification type from ",
             "listAvailableModFromRMBase() for the given 'organism' and ",
             "'genome'.")
    }
    files <- .get_RMBase_files(organism)
    f_genome <- vapply(strsplit(files,"_"),"[",character(1),2L) == genome
    f_mod <- vapply(strsplit(files[f_genome],"_"),"[",character(1),4L) == type
    files <- files[f_genome][f_mod]
    files <- .download_RMBase_files(organism, genome, type, files)
    makeEpiTxDbfromRMBaseFile(files, tx, metadata = metadata,
                              reassign.ids = reassign.ids)
}


.get_RMBase_header <- function(file){
    header <- readLines(file,7L)
    header
}

.get_RMBase_colnames <- function(file){
    colnames <- .get_RMBase_header(file)[7L]
    strsplit(trimws(gsub("#","",colnames)), "\\p{Zs}+", perl=TRUE)[[1L]]
}

.read_RMBase_file <- function(file){
    rmb <- read.table(file)
    colnames(rmb) <- .get_RMBase_colnames(file)
    rmb
}

.extract_GRanges_from_RMBase <- function(rmb){
    # extract position data
    pos <- rmb$modStart
    pos[rmb$strand == "-"] <- rmb$modEnd[rmb$strand == "-"]
    # extract modification information
    mod_id <- rmb$modId
    mod_type <- rmb$modType
    mod_name <- rmb$modName
    transcript_id <- as.integer(factor(rmb$geneName,unique(rmb$geneName)))
    # extract intermediate transcript data
    chromosome <- rmb$chromosome
    gene_name <- rmb$geneName
    # extract references
    reference_type <- pc(IRanges::CharacterList(as.list(rep("PMID",nrow(rmb)))),
                         IRanges::CharacterList(as.list(rep("RMBase_SUPPORT",nrow(rmb)))))
    reference <- pc(IRanges::CharacterList(as.list(rmb$pubmedIds)),
                    IRanges::CharacterList(as.list(rmb$supportList)))
    # assemble metadata columns
    mcols <- DataFrame(mod_id = mod_id,
                       mod_type = mod_type,
                       mod_name = mod_name,
                       transcript_id = transcript_id,
                       reference_type = reference_type,
                       reference = reference,
                       chromosome = chromosome,
                       gene_name = gene_name)
    # create GRanges result
    GenomicRanges::GRanges(seqnames = transcript_id,
                           ranges = IRanges::IRanges(pos, width = 1L),
                           strand = "*",
                           mcols)
}

.add_transcript_information <- function(gr, tx){
    msg <- "Couldn't associate transcript information  to results."
    if(missing(tx) || is.null(tx) || length(tx) == 0L){
        stop(msg, call. = FALSE)
    }
    tx <- tx[!is.na(mcols(tx)$tx_name) & !is.na(mcols(tx)$gene_name)]
    if(length(tx) == 0L){
        stop("No information available from the 'tx_name' or 'gene_name'.")
    }
    m_tx <- match(mcols(gr)$gene_name,mcols(tx)$gene_name)
    f_na <- is.na(m_tx)
    tx <- tx[m_tx[!f_na]]
    # drop problematic modifications
    if(any(!f_na)){
        warning("Dropping modifications which could not be associated with a ",
                "transcript ...", call. = FALSE)
    }
    gr <- gr[!f_na]
    mcols(gr)$transcript_ensembltrans <- mcols(tx)$transcript_ensembltrans
    mcols(gr)$transcript_name <- mcols(tx)$transcript_name
    pos <- start(gr)
    plus_strand <- as.logical(strand(tx) == "+")
    pos[plus_strand] <- pos[plus_strand] - start(tx[plus_strand]) + 1L
    pos[!plus_strand] <- end(tx[!plus_strand]) - pos[!plus_strand] + 1L
    # drop problematic modifications
    invalid_pos <- pos < 1L
    if(any(invalid_pos)){
        warning("Dropping modifications which could not be repositioned on the",
                " transcript (negative position) ...", call. = FALSE)
    }
    gr <- gr[!invalid_pos]
    pos <- pos[!invalid_pos]
    ranges <- IRanges::IRanges(pos, width = 1L)
    ranges(gr) <- ranges
    # reformat GRanges results
    mcols(gr)$mod_id <- as.integer(factor(mcols(gr)$mod_id,unique(mcols(gr)$mod_id)))
    GenomicRanges::GRanges(seqnames = mcols(gr)$transcript_id,
                           ranges = ranges(gr),
                           strand = strand(gr),
                           mcols(gr))
}

#' @rdname makeEpiTxDbfromRMBase
#' @export
makeEpiTxDbfromRMBaseFile <- function(files, tx, metadata = NULL,
                                      reassign.ids = FALSE){
    tx <- .norm_tx(tx)
    grl <- lapply(files,
                  function(file){
                      rmb <- .read_RMBase_file(file)
                      .extract_GRanges_from_RMBase(rmb)
                  })
    gr <- unlist(GRangesList(grl))
    gr <- .add_transcript_information(gr, txdb)
    makeEpiTxDbfromGRanges(gr, metadata = metadata, reassign.ids = reassign.ids)
}


#' @rdname makeEpiTxDbfromRMBase
#' @importFrom curl curl
#' @importFrom xml2 xml_find_all xml_attr
#' @export
listAvailableOrganismsFromRMBase <- function(){
    con <- curl::curl(EPITXDB_RMBASE_URL)
    page <- xml2::read_html(con)
    organisms <- xml2::xml_attr(xml2::xml_find_all(page,'//img[@alt="[DIR]"]//../following::a'),"href")
    organisms <- gsub("/","",organisms)
    organisms[!(organisms %in% c("ajax","otherspecies"))]
}


.get_RMBase_files <- function(organism){
    con <- curl::curl(paste0(EPITXDB_RMBASE_URL,organism,"/zip/"))
    page <- xml2::read_html(con)
    files <- xml2::xml_attr(xml2::xml_find_all(page,'//img[@alt="[   ]"]//../following::a'),"href")
    files[!grepl("^old",files)]
}

.get_RMBase_genomes <- function(files){
    unique(vapply(strsplit(files,"_"),"[",character(1),2L))
}

.listAvailableGenomesFromRMBase <- function(organism){
    files <- .get_RMBase_files(organism)
    .get_RMBase_genomes(files)
}

#' @rdname makeEpiTxDbfromRMBase
#' @importFrom curl curl
#' @export
listAvailableGenomesFromRMBase <- function(organism){
    if(!(organism %in% listAvailableOrganismsFromRMBase())){
        stop("'organism' must be a valid organism from ",
             "listAvailableOrganismsfromRMBase()")
    }
    .listAvailableGenomesFromRMBase(organism)
}

.get_RMBase_mod <- function(files, genome){
    f_genome <- vapply(strsplit(files,"_"),"[",character(1),2L) == genome
    unique(vapply(strsplit(files[f_genome],"_"),"[",character(1),4L))
}

#' @rdname makeEpiTxDbfromRMBase
#' @importFrom curl curl
#' @export
listAvailableModFromRMBase <- function(organism, genome){
    if(!(organism %in% listAvailableOrganismsFromRMBase())){
        stop("'organism' must be a valid organism from ",
             "listAvailableOrganismsfromRMBase()")
    }
    files <- .get_RMBase_files(organism)
    genomes <- .get_RMBase_genomes(files)
    if(!(genome %in% genomes)){
        stop("'genome' must be a valid genome for the fiven 'organism' from ",
             "listAvailableGenomesFromRMBase()")
    }
    .get_RMBase_mod(files,genome)
}