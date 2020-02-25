#' @include EpiTxDb-class.R
#' @include makeEpiTxDb.R
#' @include shiftGenomicToTranscript.R
NULL

#' @name makeEpiTxDbFromRMBase
#'
#' @title Create a \code{EpiTxDb} object from RMBase v2.0 online resources
#'
#' @description
#' \code{makeEpiTxDbFromRMBase} will make use of the RMBase v2.0 online
#' resources.
#'
#' @param organism A \code{character} value, which must match an organism
#'   descriptor on the RMBase download website.
#' @param genome A \code{character} value, which must match a genome
#'   descriptor on the RMBase download website.
#' @param type A \code{character} value, which must match one or more
#'   modification descriptors on the RMBase download website.
#' @param files From \code{organism}, \code{genome} and \code{type} the
#'   available files will be downloaded using the
#'   \code{\link[BiocFileCache:BiocFileCache-class]{BiocFileCache}} interface
#'   and passed on to \code{makeEpiTxDbFromRMBaseFiles}. However, individual
#'   files can be provided as well.
#' @param tx A \code{\link[GenomicRanges:GRangesList-class]{GRangesList}} object
#'   which will be used to shift the genomic coordinates to transcript
#'   coordinates. This is optional, but highly recommended. (default:
#'   \code{tx = NULL}).
#' @param sequences A named \code{DNAStringSet} or \code{RNAStringSet}, which
#'   will be used to check whether the defined modifications are compatible with
#'   the original base. This uses
#'   \code{\link[Modstrings:separate]{removeIncompatibleModifications()}}
#'   function from the \code{Modstrings} package.
#' @param metadata,reassign.ids See \code{\link[=makeEpiTxDb]{makeEpiTxDb}}
#'
#' @export
NULL

#' @rdname makeEpiTxDbFromRMBase
#' @export
EPITXDB_RMBASE_URL <- "http://rna.sysu.edu.cn/rmbase/download/"

# makeEpiTxDbFromRMBase --------------------------------------------------------

.get_RMBase_rnames <- function(organism, genome, type){
    paste0("RMBase_",organism,"_",genome,"_",type)
}

.check_RMBase_files_available <- function(bfc, organism, genome, type){
    # get BiocFileCache information
    rnames <- .get_RMBase_rnames(organism, genome, type)
    bfci <- BiocFileCache::bfcinfo(bfc)
    m <- match(rnames,bfci$rname)
    m <- m[!is.na(m)]
    res <- bfci[m,]
    #
    nrow(res) == length(rnames)
}

.get_RMBase_files_available <- function(bfc, organism, genome, type){
    # get BiocFileCache information
    rnames <- .get_RMBase_rnames(organism, genome, type)
    bfci <- BiocFileCache::bfcinfo(bfc)
    m <- match(rnames,bfci$rname)
    m <- m[!is.na(m)]
    res <- bfci[m,]
    #
    files <- as.list(res$rpath)
    if(length(rnames) != length(files)){
        stop(".")
    }
    #
    m <- match(bfci$rname, rnames)
    m <- m[!is.na(m)]
    # update resource if expired
    check_update_required <- BiocFileCache::bfcneedsupdate(bfc, res$rid)
    check_update_required <- check_update_required &
        !is.na(check_update_required)
    files[check_update_required] <-
        Map(function(rid){
            BiocFileCache::bfcdownload(bfc,rid)
        },
        res$rid[check_update_required])
    unname(unlist(files))
}

#' @importFrom BiocFileCache bfcquery bfcadd BiocFileCache
.download_RMBase_files <- function(bfc, organism, genome, type){
    # get file names
    files <- .get_RMBase_files(organism)
    f_genome <- vapply(strsplit(files,"_"),"[",character(1),2L) == genome
    f_mod <- vapply(strsplit(files[f_genome],"_"),"[",character(1),4L) %in% type
    files <- files[f_genome][f_mod]
    urls <- paste0(EPITXDB_RMBASE_URL,organism,"/zip/",files)
    #
    rnames <- .get_RMBase_rnames(organism, genome, type)
    if(length(rnames) != length(urls)){
        stop(".")
    }
    # get BiocFileCache info
    bfci <- BiocFileCache::bfcinfo(bfc)
    # check which rnames are available
    rnames_available <- rnames %in% bfci$rname
    # get file name if available
    files <- as.list(rep("",length(rnames)))
    files[rnames_available] <-
        Map(function(rname){
            BiocFileCache::bfcquery(bfc,rname)[,"rpath"]
        },
        rnames[rnames_available])
    # download if not available and return file name
    files[!rnames_available] <-
        Map(function(url, rname){
            BiocFileCache::bfcadd(bfc,rname,url)
        },
        urls[!rnames_available],
        rnames[!rnames_available])
    unname(unlist(files))
}

#' @rdname makeEpiTxDbFromRMBase
#' @export
downloadRMBaseFiles <- function(organism, genome, type){
    bfc <- BiocFileCache::BiocFileCache()
    if(!.check_RMBase_files_available(bfc, organism, genome, type)){
        types <- listAvailableModFromRMBase(organism, genome)
        if(!all(type %in% types)){
            stop("'type' must be a valid modification type from ",
                 "listAvailableModFromRMBase() for the given 'organism' and ",
                 "'genome'.")
        }
        files <- .download_RMBase_files(bfc, organism, genome, type)
    } else {
        files <- .get_RMBase_files_available(bfc, organism, genome, type)
    }
    files
}

#' @rdname makeEpiTxDbFromRMBase
#' @export
makeEpiTxDbFromRMBase <- function(organism, genome, type, tx = NULL,
                                  sequences = NULL, metadata = NULL,
                                  reassign.ids = FALSE){
    message("Loading RMBase files ...")
    files <- downloadRMBaseFiles(organism, genome, type)
    makeEpiTxDbFromRMBaseFiles(files, tx = tx, sequences = sequences,
                               metadata = metadata, reassign.ids = reassign.ids)
}

# makeEpiTxDbFromRMBaseFiles ---------------------------------------------------

.get_RMBase_header <- function(file){
    header <- readLines(file,7L)
    header
}

.get_RMBase_colnames <- function(file){
    colnames <- .get_RMBase_header(file)[7L]
    strsplit(trimws(gsub("#","",colnames)), "\\p{Zs}+|\\t", perl=TRUE)[[1L]]
}

EPITXDB_RMBASE_REQ_COLUMS <- c("chromosome", "modStart", "modEnd", "modId",
                               "score", "strand", "modName", "modType",
                               "supportList", "pubmedIds", "geneName",
                               "sequence")

#' @importFrom utils read.delim
.read_RMBase_file <- function(file){
    rmb <- try(read.delim(file, skip = 7L), silent = TRUE)
    if(is(rmb,"try-error")){
        stop("Malformated input file: ", file,".\nError message: ",
             as.character(rmb),
             call. = FALSE)
    }
    if(ncol(rmb) < 15L){
        stop(".")
    }
    colnames <- .get_RMBase_colnames(file)
    if(length(colnames) < ncol(rmb)){
        stop("Header of input could not read correctly. Mismatching number of ",
             "data columns and column labels.",
             call. = FALSE)
    }
    colnames(rmb) <- colnames[seq_len(ncol(rmb))]
    f <- EPITXDB_RMBASE_REQ_COLUMS %in% colnames(rmb)
    if(!all(f)){
        stop("Missing required column from RMBase data. (File: ",
             file, " / Missing columns: ",
             EPITXDB_RMBASE_REQ_COLUMS[!(f)],")",
             call. = FALSE)
    }
    rmb <- rmb[,colnames(rmb) %in% EPITXDB_RMBASE_REQ_COLUMS]
    rmb
}

# if a prefix in the chromosome identifier is present, try to remove it
# only if it simplifies the result
.simplify_chromosome_identifiers <- function(chromosome, seqlevels = NA){
    chromosome <- as.character(chromosome)
    # fix some weird mito chromosome annotation
    f <- chromosome %in% "Mito"
    if(any(f)){
        chromosome[f] <- gsub("Mito","M",chromosome[f])
    }
    # mismatching chromosomes
    if(all(!is.na(seqlevels))){
        mm_chr <- !(chromosome %in% seqlevels)
    } else {
        return(factor(chromosome,unique(chromosome)))
    }
    f_chr_remove <- grepl("Chr|chr",chromosome[mm_chr])
    if(any(f_chr_remove)){
        if(all(levels(factor(gsub("Chr|chr","",chromosome[mm_chr][f_chr_remove]))) %in%
               levels(chromosome[mm_chr][f_chr_remove])) ||
           all(levels(chromosome[mm_chr][f_chr_remove]) %in%
               levels(factor(gsub("Chr|chr","",chromosome[mm_chr][f_chr_remove]))))){
            tmp <- gsub("Chr|chr","",chromosome[mm_chr])
            chromosome[mm_chr] <- tmp
        }
    } else {
        chr_add1 <- paste0("chr",chromosome[mm_chr])
        if(all(chr_add1 %in% seqlevels)){
            chromosome[mm_chr] <- chr_add1
        } else {
            chr_add2 <- paste0("Chr",chromosome[mm_chr])
            if(all(chr_add2 %in% seqlevels)){
                chromosome[mm_chr] <- chr_add2
            }
        }
    }
    factor(chromosome,unique(chromosome))
}

.fix_non_standard_mod_types <- function(mod_type){
    mod_type <- gsub("m22","m2,2",mod_type)
    mod_type <- gsub("m42","m4,4",mod_type)
    mod_type <- gsub("m62","m6,6",mod_type)
    mod_type <- gsub("Tm","Um",mod_type)
    if("m5C" %in% unique(mod_type)){
        mod_type <- gsub("N","m5C",mod_type)
    }
    mod_type
}

#' @importFrom Biostrings DNAStringSet subseq
.extract_GRanges_from_RMBase <- function(rmb, seqtype = "RNA"){
    ############################################################################
    ### check modification information on correct base
    seq <- Biostrings::DNAStringSet(rmb$sequence)
    if(unique(width(seq)) != 41L){
        stop(".")
    }
    seqtype(seq) <- seqtype
    seqtype <- paste0("Mod",seqtype(seq))
    # get the type of modification annotation
    mod_type <- .fix_non_standard_mod_types(as.character(rmb$modType))
    nc_type <- Modstrings:::.get_nc_type(mod_type, seqtype)
    if(nc_type != "short"){
        stop(".")
    }
    # get original base
    codec <- Modstrings:::modscodec(seqtype)
    modValues <- Modstrings:::.norm_seqtype_modtype(mod_type, seqtype,
                                                    nc_type)
    f <- Modstrings:::values(codec)[match(modValues,
                                          Modstrings:::values(codec))]
    # check if reported base matches original base
    subseq <- Biostrings::subseq(seq,21L,21L)
    base_mm <- Modstrings:::originatingBase(codec)[f] != subseq
    # if not delete the modifications with a warning
    if(any(base_mm)){
        warning("Detected mismatch of modification and originating base in ",
                "RMBase data for '", paste(unique(mod_type), collapse = "', '"),
                "'. Removing them ... ",
                call. = FALSE)
        rmb <- rmb[!base_mm,]
    }
    ############################################################################
    # extract position data
    strand <- rmb$strand
    # pos <- rmb$modStart
    # pos[strand == "-"] <- rmb$modEnd[strand == "-"]
    pos <- rmb$modEnd
    # extract modification information
    mod_id <- rmb$modId
    mod_type <- .fix_non_standard_mod_types(as.character(rmb$modType))
    mod_name <- rmb$modName
    # extract intermediate transcript data
    chromosome <- rmb$chromosome
    gene_name <- rmb$geneName
    if(any(is.na(gene_name)) || any(is.null(gene_name))){
        stop(".")
    }
    # extract references
    f_pmid <- rmb$pubmedIds != "" & rmb$pubmedIds != "-"
    f_ref <- rmb$supportList != "" & rmb$supportList != "-"
    ref_type <- list(IRanges::CharacterList(as.list(rep("RMBASE_MOD_ID",
                                                        nrow(rmb)))),
                     IRanges::CharacterList(as.list(rep("PMID",
                                                        nrow(rmb)))),
                     IRanges::CharacterList(as.list(rep("RMBase_REF",
                                                        nrow(rmb)))))
    ref_type[[2]][!f_pmid] <- ""
    ref_type[[3]][!f_ref] <- ""
    ref_type <- do.call(pc,ref_type)
    ref_type <- ref_type[ref_type != ""]
    reference <- list(relist(as.character(rmb$modId),
                             IRanges::PartitioningByWidth(lengths(rmb$modId))),
                      relist(as.character(rmb$pubmedIds),
                             IRanges::PartitioningByWidth(lengths(rmb$pubmedIds))),
                      relist(as.character(rmb$supportList),
                             IRanges::PartitioningByWidth(lengths(rmb$supportList))))
    reference[[2]][!f_pmid] <- ""
    reference[[3]][!f_ref] <- ""
    reference <- do.call(pc,reference)
    reference <- reference[reference != ""]
    # assemble metadata columns
    gene_name <-
        IRanges::CharacterList(strsplit(as.character(gene_name),","))
    mcols <- DataFrame(mod_id = mod_id,
                       mod_type = mod_type,
                       mod_name = mod_name,
                       chromosome = chromosome,
                       gene_name = gene_name,
                       score = rmb$score,
                       ref_type = ref_type,
                       ref = reference)
    # create GRanges result
    ans <- GenomicRanges::GRanges(seqnames = chromosome,
                                  ranges = IRanges::IRanges(pos, width = 1L),
                                  strand = strand,
                                  mcols)
    #
    ans
}

.get_RMBase_data <- function(files){
    grl <- lapply(files,
                  function(file){
                      rmb <- .read_RMBase_file(file)
                      .extract_GRanges_from_RMBase(rmb)
                  })
    grl
}

.gene_names_to_chromosome <- function(gr){
    chromosome <- mcols(gr)$gene_name
    GenomicRanges::GRanges(seqnames = chromosome,
                           ranges = IRanges::ranges(gr),
                           strand = strand(gr),
                           mcols(gr))
}



.check_tx_sequences <- function(tx, sequences){
    if(is.null(sequences)){
        return()
    }
    if(!is(sequences,"DNAStringSet") & !is(sequences,"RNAStringSet")){
        stop("'sequences' must be a DNA/RNAStringSet.", call. = FALSE)
    }
    if(is.null(names(sequences)) || !all(names(sequences) == names(tx) )){
        stop("'sequences' must be named and names have to match 'tx'",
             call. = FALSE)
    }
    if(!all(sum(width(tx)) == width(sequences))){
        stop("'sequences' and 'tx' must have the same dimensions, e.g. the ",
             "length of sequences and annotation must match.",
             call. = FALSE)
    }
}

#' @rdname makeEpiTxDbFromRMBase
#' @export
getRMBaseDataAsGRanges <- function(files){
    message("Assembling RMBase data ...")
    # getting raw data from RMBase files
    grl <- .get_RMBase_data(files)
    gr <- unlist(GenomicRanges::GRangesList(grl))
    gr
}

.remove_incompatible_mod_by_seq <- function(gr, seq){
    seq <- as(seq,"RNAStringSet")
    colnames(mcols(gr)) <- gsub("mod_type","mod",colnames(mcols(gr)))
    # do plus and minus strand separatly, since removeIncompatibleModifications
    # only accepts plus strand
    gr_plus <-
        Modstrings::removeIncompatibleModifications(gr[strand(gr) == "+"], seq)
    gr_minus <- gr[strand(gr) == "-"]
    strand(gr_minus) <- "+"
    gr_minus <-
        Modstrings::removeIncompatibleModifications(gr_minus,
                                                    Biostrings::complement(seq))
    strand(gr_minus) <- "-"
    gr <- c(gr_plus,gr_minus)
    gr <- gr[order(seqnames(gr),start(gr),strand(gr))]
    colnames(mcols(gr)) <- gsub("^mod$","mod_type",colnames(mcols(gr)))
    gr
}

.add_sequence_check_to_metadata <- function(metadata){
    metadata <- rbind(metadata,
                      data.frame(name = "Checked against sequence",
                                 value = "Yes",
                                 check.names = FALSE,
                                 stringsAsFactors = FALSE))
    metadata
}

#' @rdname makeEpiTxDbFromRMBase
#' @export
makeEpiTxDbFromRMBaseFiles <- function(files, tx = NULL, sequences = NULL,
                                       metadata = NULL, reassign.ids = FALSE){
    gr <- getRMBaseDataAsGRanges(files)
    if(!is.null(tx)){
        sl <- GenomeInfoDb::seqlevels(tx)
    } else if(!is.null(sequences)) {
        sl <- names(sequences)
    } else {
        sl <- GenomeInfoDb::seqlevels(gr)
    }
    chromosome <- .simplify_chromosome_identifiers(seqnames(gr), sl)
    gr <- GenomicRanges::GRanges(seqnames = chromosome,
                                 ranges = IRanges::ranges(gr),
                                 strand = BiocGenerics::strand(gr),
                                 S4Vectors::mcols(gr))
    # shift position to transcript coordinates
    if(!is.null(tx)){
        tx <- .norm_tx(tx)
        .check_tx_sequences(tx, sequences)
        message("Shifting RMBase result's coordinates based on transcript ",
                "data ...")
        gr <- shiftGenomicToTranscript(gr, tx)
        f <- !duplicated(paste0(as.character(gr),"-",mcols(gr)$mod_type))
        gr <- gr[f]
    }
    # check if all the referenced modification match their originating base
    if(!is.null(sequences)){
        # remove any positions which do not match the originating base
        gr <- .remove_incompatible_mod_by_seq(gr, sequences)
        metadata <- .add_sequence_check_to_metadata(metadata)
    }
    #
    mcols(gr)$mod_id <- seq_along(gr)
    makeEpiTxDbFromGRanges(gr, metadata = metadata, reassign.ids = reassign.ids)
}


#' @rdname makeEpiTxDbFromRMBase
#' @importFrom curl curl
#' @export
listAvailableOrganismsFromRMBase <- function(){
    # con <- curl::curl(EPITXDB_RMBASE_URL)
    # page <- xml2::read_html(con)
    # organisms <- xml2::xml_attr(xml2::xml_find_all(page,'//img[@alt="[DIR]"]//../following::a'),"href")
    # organisms <- gsub("/","",organisms)
    # organisms[!(organisms %in% c("ajax","otherspecies"))]
    rmbase_data <- NULL
    utils::data("rmbase_data", envir = environment(), package = "EpiTxDb")
    as.character(unique(rmbase_data$organism))
}

#' @importFrom curl curl
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
    # files <- .get_RMBase_files(organism)
    # .get_RMBase_genomes(files)
    rmbase_data <- NULL
    utils::data("rmbase_data", envir = environment(), package = "EpiTxDb")
    as.character(unique(rmbase_data[rmbase_data$organism == organism,]$genome))
}

#' @rdname makeEpiTxDbFromRMBase
#' @export
listAvailableGenomesFromRMBase <- function(organism){
    if(!(organism %in% listAvailableOrganismsFromRMBase())){
        stop("'organism' must be a valid organism from ",
             "listAvailableOrganismsfromRMBase()")
    }
    .listAvailableGenomesFromRMBase(organism)
}

.get_RMBase_mod <- function(organism, genome){
    # f_genome <- vapply(strsplit(files,"_"),"[",character(1),2L) == genome
    # ans <- unique(vapply(strsplit(files[f_genome],"_"),"[",character(1),4L))
    # ans[!grepl("mod",ans)]
    rmbase_data <- NULL
    utils::data("rmbase_data", envir = environment(), package = "EpiTxDb")
    as.character(rmbase_data[rmbase_data$organism == organism &
                                 rmbase_data$genome == genome,]$mod)
}

#' @rdname makeEpiTxDbFromRMBase
#' @export
listAvailableModFromRMBase <- function(organism, genome){
    if(!(organism %in% listAvailableOrganismsFromRMBase())){
        stop("'organism' must be a valid organism from ",
             "listAvailableOrganismsfromRMBase()")
    }
    if(!(genome %in% listAvailableGenomesFromRMBase(organism))){
        stop("'genome' must be a valid genome for the given 'organism' from ",
             "listAvailableGenomesFromRMBase()")
    }
    .get_RMBase_mod(organism, genome)
}