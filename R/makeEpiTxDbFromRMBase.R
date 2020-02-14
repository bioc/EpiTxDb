#' @include EpiTxDb-class.R
#' @include makeEpiTxDb.R
#' @include shiftGenomicToTranscript.R
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

#' @rdname makeEpiTxDbfromRMBase
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

#' @rdname makeEpiTxDbfromRMBase
#' @export
makeEpiTxDbfromRMBase <- function(organism, genome, type, tx, sequences = NULL,
                                  metadata = NULL, reassign.ids = FALSE){
    message("Loading RMBase files ...")
    files <- downloadRMBaseFiles(organism, genome, type)
    makeEpiTxDbfromRMBaseFiles(files, tx = tx, sequences = sequences,
                               metadata = metadata, reassign.ids = reassign.ids)
}

# makeEpiTxDbfromRMBaseFiles ---------------------------------------------------

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
.simplify_chromosome_identifiers <- function(chromosome, seqlevels){
    chromosome <- as.character(chromosome)
    # fix some weird mito chromosome annotation
    f <- chromosome %in% "Mito"
    if(any(f)){
        chromosome[f] <- gsub("Mito","M",chromosome[f])
    }
    # mismatching chromosomes
    mm_chr <- !(chromosome %in% seqlevels)
    f_chr_remove <- grepl("Chr|chr",chromosome[mm_chr])
    if(any(f_chr_remove)){
        if(all(levels(factor(gsub("Chr|chr","",chromosome[mm_chr][f_chromosome]))) %in%
               levels(chromosome[mm_chr][f_chromosome])) ||
           all(levels(chromosome[mm_chr][f_chromosome]) %in%
               levels(factor(gsub("Chr|chr","",chromosome[mm_chr][f_chromosome]))))){
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

.extract_GRanges_from_RMBase <- function(rmb, seqlevels = NA, seqtype = "RNA"){
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
    base_mm <- Modstrings:::originatingBase(codec)[f] != subseq(seq,21L,21L)
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
    chromosome <- .simplify_chromosome_identifiers(chromosome, seqlevels)
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
    mcols <- DataFrame(mod_id = mod_id,
                       mod_type = mod_type,
                       mod_name = mod_name,
                       chromosome = chromosome,
                       gene_name = gene_name,
                       score = rmb$score,
                       reference_type = ref_type,
                       reference = reference)
    # create GRanges result
    ans <- GenomicRanges::GRanges(seqnames = chromosome,
                                  ranges = IRanges::IRanges(pos, width = 1L),
                                  strand = strand,
                                  mcols)
    # expand by gene names
    gene_name <-
        IRanges::CharacterList(strsplit(as.character(mcols(ans)$gene_name),","))
    ans <- ans[unlist(Map(rep,seq_along(gene_name),lengths(gene_name)))]
    mcols(ans)$gene_name <- unlist(gene_name)
    ans <- ans[mcols(ans)$gene_name != ""]
    ans <- ans[!duplicated(ans) | !duplicated(mcols(ans))]
    #
    ans
}

.get_RMBase_data <- function(files, seqlevels = NA){
    grl <- lapply(files,
                  function(file){
                      rmb <- .read_RMBase_file(file)
                      .extract_GRanges_from_RMBase(rmb, seqlevels)
                  })
    grl
}

.gene_names_to_chromosome <- function(gr){
    chromosome <- mcols(gr)$gene_name
    GenomicRanges::GRanges(seqnames = chromosome,
                           ranges = ranges(gr),
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

#' @rdname makeEpiTxDbfromRMBase
#' @export
getRMBaseDataAsGRanges <- function(files, tx = NULL, sequences = NULL,
                                   shift.to.transcript = FALSE,
                                   check.vs.sequence = FALSE){
    message("Assembling data ...")
    # getting raw data from RMBase files
    sl <- NA
    if(!is.null(tx)){
        sl <- GenomeInfoDb::seqlevels(tx)
    }
    grl <- .get_RMBase_data(files, seqlevels = sl)
    gr <- unlist(GenomicRanges::GRangesList(grl))
    # shoft position to transcript coordinates
    if(shift.to.transcript){
        tx <- .norm_tx(tx)
        .check_tx_sequences(tx, sequences)
        message("Shifting RMBase result's coordinates based on transcript data ...")
        gr <- shiftGenomicToTranscript(gr, tx)
        f <- !duplicated(paste0(as.character(gr),"-",mcols(gr)$mod_type))
        gr <- gr[f]
        # check if all the referenced modification match their originating base
        if(!is.null(sequences) & check.vs.sequence){
            # remove any positions which do not match the originating base
            colnames(mcols(gr)) <- gsub("mod_type","mod",colnames(mcols(gr)))
            gr <- Modstrings::removeIncompatibleModifications(gr, sequences)
            colnames(mcols(gr)) <- gsub("^mod$","mod_type",colnames(mcols(gr)))
        }
    }
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

#' @rdname makeEpiTxDbfromRMBase
#' @export
makeEpiTxDbfromRMBaseFiles <- function(files, tx, sequences = NULL,
                                       metadata = NULL, reassign.ids = FALSE){
    gr <- getRMBaseDataAsGRanges(files, tx = tx, sequences = sequences,
                                 shift.to.transcript = TRUE,
                                 check.vs.sequence = TRUE)
    #
    mcols(gr)$mod_id <- seq_along(gr)
    mcols(gr)$transcript_name <- mcols(gr)$tx_id
    mcols(gr)$transcript_id <-
        as.integer(factor(mcols(gr)$transcript_name,
                          unique(mcols(gr)$transcript_name)))
    mcols(gr)$tx_id <- NULL
    if(!is.null(sequences)){
        metadata <- .add_sequence_check_to_metadata(metadata)
    }
    makeEpiTxDbfromGRanges(gr, metadata = metadata, reassign.ids = reassign.ids)
}


#' @rdname makeEpiTxDbfromRMBase
#' @export
listAvailableOrganismsFromRMBase <- function(){
    # con <- curl::curl(EPITXDB_RMBASE_URL)
    # page <- xml2::read_html(con)
    # organisms <- xml2::xml_attr(xml2::xml_find_all(page,'//img[@alt="[DIR]"]//../following::a'),"href")
    # organisms <- gsub("/","",organisms)
    # organisms[!(organisms %in% c("ajax","otherspecies"))]
    load(system.file("extdata", "rmbase_data.rda", package = "EpiTxDb"))
    as.character(unique(rmbase_data$organism))
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
    # files <- .get_RMBase_files(organism)
    # .get_RMBase_genomes(files)
    load(system.file("extdata", "rmbase_data.rda", package = "EpiTxDb"))
    as.character(unique(rmbase_data[rmbase_data$organism == organism,]$genome))
}

#' @rdname makeEpiTxDbfromRMBase
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
    load(system.file("extdata", "rmbase_data.rda", package = "EpiTxDb"))
    as.character(rmbase_data[rmbase_data$organism == organism &
                                 rmbase_data$genome == genome,]$mod)
}

#' @rdname makeEpiTxDbfromRMBase
#' @export
listAvailableModFromRMBase <- function(organism, genome){
    if(!(organism %in% listAvailableOrganismsFromRMBase())){
        stop("'organism' must be a valid organism from ",
             "listAvailableOrganismsfromRMBase()")
    }
    if(!(genome %in% listAvailableGenomesFromRMBase(organism))){
        stop("'genome' must be a valid genome for the fiven 'organism' from ",
             "listAvailableGenomesFromRMBase()")
    }
    .get_RMBase_mod(organism, genome)
}