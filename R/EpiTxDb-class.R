#' @include EpiTxDb.R
#' @include AllGenerics.R
#' @include EpiTxDb-schema.R
#' @include EpiTxDb-SELECT-helpers.R
NULL

#' @name EpiTxDb-class
#'
#' @title EpiTxDb objects
#'
#' @description
#' The \code{EpiTxDb} class is a
#' \code{\link[AnnotationDbi:AnnotationDb-class]{AnnotationDb}} type container
#' for storing Epitranscriptomic information.
#'
#' The information are typically stored on a per transcript and not as genomic
#' coordinates, but the \code{EpiTxDb} class is agnostic to this. In case of
#' genomic coordinates \code{transcriptsBy} will return modifications per
#' chromosome.
#'
#' @param x,object a \code{EpiTxDb} object
#'
#' @seealso
#'   \itemize{
#'     \item{\code{\link[=makeEpiTxDbFromGRanges]{makeEpiTxDbFromGRanges}} for
#'       creating a \code{EpiTxDb} object from a
#'       \code{\link[GenomicRanges:GRanges-class]{GRanges}} object and it's
#'        metadata columns}
#'     \item{\code{\link[=makeEpiTxDbFromRMBase]{makeEpiTxDbFromRMBase}} for
#'       creating a \code{EpiTxDb} object from RMBase online resources}
#'     \item{\code{\link[=makeEpiTxDbFromtRNAdb]{makeEpiTxDbFromtRNAdb}} for
#'       creating a \code{EpiTxDb} object from tRNAdb online resources}
#'     \item{\code{\link[=makeEpiTxDb]{makeEpiTxDb}} for creating a
#'       \code{EpiTxDb} object from \code{data.frames}}
#'     \item{\code{\link[=modifications]{modifications}},
#'       \code{\link[=modifications]{modificationsBy}} for getting
#'       epitranscriptomic modification locations}
#'     \item{\code{\link[=select]{select}} for using the default interface of
#'       \code{\link[AnnotationDbi:AnnotationDb-class]{AnnotationDb}} objects.}
#'     \item{\code{\link[=shiftGenomicToTranscript]{shiftGenomicToTranscript}}
#'       and \code{\link[=shiftGenomicToTranscript]{shiftTranscriptToGenomic}}
#'       for transfering genomic to transcript coordinates and back again.}
#' }
#'
#' @export
#'
#' @examples
#' etdb_file <- system.file("extdata", "EpiTxDb.Hs.hg38.snoRNAdb.sqlite",
#'                         package="EpiTxDb")
#' etdb <- loadDb(etdb_file)
#' etdb
#'
#' # general methods
#' seqinfo(etdb) #
#' seqlevels(etdb) # easy access to all transcript names
.EpiTxDb <- setRefClass("EpiTxDb", contains = "AnnotationDb",
                        fields = list( ),
                        methods = list(
                            initialize = function(...) {
                                callSuper(...)
                                .self
                            }))


# Low-level data loaders -------------------------------------------------------

makeZeroRowDataFrame <- GenomicFeatures:::makeZeroRowDataFrame
setDataFrameColClass <- GenomicFeatures:::setDataFrameColClass

.format_modifications <- function(modifications, set.col.class = FALSE){
    COL2CLASS <- c(
        mod_id = "integer",
        mod_type = "character",
        mod_name = "character",
        mod_start = "integer",
        mod_end = "integer",
        tx_id = "integer"
    )
    if (is.null(modifications)) {
        modifications <- makeZeroRowDataFrame(COL2CLASS)
    } else {
        if (!is.data.frame(modifications))
            stop("'modifications' must be a data frame")
        if (!identical(names(modifications), names(COL2CLASS)))
            modifications <- modifications[names(COL2CLASS)]
        if (set.col.class)
            modifications <- setDataFrameColClass(modifications, COL2CLASS)
    }
    modifications
}

load_modifications <- function(epitxdb, set.col.class = FALSE){
    modifications <- EpiTxDb_SELECT_from_modification(epitxdb)
    colnames(modifications) <- sub("^_", "", colnames(modifications))
    .format_modifications(modifications, set.col.class = set.col.class)
}

.format_transcripts <- function(transcripts, set.col.class = FALSE){
    COL2CLASS <- c(
        tx_id = "integer",
        tx_name = "character",
        tx_ensembl = "character"
    )
    if (is.null(transcripts)) {
      transcripts <- makeZeroRowDataFrame(COL2CLASS)
    } else {
        if (!is.data.frame(transcripts))
            stop("'transcripts' must be a data frame")
        if (!identical(names(transcripts), names(COL2CLASS)))
          transcripts <- transcripts[names(COL2CLASS)]
        if (set.col.class)
            transcripts <- setDataFrameColClass(transcripts, COL2CLASS)
    }
    transcripts
}

load_transcripts <- function(epitxdb, set.col.class = FALSE){
    transcripts <- EpiTxDb_SELECT_from_transcript(epitxdb)
    colnames(transcripts) <- sub("^_", "", colnames(transcripts))
    .format_transcripts(transcripts, set.col.class = set.col.class)
}

.format_reactions <- function(reactions, set.col.class = FALSE){
    COL2CLASS <- c(
        rx_id = "integer",
        rx_genename = "character",
        rx_rank = "integer",
        rx_ensembl = "character",
        rx_ensembltrans = "character",
        rx_entrezid = "character"
    )
    if (is.null(reactions)) {
        reactions <- makeZeroRowDataFrame(COL2CLASS)
    } else {
        if (!is.data.frame(reactions))
            stop("'reactions' must be a data frame")
        if (!identical(names(reactions), names(COL2CLASS)))
            reactions <- reactions[names(COL2CLASS)]
        if (set.col.class)
            reactions <- setDataFrameColClass(reactions, COL2CLASS)
    }
    reactions
}

load_reactions <- function(epitxdb, set.col.class = FALSE){
    reactions <- EpiTxDb_SELECT_from_reaction(epitxdb)
    colnames(reactions) <- sub("^_", "", colnames(reactions))
    .format_reactions(reactions,  set.col.class = set.col.class)
}

.format_specifiers <- function(specifiers, set.col.class = FALSE){
    COL2CLASS <- c(
        spec_id="integer",
        spec_type="factor",
        spec_genename="character",
        spec_ensembl="character",
        spec_ensembltrans="character",
        spec_entrezid="character"
    )
    if (is.null(specifiers)) {
        specifiers <- makeZeroRowDataFrame(COL2CLASS)
    } else {
        if (!is.data.frame(specifiers))
            stop("'specifiers' must be a data frame")
        if (!identical(names(specifiers), names(COL2CLASS)))
            specifiers <- specifiers[names(COL2CLASS)]
        if (set.col.class)
            specifiers <- setDataFrameColClass(specifiers, COL2CLASS)
    }
    specifiers
}

load_specifiers <- function(epitxdb, set.col.class = FALSE) {
    specifiers <- EpiTxDb_SELECT_from_specifier(epitxdb)
    colnames(specifiers) <- sub("^_", "", colnames(specifiers))
    .format_specifiers(specifiers, set.col.class = set.col.class)
}

.format_references <- function(references, set.col.class = FALSE){
    COL2CLASS <- c(
        ref_id="integer",
        ref_type="factor",
        ref="character"
    )
    if (is.null(references)) {
       references <- makeZeroRowDataFrame(COL2CLASS)
    } else {
        if (!is.data.frame(references))
            stop("'references' must be a data frame")
        if (!identical(names(references), names(COL2CLASS)))
            references <- references[names(COL2CLASS)]
        if (set.col.class)
            references <- setDataFrameColClass(references, COL2CLASS)
    }
    references
}

load_references <- function(epitxdb, set.col.class = FALSE) {
    references <- EpiTxDb_SELECT_from_reference(epitxdb)
    colnames(references) <- sub("^_", "", colnames(references))
    .format_references(references, set.col.class = set.col.class)
}


# validity EpiTxDb -------------------------------------------------------------

.valid.modification.table <- function(conn){
    schema_version <- EpiTxDb_schema_version(conn)
    colnames <- EPITXDB_table_columns("modification",
                                      schema_version=schema_version)
    msg <- AnnotationDbi:::.valid.table.colnames(conn, "modification", colnames)
    if (!is.null(msg))
      return(msg)
    NULL
}

.valid.transcript.table <- function(conn){
    schema_version <- EpiTxDb_schema_version(conn)
    colnames <- EPITXDB_table_columns("transcript",
                                      schema_version=schema_version)
    msg <- AnnotationDbi:::.valid.table.colnames(conn, "transcript", colnames)
    if (!is.null(msg))
      return(msg)
    NULL
}

.valid.reaction.table <- function(conn){
    schema_version <- EpiTxDb_schema_version(conn)
    colnames <- EPITXDB_table_columns("reaction",
                                      schema_version=schema_version)
    msg <- AnnotationDbi:::.valid.table.colnames(conn, "reaction", colnames)
    if (!is.null(msg))
      return(msg)
    NULL
}

.valid.specifier.table <- function(conn){
    schema_version <- EpiTxDb_schema_version(conn)
    colnames <- EPITXDB_table_columns("specifier",
                                      schema_version=schema_version)
    msg <- AnnotationDbi:::.valid.table.colnames(conn, "specifier", colnames)
    if (!is.null(msg))
      return(msg)
    NULL
}

.valid.reference.table <- function(conn){
    schema_version <- EpiTxDb_schema_version(conn)
    colnames <- EPITXDB_table_columns("reference",
                                      schema_version=schema_version)
    msg <- AnnotationDbi:::.valid.table.colnames(conn, "reference", colnames)
    if (!is.null(msg))
      return(msg)
    NULL
}

.valid.EpiTxDb <- function(object){
    conn <- dbconn(object)
    c(AnnotationDbi:::.valid.metadata.table(conn, DB_TYPE_NAME, DB_TYPE_VALUE),
        .valid.modification.table(conn),
        .valid.transcript.table(conn),
        .valid.reaction.table(conn),
        .valid.specifier.table(conn),
        .valid.reference.table(conn))
}

setValidity("EpiTxDb", .valid.EpiTxDb)

# constructor ------------------------------------------------------------------

EpiTxDb <- function(conn) .EpiTxDb$new(conn = conn)

# methods ----------------------------------------------------------------------

#' @rdname EpiTxDb-class
#' @importFrom BiocGenerics organism
#' @export
setMethod("organism", "EpiTxDb",
    function(object)
    {
        metadata <- metadata(object)
        names <- tolower(metadata[ , "name"])
        metadata <- metadata[ , "value"]
        names(metadata) <- names
        unname(metadata["organism"])
    }
)

# seqinfo ----------------------------------------------------------------------

.get_seqnames <- function(seqnames_db){
    if(!all(lengths(seqnames_db$`_tx_id`) == 1L) ||
       !all(lengths(seqnames_db$tx_name) == 1L) ||
       !all(lengths(seqnames_db$tx_ensembl) == 1L)){
        stop(".")
    }
    tx_name <- seqnames_db$tx_name
    tx_ensembl <- seqnames_db$tx_ensembl
    tx_id <- seqnames_db$`_tx_id`
    if(is(tx_name,"List")){
        tx_name <- unlist(tx_name)
    }
    if(is(tx_ensembl,"List")){
        tx_ensembl <- unlist(tx_ensembl)
    }
    if(is(tx_id,"List")){
        tx_id <- unlist(tx_id)
    }
    name_not_empty <- tx_name != ""
    ensembl_not_empty <- tx_ensembl != ""
    both_empty <- !ensembl_not_empty & !name_not_empty
    if(all(name_not_empty)){
        return(tx_name)
    }
    if(all(ensembl_not_empty)){
        return(tx_ensembl)
    }
    if(all(both_empty)){
        return(tx_id)
    }
    seqnames <- tx_id
    if(length(intersect(seqnames, tx_name)) != 0L){
        return(seqnames)
    }
    seqnames[name_not_empty] <- tx_name[name_not_empty]
    if(length(intersect(seqnames, tx_ensembl)) != 0L){
        return(seqnames)
    }
    seqnames[ensembl_not_empty] <- tx_ensembl[ensembl_not_empty]
    seqnames
}

.get_EpiTxDb_seqinfo <- function(x){
    sql <- paste0("SELECT DISTINCT _tx_id, tx_name, tx_ensembl FROM ",
                  "transcript")
    seqnames_db <- queryAnnotationDb(x, sql)
    ans <- Seqinfo(seqnames = .get_seqnames(seqnames_db))
    sql <- "SELECT value FROM metadata WHERE name='Genome'"
    genome <- unlist(queryAnnotationDb(x, sql))
    names(genome) <- NULL
    if (length(genome) != 0L)
      genome(ans) <- genome
    ans
}

#' @rdname EpiTxDb-class
#' @importFrom GenomeInfoDb seqinfo
#' @export
setMethod("seqinfo", "EpiTxDb", .get_EpiTxDb_seqinfo)

# as.list and comparison -------------------------------------------------------

.format_txdb_dump <- function(modifications = NULL, transcripts = NULL,
                              reactions = NULL, specifiers = NULL,
                              references = NULL){
    modifications <- .format_modifications(modifications, set.col.class = TRUE)
    transcripts <- .format_transcripts(transcripts, set.col.class = TRUE)
    reactions <- .format_reactions(reactions, set.col.class = TRUE)
    specifiers <- .format_specifiers(specifiers, set.col.class = TRUE)
    references <- .format_references(references, set.col.class = TRUE)
    list(modifications = modifications, transcripts = transcripts,
         reactions = reactions, specifiers = specifiers,
         references = references)
}

#' @rdname EpiTxDb-class
#' @export
setMethod("as.list", "EpiTxDb",
    function(x)
    {
        modifications <- load_modifications(x)
        transcripts <- load_transcripts(x)
        reactions <- load_reactions(x)
        specifiers <- load_specifiers(x)
        references <- load_references(x)
        .format_txdb_dump(modifications, transcripts, reactions, specifiers,
                          references)
    }
)

compareEpiTxDbs <- function(epitxdb1, epitxdb2)
{
    if (!is(epitxdb1, "EpiTxDb") || !is(epitxdb2, "EpiTxDb"))
        stop("'epitxdb1' and 'epitxdb2' must be EpiTxDb objects")
    epitxdb1_dump <- as.list(epitxdb1)
    epitxdb2_dump <- as.list(epitxdb2)
    identical(epitxdb1_dump, epitxdb2_dump)
}

# coercion ---------------------------------------------------------------------


# show -------------------------------------------------------------------------

setMethod("show", "EpiTxDb",
          function(object)
          {
              cat(class(object), "object:\n")
              metadata <- metadata(object)
              for (i in seq_len(nrow(metadata))) {
                  cat("# ", metadata[i, "name"], ": ", metadata[i, "value"],
                      "\n", sep ="")
              }
          }
)
