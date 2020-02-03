#' @include TxModDb.R
#' @include AllGenerics.R
#' @include TxModDb-schema.R
#' @include TxModDb-SELECT-helpers.R
NULL

#' @name TxModDb-class
#'
#' @title TxModDb objects
#'
#' @description
#' title
#'
#' @return
#' @export
#'
#' @examples
.TxModDb <- setRefClass("TxModDb", contains = "AnnotationDb",
                        fields = list( ),
                        methods = list(
                            initialize = function(...) {
                                callSuper(...)
                                if (length(dbListTables(conn) != 0L)) {


                                }
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
        mod_end = "integer"
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

load_modifications <- function(txdb, set.col.class = FALSE){
    modifications <- TxModDb_SELECT_from_modification(txdb)
    colnames(modifications) <- sub("^_", "", colnames(modifications))
    .format_modifications(modifications, set.col.class = set.col.class)
}

.format_reactions <- function(reactions, set.col.class = FALSE){
    COL2CLASS <- c(
        mod_id = "integer",
        mod_rank = "factor",
        reaction_genename = "factor",
        reaction_ensembl = "factor",
        reaction_ensembltrans = "factor",
        reaction_entrezid = "factor",
        reaction_enzyme = "factor"
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

load_reactions <- function(txdb, set.col.class = FALSE){
    reactions <- TxModDb_SELECT_from_reaction(txdb)
    colnames(reactions) <- sub("^_", "", colnames(reactions))
    .format_reactions(reactions,  set.col.class = set.col.class)
}

.format_specifiers <- function(specifiers, set.col.class = FALSE){
    COL2CLASS <- c(
        mod_id="integer",
        specifier_type="factor",
        specifier_genename="character",
        specifier_entrezid="character",
        specifier_ensembl="character"
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

load_specifiers <- function(txdb, set.col.class = FALSE) {
    specifiers <- TxModDb_SELECT_from_specifier(txdb)
    colnames(specifiers) <- sub("^_", "", colnames(specifiers))
    .format_specifiers(specifiers, set.col.class = set.col.class)
}

.format_transcripts <- function(transcripts, set.col.class = FALSE){
    COL2CLASS <- c(
        mod_id="integer",
        entrezid="character"
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

load_transcripts <- function(txdb, set.col.class = FALSE)
{
  transcripts <- TxModDb_SELECT_from_transcript(txdb)
  colnames(transcripts) <- sub("^_", "", colnames(transcripts))
  .format_transcripts(transcripts, set.col.class = set.col.class)
}

# validity TxModDb -------------------------------------------------------------

.valid.modification.table <- function(conn){
    schema_version <- TxModDb_schema_version(conn)
    colnames <- TXMODDB_table_columns("modification",
                                      schema_version=schema_version)
    msg <- AnnotationDbi:::.valid.table.colnames(conn, "modification", colnames)
    if (!is.null(msg))
      return(msg)
    NULL
}

.valid.reaction.table <- function(conn){
    schema_version <- TxModDb_schema_version(conn)
    colnames <- TXMODDB_table_columns("reaction",
                                      schema_version=schema_version)
    msg <- AnnotationDbi:::.valid.table.colnames(conn, "reaction", colnames)
    if (!is.null(msg))
      return(msg)
    NULL
}

.valid.specifier.table <- function(conn){
    schema_version <- TxModDb_schema_version(conn)
    colnames <- TXMODDB_table_columns("specifier",
                                      schema_version=schema_version)
    msg <- AnnotationDbi:::.valid.table.colnames(conn, "specifier", colnames)
    if (!is.null(msg))
      return(msg)
    NULL
}

.valid.transcript.table <- function(conn){
    schema_version <- TxModDb_schema_version(conn)
    colnames <- TXMODDB_table_columns("transcript",
                                      schema_version=schema_version)
    msg <- AnnotationDbi:::.valid.table.colnames(conn, "transcript", colnames)
    if (!is.null(msg))
      return(msg)
    NULL
}

.valid.TxModDb <- function(object){
    conn <- dbconn(object)
    c(AnnotationDbi:::.valid.metadata.table(conn, DB_TYPE_NAME, DB_TYPE_VALUE),
      .valid.modification.table(conn),
      .valid.reaction.table(conn),
      .valid.specifier.table(conn),
      .valid.transcript.table(conn))
}

setValidity("TxModDb", .valid.TxModDb)

# constructor ------------------------------------------------------------------

TxModDb <- function(conn) .TxModDb$new(conn =conn)

# methods ----------------------------------------------------------------------

#' @importFrom BiocGenerics organism
setMethod("organism", "TxModDb",
    function(object)
    {
        metadata <- metadata(object)
        metadata <- setNames(metadata[ , "value"],
                             tolower(metadata[ , "name"]))
        unname(metadata["organism"])
    }
)

# as.list and comparison -------------------------------------------------------

.format_txdb_dump <- function(modifications = NULL, reactions = NULL,
                              specifiers = NULL, transcripts = NULL){
    modifications <- .format_modifications(modifications, set.col.class = TRUE)
    reactions <- .format_reactions(reactions, set.col.class = TRUE)
    specifiers <- .format_specifiers(specifiers, set.col.class = TRUE)
    transcripts <- .format_transcripts(transcripts, set.col.class = TRUE)
    list(modifications = modifications, reactions = reactions,
         specifiers = specifiers, transcripts = transcripts)
}

setMethod("as.list", "TxModDb",
    function(x, ...)
    {
        modifications <- load_modifications(x)
        reactions <- load_reactions(x)
        specifiers <- load_specifiers(x)
        transcripts <- load_transcripts(x)
        .format_txdb_dump(modifications, reactions, specifiers, transcripts)
    }
)

compareTxModDbs <- function(txmoddb1, txmoddb2)
{
    if (!is(txmoddb1, "TxModDb") || !is(txmoddb2, "TxModDb"))
      stop("'txmoddb1' and 'txmoddb2' must be TxModDb objects")
    txmoddb1_dump <- as.list(txmoddb1)
    txmoddb2_dump <- as.list(txmoddb2)
    identical(txmoddb1_dump, txmoddb2_dump)
}

# coercion ---------------------------------------------------------------------


# show -------------------------------------------------------------------------

setMethod("show", "TxModDb",
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
