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
#' title
#'
#' @return
#' @export
#'
#' @examples
.EpiTxDb <- setRefClass("EpiTxDb", contains = "AnnotationDb",
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
        mod_end = "integer",
        transcript_id = "character",
        transcript_name = "character",
        transcript_ensembltrans = "character"
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

load_reactions <- function(epitxdb, set.col.class = FALSE){
    reactions <- EpiTxDb_SELECT_from_reaction(epitxdb)
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

load_specifiers <- function(epitxdb, set.col.class = FALSE) {
  specifiers <- EpiTxDb_SELECT_from_specifier(epitxdb)
  colnames(specifiers) <- sub("^_", "", colnames(specifiers))
  .format_specifiers(specifiers, set.col.class = set.col.class)
}

.format_references <- function(references, set.col.class = FALSE){
  COL2CLASS <- c(
    mod_id="integer",
    reference_type="factor",
    reference="character"
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
      .valid.reaction.table(conn),
      .valid.specifier.table(conn),
      .valid.reference.table(conn))
}

setValidity("EpiTxDb", .valid.EpiTxDb)

# constructor ------------------------------------------------------------------

EpiTxDb <- function(conn) .EpiTxDb$new(conn =conn)

# methods ----------------------------------------------------------------------

#' @importFrom BiocGenerics organism
setMethod("organism", "EpiTxDb",
    function(object)
    {
        metadata <- metadata(object)
        metadata <- setNames(metadata[ , "value"],
                             tolower(metadata[ , "name"]))
        unname(metadata["organism"])
    }
)

# seqinfo ----------------------------------------------------------------------

.get_seqnames <- function(seqnames_db){
    name_not_empty <- seqnames_db$transcript_name != ""
    ensembl_not_empty <- seqnames_db$transcript_ensembltrans != ""
    both_empty <- !ensembl_not_empty & !name_not_empty
    if(all(name_not_empty)){
      return(seqnames_db$transcript_name)
    }
    if(all(ensembl_not_empty)){
      return(seqnames_db$transcript_ensembltrans)
    }
    if(all(both_empty)){
      return(seqnames_db$transcript_id)
    }
    seqnames <- seqnames_db$transcript_id
    if(length(intersect(seqnames,
                        seqnames_db$transcript_name)) != 0L){
        return(seqnames)
    }
    seqnames[name_not_empty] <-
      seqnames_db$transcript_name[name_not_empty]
    if(length(intersect(seqnames,
                        seqnames_db$transcript_ensembltrans)) != 0L){
        return(seqnames)
    }
    seqnames[ensembl_not_empty] <-
      seqnames_db$transcript_ensembltrans[ensembl_not_empty]
    seqnames
}

.get_EpiTxDb_seqinfo <- function(x){
    sql <- paste0("SELECT DISTINCT transcript_id, transcript_name, ",
                  "transcript_ensembltrans FROM modification")
    seqnames_db <- queryAnnotationDb(x, sql)
    ans <- Seqinfo(seqnames = .get_seqnames(seqnames_db))
    sql <- "SELECT value FROM metadata WHERE name='Genome'"
    genome <- unlist(queryAnnotationDb(x, sql))
    names(genome) <- NULL
    if (length(genome) != 0L)
      genome(ans) <- genome
    ans
}

setMethod("seqinfo", "EpiTxDb", .get_EpiTxDb_seqinfo)

# as.list and comparison -------------------------------------------------------

.format_txdb_dump <- function(modifications = NULL, reactions = NULL,
                              specifiers = NULL, references = NULL){
    modifications <- .format_modifications(modifications, set.col.class = TRUE)
    reactions <- .format_reactions(reactions, set.col.class = TRUE)
    specifiers <- .format_specifiers(specifiers, set.col.class = TRUE)
    references <- .format_references(references, set.col.class = TRUE)
    list(modifications = modifications, reactions = reactions,
         specifiers = specifiers, references = references)
}

setMethod("as.list", "EpiTxDb",
    function(x, ...)
    {
        modifications <- load_modifications(x)
        reactions <- load_reactions(x)
        specifiers <- load_specifiers(x)
        references <- load_references(x)
        .format_txdb_dump(modifications, reactions, specifiers, references)
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
