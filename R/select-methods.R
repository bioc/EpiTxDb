#' @include EpiTxDb-class.R
NULL

#' @name select
#' @aliases select columns keys keytypes
#'
#' @title Using the "select" interface on \code{EpiTxDb} objects
#'
#' @description
#'   As expected for a \code{AnnotationDb} object, the general accessors
#'   \code{select}, \code{keys}, \code{columns} and \code{keytypes} can be used
#'   to get information from a \code{\link[=EpiTxDb-class]{EpiTxDb}} object.
#'
#' @param x a \code{\link[=EpiTxDb-class]{EpiTxDb}} object
#' @param keys,columns,keytype,... See
#'   \code{\link[AnnotationDbi:AnnotationDb-class]{AnnotationDb}} for more
#'   comprehensive description of the methods \code{select}, \code{keys},
#'   \code{columns} and \code{keytypes} and their arguments.
#'
#' @return a \code{data.frame} object for \code{select()} and a \code{character}
#'   vecor for \code{keytypes()}, \code{keys()} and \code{columns()}.
#'
#' @examples
#' etdb_file <- system.file("extdata", "EpiTxDb.Hs.hg38.snoRNAdb.sqlite",
#'                          package="EpiTxDb")
#' etdb <- loadDb(etdb_file)
#' etdb
NULL

# helper functions for column name selection/conversion ------------------------

# Helpers to access/process the table names and columns
#' @importFrom DBI dbListTables dbListFields
.getTableColMapping <- function(x){
    conn <- dbconn(x)
    tables <- DBI::dbListTables(conn)
    tCols <- sapply(tables, DBI::dbListFields, conn = conn)
    ## right up front we are getting rid of metadata tables...
    tCols[names(tCols) != "metadata"]
}

# used to match and generate abbreviated column names
.makeColAbbreviations <- function(x){
    tCols <- .getTableColMapping(x)
    longNames <- unique(unlist(tCols, use.names = FALSE))
    abbrev <- unique(toupper((gsub("_","",unlist(tCols, use.names = FALSE)))))
    abbrev <- gsub("^SEQNAMES", "SN", abbrev)
    abbrev <- gsub("^REACTION", "RX", abbrev)
    abbrev <- gsub("^SPECIFIER", "SPEC", abbrev)
    abbrev <- gsub("^REFERENCE", "REF", abbrev)
    names(abbrev) <- longNames
    abbrev
}

# For when you need to get the true table names back from the abbrev's
.reverseColAbbreviations <- function(x, cnames){
    abr <- .makeColAbbreviations(x)
    names(abr)[match(cnames,abr)]
}

# used to retrieve vector of table names that go with vector of col names
.getTableNames <- function(x, cnames){
    realColNames <- .reverseColAbbreviations(x, cnames)
    tCols <- .getTableColMapping(x)
    ## Translate all names to one unique vector.
    getTabNames <- function(name, tCols){names(tCols[grep(name, tCols)])}
    tabNames <- lapply(realColNames, getTabNames, tCols)
    names(tabNames) <- realColNames
    tabNames
}

.getSimpleTableNames <- function(x, cnames){
    unique(unlist(.getTableNames(x, cnames)))
}

# helper functions for generating SQL statements from selected columns ---------

# this just takes the 1 letter abrevs and makes them into a sorted string
# that can be used as a key below
.encodeSortedTableKey <- function(sTNames){
    prefSort <- c("mod","seq","rea","spe","ref")
    res <- sTNames[match(prefSort, sTNames)]
    paste(res[!is.na(res)], collapse="")
}

.makeTableKey <- function(x,cnames){
    sTNames <- substr(.getSimpleTableNames(x, cnames),1,3)
    .encodeSortedTableKey(sTNames)
}

# real joins for likely combos
.tableJoinSelector <- function(tName){
    #
    modseq <- paste("(SELECT * FROM modification ",
                     "JOIN seqnames ON (modification._sn_id = seqnames._sn_id) ) ")

    modrea <- paste("(SELECT * FROM modification ",
                    "LEFT JOIN modification2reaction ON (modification._mod_id = modification2reaction._mod_id) ",
                    "LEFT JOIN reaction ON (modification2reaction._rx_id = reaction._rx_id) )")

    modspe <- paste("(SELECT * FROM modification ",
                    "LEFT JOIN modification2specifier ON (modification._mod_id = modification2specifier._mod_id) ",
                    "LEFT JOIN specifier ON (modification2specifier._spec_id = specifier._spec_id) )")

    modref <- paste("(SELECT * FROM modification ",
                    "LEFT JOIN modification2reference ON (modification._mod_id = modification2reference._mod_id) ",
                    "LEFT JOIN reference ON (modification2reference._ref_id = reference._ref_id) )")
    #
    modseqrea <- paste("(SELECT * FROM modification ",
                    "JOIN seqnames ON (modification._sn_id = seqnames._sn_id) ",
                    "LEFT JOIN modification2reaction ON (modification._mod_id = modification2reaction._mod_id) ",
                    "LEFT JOIN reaction ON (modification2reaction._rx_id = reaction._rx_id) )")

    modseqspe <- paste("(SELECT * FROM modification ",
                    "JOIN seqnames ON (modification._sn_id = seqnames._sn_id) ",
                    "LEFT JOIN modification2specifier ON (modification._mod_id = modification2specifier._mod_id) ",
                    "LEFT JOIN specifier ON (modification2specifier._spec_id = specifier._spec_id) )")

    modseqref <- paste("(SELECT * FROM modification ",
                    "JOIN seqnames ON (modification._sn_id = seqnames._sn_id) ",
                    "LEFT JOIN modification2reference ON (modification._mod_id = modification2reference._mod_id) ",
                    "LEFT JOIN reference ON (modification2reference._ref_id = reference._ref_id) )")
    #
    modreaspe <- paste("(SELECT * FROM modification ",
                       "LEFT JOIN modification2reaction ON (modification._mod_id = modification2reaction._mod_id) ",
                       "LEFT JOIN reaction ON (modification2reaction._rx_id = reaction._rx_id)  ",
                       "LEFT JOIN modification2specifier ON (modification._mod_id = modification2specifier._mod_id) ",
                       "LEFT JOIN specifier ON (modification2specifier._spec_id = specifier._spec_id) ) ")

    modrearef <- paste("(SELECT * FROM modification ",
                       "LEFT JOIN modification2reaction ON (modification._mod_id = modification2reaction._mod_id) ",
                       "LEFT JOIN reaction ON (modification2reaction._rx_id = reaction._rx_id)  ",
                       "LEFT JOIN modification2reference ON (modification._mod_id = modification2reference._mod_id) ",
                       "LEFT JOIN reference ON (modification2reference._ref_id = reference._ref_id) ) ")
    #
    modsperef <- paste("(SELECT * FROM modification ",
                       "LEFT JOIN modification2specifier ON (modification._mod_id = modification2specifier._mod_id) ",
                       "LEFT JOIN specifier ON (modification2specifier._spec_id = specifier._spec_id) ",
                       "LEFT JOIN modification2reference ON (modification._mod_id = modification2reference._mod_id) ",
                       "LEFT JOIN reference ON (modification2reference._ref_id = reference._ref_id) ) ")
    #
    modreasperef <- paste("(SELECT * FROM modification ",
                          "JOIN seqnames ON (modification._sn_id = seqnames._sn_id) ",
                          "LEFT JOIN modification2reaction ON (modification._mod_id = modification2reaction._mod_id) ",
                          "LEFT JOIN reaction ON (modification2reaction._rx_id = reaction._rx_id) ",
                          "LEFT JOIN modification2specifier ON (modification._mod_id = modification2specifier._mod_id) ",
                          "LEFT JOIN specifier ON (modification2specifier._spec_id = specifier._spec_id) ",
                          "LEFT JOIN modification2reference ON (modification._mod_id = modification2reference._mod_id) ",
                          "LEFT JOIN reference ON (modification2reference._ref_id = reference._ref_id) ) ")
    #
    modseqreaspe <- paste("(SELECT * FROM modification ",
                          "JOIN seqnames ON (modification._sn_id = seqnames._sn_id) ",
                          "LEFT JOIN modification2reaction ON (modification._mod_id = modification2reaction._mod_id) ",
                          "LEFT JOIN reaction ON (modification2reaction._rx_id = reaction._rx_id) ",
                          "LEFT JOIN modification2specifier ON (modification._mod_id = modification2specifier._mod_id) ",
                          "LEFT JOIN specifier ON (modification2specifier._spec_id = specifier._spec_id) ) ")

    modseqrearef <- paste("(SELECT * FROM modification ",
                          "JOIN seqnames ON (modification._sn_id = seqnames._sn_id) ",
                          "LEFT JOIN modification2reaction ON (modification._mod_id = modification2reaction._mod_id) ",
                          "LEFT JOIN reaction ON (modification2reaction._rx_id = reaction._rx_id) ",
                          "LEFT JOIN modification2reference ON (modification._mod_id = modification2reference._mod_id) ",
                          "LEFT JOIN reference ON (modification2reference._ref_id = reference._ref_id) ) ")

    modseqsperef <- paste("(SELECT * FROM modification ",
                          "JOIN seqnames ON (modification._sn_id = seqnames._sn_id) ",
                          "LEFT JOIN modification2specifier ON (modification._mod_id = modification2specifier._mod_id) ",
                          "LEFT JOIN specifier ON (modification2specifier._spec_id = specifier._spec_id) ",
                          "LEFT JOIN modification2reference ON (modification._mod_id = modification2reference._mod_id) ",
                          "LEFT JOIN reference ON (modification2reference._ref_id = reference._ref_id) ) ")

    modreasperef <- paste("(SELECT * FROM modification ",
                          "LEFT JOIN modification2reaction ON (modification._mod_id = modification2reaction._mod_id) ",
                          "LEFT JOIN reaction ON (modification2reaction._rx_id = reaction._rx_id) ",
                          "LEFT JOIN modification2specifier ON (modification._mod_id = modification2specifier._mod_id) ",
                          "LEFT JOIN specifier ON (modification2specifier._spec_id = specifier._spec_id) ",
                          "LEFT JOIN modification2reference ON (modification._mod_id = modification2reference._mod_id) ",
                          "LEFT JOIN reference ON (modification2reference._ref_id = reference._ref_id) ) ")
    #
    modseqreasperef <- paste("(SELECT * FROM modification ",
                             "JOIN seqnames ON (modification._sn_id = seqnames._sn_id) ",
                             "LEFT JOIN modification2reaction ON (modification._mod_id = modification2reaction._mod_id) ",
                             "LEFT JOIN reaction ON (modification2reaction._rx_id = reaction._rx_id) ",
                             "LEFT JOIN modification2specifier ON (modification._mod_id = modification2specifier._mod_id) ",
                             "LEFT JOIN specifier ON (modification2specifier._spec_id = specifier._spec_id) ",
                             "LEFT JOIN modification2reference ON (modification._mod_id = modification2reference._mod_id) ",
                             "LEFT JOIN reference ON (modification2reference._ref_id = reference._ref_id) ) ")

    sql <- switch(EXPR = tName,
                  "mod" = "modification",
                  "seq" = "seqnames",
                  "rea" = "reaction",
                  "spe" = "specifier",
                  "ref" = "reference",
                  "modseq" = modseq,
                  "modrea" = modrea,
                  "modspe" = modspe,
                  "modref" = modref,
                  "seqrea" = modseqrea,
                  "seqspe" = modseqspe,
                  "seqref" = modseqref,
                  "reaspe" = modreaspe,
                  "rearef" = modrearef,
                  "speref" = modsperef,
                  "modseqrea" = modseqrea,
                  "modseqspe" = modseqspe,
                  "modseqref" = modseqref,
                  "modreaspe" = modreaspe,
                  "modrearef" = modrearef,
                  "modsperef" = modsperef,
                  "seqreaspe" = modseqreaspe,
                  "seqrearef" = modseqrearef,
                  "seqsperef" = modseqsperef,
                  "reasperef" = modreasperef,
                  "modseqreaspe" = modseqreaspe,
                  "modseqrearef" = modseqrearef,
                  "modreasperef" = modreasperef,
                  "modseqreasperef" = modseqreasperef,
                  stop(paste("No query for this combination of tables.",
                             "Please add ",tName," to the interpolator")))
    sql
}

# helper function to generate SQL statements -----------------------------------

#
.makeJoinSQL <- function(x, cnames){
    tKey <- .makeTableKey(x,cnames)
    .tableJoinSelector(tKey)
}

.makeSelectList <- function(x, cnames, abbrev = TRUE){
    tNames <- .getTableNames(x, cnames)
    ## Here is where we only grab the 1st one...
    tNames <- lapply(tNames,function(x){x[1]})
    ## then continue on...
    tabAbbrevs <- substr(unlist(tNames),1,1)
    names(tabAbbrevs) <- rep(names(tNames),elementNROWS(tNames))
    if(abbrev==TRUE){
        paste(paste0(tabAbbrevs, ".", names(tabAbbrevs)), collapse=", ")
    } else {
        paste(names(tabAbbrevs), collapse=", ")
    }
}

.makeKeyList <- function(x, keys, keytype, abbrev = TRUE){
    colType <- .makeSelectList(x, keytype, abbrev)
    keys <- paste(paste0("'", keys, "'"), collapse=",")
    paste(colType, "IN (", keys,")")
}


# select -----------------------------------------------------------------------

.select <- function(x, keys, columns, keytype, ...){
    extraArgs <- list(...)
    if(missing(keys) || !is.character(keys))
        stop("'keys' must be a character vector")
    if(missing(columns) || !is.character(columns))
        stop("'columns' must be a character vector")
    # Some argument checking
    if('skipValidKeysTest' %in% names(extraArgs)){
        skipValidKeysTest < -extraArgs[["skipValidKeysTest"]]
    } else {
        skipValidKeysTest <- FALSE
    }
    testSelectArgs(x, keys = keys, cols = columns, keytype = keytype,
                   skipValidKeysTest = skipValidKeysTest)

    # 1st we check the keytype to see if it is valid:
    if(is.na(keys(x, keytype)[1]) & length(keys(x, keytype)) == 1){
        stop(paste("There do not appear to be any keys",
                   "for the keytype you have specified."))
    }

    # generate column names
    cnames <- unique(c(keytype, columns))
    # create thejoin SQL
    majorJoin <- .makeJoinSQL(x, cnames)

    if(length(keys) <= 1000){  ##if more than about this, prolly faster to get all
        sql <- paste("SELECT DISTINCT",
                     .makeSelectList(x, cnames, abbrev = FALSE),
                     "FROM",
                     majorJoin,
                     "WHERE",
                     .makeKeyList(x, keys, keytype, abbrev = FALSE))
    } else {
        sql <- paste("SELECT DISTINCT",
                     .makeSelectList(x, cnames, abbrev = FALSE),
                     "FROM",
                     majorJoin)
    }

    res <- AnnotationDbi:::dbQuery(dbconn(x), sql)

    if(length(keys) > 1000){ ##Then drop the extras now(in event there are some)
        ktColId <- .reverseColAbbreviations(x, keytype)
        res <-  res[res[[ktColId]] %in% keys,,drop = FALSE]
    }

    ## Then drop any columns that were not explicitely requested but that may have
    ## been appended to make a joind (like TXID)
    res <- res[,.reverseColAbbreviations(x,cnames),drop = FALSE]

    ## Then sort rows and columns and drop the filtered rows etc. using resort_base
    ## from AnnotationDbi
    joinType <- .reverseColAbbreviations(x, keytype)
    if(dim(res)[1]>0){
        res <- resort_base(res, keys, joinType,
                           .reverseColAbbreviations(x, cnames))
    }

    ## Then put the user preferred headers onto the table
    fcNames <- .makeColAbbreviations(x)
    colnames(res) <- fcNames[match(colnames(res), names(fcNames))]
    res
}

#' @rdname select
#' @export
setMethod("select", "EpiTxDb",
          function(x, keys, columns, keytype, ...) {
              .select(x, keys, columns, keytype, ...)
          }
)

# columns ----------------------------------------------------------------------

.columns <- function(x){
    res <- .makeColAbbreviations(x)
    unname(res)
}

#' @rdname select
#' @export
setMethod("columns", "EpiTxDb",
          function(x) .columns(x)
)

# keys -------------------------------------------------------------------------

.keys <- function(x, keytype){
    testForValidKeytype(x, keytype)
    res <- switch(
        EXPR = keytype,
        "MODID" = AnnotationDbi:::dbQuery(
            dbconn(x),
            "SELECT DISTINCT _mod_id FROM modification", 1L),
        "MODTYPE" = AnnotationDbi:::dbQuery(
            dbconn(x),
            "SELECT DISTINCT mod_type FROM modification", 1L),
        "MODNAME" = AnnotationDbi:::dbQuery(
            dbconn(x),
            "SELECT DISTINCT mod_name FROM modification", 1L),
        "MODSTRAND" = AnnotationDbi:::dbQuery(
            dbconn(x),
            "SELECT DISTINCT mod_strand FROM modification", 1L),
        "SNID" = AnnotationDbi:::dbQuery(
            dbconn(x),
            "SELECT DISTINCT _sn_id FROM seqnames", 1L),
        "SNNAME" = AnnotationDbi:::dbQuery(
            dbconn(x),
            "SELECT DISTINCT sn_name FROM seqnames", 1L),
        "RXGENENAME" = AnnotationDbi:::dbQuery(
            dbconn(x),
            "SELECT DISTINCT rx_genename FROM reaction", 1L),
        "RXENSEMBL" = AnnotationDbi:::dbQuery(
            dbconn(x),
            "SELECT DISTINCT rx_ensembl FROM reaction", 1L),
        "RXENSEMBLTRANS" = AnnotationDbi:::dbQuery(
            dbconn(x),
            "SELECT DISTINCT rx_ensembltrans FROM reaction", 1L),
        "RXENTREZID" = AnnotationDbi:::dbQuery(
            dbconn(x),
            "SELECT DISTINCT rx_entrezid FROM reaction", 1L),
        "SPECTYPE" = AnnotationDbi:::dbQuery(
            dbconn(x),
            "SELECT DISTINCT spec_type FROM specifier", 1L),
        "SPECGENENAME" = AnnotationDbi:::dbQuery(
            dbconn(x),
            "SELECT DISTINCT spec_genename FROM specifier", 1L),
        "SPECENSEMBL" = AnnotationDbi:::dbQuery(
            dbconn(x),
            "SELECT DISTINCT spec_entrezid FROM specifier", 1L),
        "SPECENTREZID" = AnnotationDbi:::dbQuery(
            dbconn(x),
            "SELECT DISTINCT spec_ensembl FROM specifier", 1L),
        "REFTYPE" = AnnotationDbi:::dbQuery(
            dbconn(x),
            "SELECT DISTINCT ref_type FROM reference", 1L),
        "REF" = AnnotationDbi:::dbQuery(
            dbconn(x),
            "SELECT DISTINCT ref FROM reference", 1L),
        stop(paste(keytype, "is not a supported keytype.",
                   " Please use the keytypes",
                   "method to identify viable keytypes")))
    as.character(res)
}

.keysDispatch <- function(x, keytype, ...){
    if (missing(keytype)) keytype <- "MODID"
    AnnotationDbi:::smartKeys(x = x, keytype = keytype, ..., FUN = .keys)
}

#' @rdname select
#' @export
setMethod("keys", "EpiTxDb", .keysDispatch)

# keytypes ---------------------------------------------------------------------

#' @rdname select
#' @export
setMethod("keytypes", "EpiTxDb",
          function(x) return(c("MODID","MODTYPE","MODNAME","MODSTRAND","SNID",
                               "SNNAME","RXGENENAME","RXENSEMBL",
                               "RXENSEMBLTRANS","RXENTREZID", "SPECTYPE",
                               "SPECGENENAME","SPECENSEMBL","SPECENSEMBLTRANS",
                               "SPECENTREZID","REFTYPE","REF"))
)
