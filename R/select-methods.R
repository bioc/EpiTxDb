#' @include TxModDb-class.R
NULL

# helper functions for column name selection/conversion ------------------------

# Helpers to access/process the table names and columns
.getTableColMapping <- function(x){
    conn <- dbconn(x)
    tables <- DBI::dbListTables(conn)
    tCols <- sapply(tables, DBI::dbListFields, conn=conn)
    ## right up front we are getting rid of metadata tables...
    tCols[names(tCols)!="metadata"]
}

# used to match and generate abbreviated column names
.makeColAbbreviations <- function(x){
    tCols <- .getTableColMapping(x)
    longNames <- unique(unlist(tCols,use.names=FALSE))
    abbrev <- unique(toupper((gsub("_","",unlist(tCols,use.names=FALSE)))))
    abbrev <- gsub("^REACTION","RX",abbrev)
    abbrev <- gsub("^SPECIFIER","SPEC",abbrev)
    abbrev <- gsub("^TRANSCRIPT","TX",abbrev)
    names(abbrev) <- longNames
    abbrev
}

# For when you need to get the true table names back from the abbrev's
.reverseColAbbreviations <- function(x, cnames){
    abr <- .makeColAbbreviations(x)
    names(abr)[abr %in% cnames]
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
    prefSort <- c("m","r","s")
    res <- sTNames[match(prefSort, sTNames)]
    paste(res[!is.na(res)], collapse="")
}

.makeTableKey <- function(x,cnames){
    sTNames <- substr(.getSimpleTableNames(x, cnames),1,1)
    .encodeSortedTableKey(sTNames)
}

# real joins for likely combos
.tableJoinSelector <- function(tName){
    rs <- paste("(SELECT * FROM reaction ",
                "LEFT JOIN specifier ON (reaction._mod_id = specifier._mod_id) ) ")
    ms <- paste("SELECT * FROM modification ",
                "LEFT JOIN specifier ON (modification._mod_id = specifier._mod_id) ")
    mr <- paste("SELECT * FROM modification ",
                "LEFT JOIN reaction ON (modification._mod_id = reaction._mod_id) ")
    mrs <- paste("(SELECT * FROM modification ",
                 "LEFT JOIN reaction ON (modification._mod_id = reaction._mod_id) ",
                 "LEFT JOIN specifier ON (modification._mod_id = specifier._mod_id) ) ")

    sql <- switch(EXPR = tName,
                  "m" = "modification",
                  "r" = "reaction",
                  "s" = "specifier",
                  "rs" = rs,
                  "ms" = ms,
                  "mr" = mr,
                  "mrs" = mrs,
                  stop(paste("No query for this combination of tables.",
                             "Please add",tName,"to the interpolator")))
    sql
}

# helper function to generate SQL statements -----------------------------------

#
.makeJoinSQL <- function(x, cnames){
    tKey <- .makeTableKey(x,cnames)
    .tableJoinSelector(tKey)
}

.makeSelectList <- function(x, cnames, abbrev=TRUE){
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

.makeKeyList <- function(x, keys, keytype, abbrev=TRUE){
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
    ## Some argument checking
    if('skipValidKeysTest' %in% names(extraArgs)){
        skipValidKeysTest < -extraArgs[["skipValidKeysTest"]]
    } else {
        skipValidKeysTest <- FALSE
    }
    testSelectArgs(x, keys = keys, cols = columns, keytype = keytype,
                   skipValidKeysTest = skipValidKeysTest)

    ## 1st we check the keytype to see if it is valid:
    if(is.na(keys(x, keytype)[1]) & length(keys(x, keytype))==1){
        stop(paste("There do not appear to be any keys",
                   "for the keytype you have specified."))
    }

    # generate column names
    cnames <- unique(c(keytype, columns))

    ## the following just gets the major join and then modifies it ONLY if the
    ## keytype is a GENEID
    majorJoin <- .makeJoinSQL(x, cnames)
    # if(keytype=="GENEID"){
    #     majorJoin <- sub("FROM transcript LEFT JOIN gene",
    #                      "FROM transcript INNER JOIN gene",majorJoin)
    # }

    if(length(keys) <= 1000){  ##if more than about this, prolly faster to get all
        sql <- paste("SELECT DISTINCT",
                     .makeSelectList(x, cnames, abbrev=FALSE),
                     "FROM",
                     majorJoin,
                     "WHERE",
                     .makeKeyList(x, keys, keytype, abbrev=FALSE))
    } else {
        sql <- paste("SELECT DISTINCT",
                     .makeSelectList(x, cnames, abbrev=FALSE),
                     "FROM",
                     majorJoin)
    }

    res <- AnnotationDbi:::dbQuery(dbconn(x), sql)

    if(length(keys) > 1000){ ##Then drop the extras now(in event there are some)
        ktColId <- .reverseColAbbreviations(x, keytype)
        res <-  res[res[[ktColId]] %in% keys,,drop=FALSE]
    }

    ## Then drop any columns that were not explicitely requested but that may have
    ## been appended to make a joind (like TXID)
    res <- res[,.reverseColAbbreviations(x,cnames),drop=FALSE]

    ## Then sort rows and columns and drop the filtered rows etc. using resort_base
    ## from AnnotationDbi
    joinType <- .reverseColAbbreviations(x, keytype)
    if(dim(res)[1]>0){
        res <- resort_base(res, keys, joinType,
                           .reverseColAbbreviations(x,cnames))
    }

    ## Then put the user preferred headers onto the table
    fcNames <- .makeColAbbreviations(x)
    colnames(res) <- fcNames[match(colnames(res), names(fcNames))]
    res
}

setMethod("select", "TxModDb",
          function(x, keys, columns, keytype, ...) {
              .select(x, keys, columns, keytype, ...)
          }
)

# columns ----------------------------------------------------------------------

.columns <- function(x){
    res <- .makeColAbbreviations(x)
    names(res) <- NULL
    res
}


setMethod("columns", "TxModDb",
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
        "TXID" = AnnotationDbi:::dbQuery(
            dbconn(x),
            "SELECT DISTINCT transcript_id FROM modification", 1L),
        "TXENSEMBLTRANS" = AnnotationDbi:::dbQuery(
            dbconn(x),
            "SELECT DISTINCT transcript_ensembltrans FROM modification", 1L),
        "TXENTREZID" = AnnotationDbi:::dbQuery(
            dbconn(x),
            "SELECT DISTINCT transcript_entrezid FROM modification", 1L),
        "RXGENENAME" = AnnotationDbi:::dbQuery(
            dbconn(x),
            "SELECT DISTINCT reaction_genename FROM reaction", 1L),
        "RXENSEMBL" = AnnotationDbi:::dbQuery(
            dbconn(x),
            "SELECT DISTINCT reaction_ensembl FROM reaction", 1L),
        "RXENSEMBLTRANS" = AnnotationDbi:::dbQuery(
            dbconn(x),
            "SELECT DISTINCT reaction_ensembltrans FROM reaction", 1L),
        "RXENTREZID" = AnnotationDbi:::dbQuery(
            dbconn(x),
            "SELECT DISTINCT reaction_entrezid FROM reaction", 1L),
        "RXENZYME" = AnnotationDbi:::dbQuery(
            dbconn(x),
            "SELECT DISTINCT reaction_enzyme FROM reaction", 1L),
        "SPECTYPE" = AnnotationDbi:::dbQuery(
            dbconn(x),
            "SELECT DISTINCT specifier_type FROM specifier", 1L),
        "SPECGENENAME" = AnnotationDbi:::dbQuery(
            dbconn(x),
            "SELECT DISTINCT specifier_genename FROM specifier", 1L),
        "SPECENSEMBL" = AnnotationDbi:::dbQuery(
            dbconn(x),
            "SELECT DISTINCT specifier_entrezid FROM specifier", 1L),
        "SPECENTREZID" = AnnotationDbi:::dbQuery(
            dbconn(x),
            "SELECT DISTINCT specifier_ensembl FROM specifier", 1L),
        stop(paste(keytype, "is not a supported keytype.",
                   " Please use the keytypes",
                   "method to identify viable keytypes")))
    as.character(res)
}

.keysDispatch <- function(x, keytype, ...){
    if (missing(keytype)) keytype <- "MODID"
    AnnotationDbi:::smartKeys(x=x, keytype=keytype, ..., FUN=.keys)
}

setMethod("keys", "TxModDb",.keysDispatch)

# keytypes ---------------------------------------------------------------------

setMethod("keytypes", "TxModDb",
          function(x) return(c("MODID","MODTYPE","MODNAME","TXID",
                               "TXENSEMBLTRANS","TXENTREZID","RXGENENAME",
                               "RXENSEMBL","RXENSEMBLTRANS","RXENTREZID",
                               "RXENZYME","SPECTYPE","SPECGENENAME",
                               "SPECENSEMBL","SPECENTREZID"))
)
