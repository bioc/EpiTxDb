#' @include EpiTxDb-class.R
NULL

#' @name modifications
#' @aliases modifications modificationsBy
#'
#' @title modifications
#'
#' @description
#' title
#'
#' @param x a \code{\link[=EpiTxDb-class]{EpiTxDb}}
#' @param columns XXX
#' @param filter XXX
#' @use.names XXX
#'
#' @examples
#' \dontrun{
#' library(EpiTxDb.Hsapiens.hg38)
#' modifications(EpiTxDb.Hsapiens.hg38.snoRNAdb())
#' }
NULL

# helper functions for extracting feature data ---------------------------------

.extract_modifications <- function(epitxdb, table, mcolumns = character(0),
                                   filter = list(), core_columns){
    schema_version <- EpiTxDb_schema_version(epitxdb)
    names(mcolumns) <- EPITXDB_column2table(mcolumns, from_table = table,
                                            schema_version = schema_version)
    proxy_column <- orderby <- c(modification = "_mod_id")

    ## 1st SELECT: extract stuff from the proxy table.
    columns1 <- union(core_columns, mcolumns[names(mcolumns) == table])
    df1 <- EpiTxDb_SELECT_from_LEFT_JOIN(epitxdb, table, columns1,
                                         filter = filter, orderby = orderby)

    ## Additional SELECTs: 1 additional SELECT per satellite table
    foreign_columns <- mcolumns[names(mcolumns) != table]
    foreign_columns <- split(foreign_columns, names(foreign_columns))
    satellite_tables <- names(foreign_columns)
    names(satellite_tables) <- satellite_tables
    df_list <- lapply(satellite_tables, function(satellite_table) {
        columns2 <- foreign_columns[[satellite_table]]
        if (length(filter) == 0L) {
            filter2 <- list()
        } else {
            filter2 <- list(df1[[proxy_column]])
            names(filter2) <- proxy_column
        }
        if (satellite_table %in% c("transcript","reaction","specifier",
                                   "reference")) {
            columns2 <- c(proxy_column, columns2)
            EpiTxDb_SELECT_from_LEFT_JOIN(epitxdb, satellite_table, columns2,
                                          filter = filter2, orderby = orderby)
        } else {
            stop(satellite_table, ": unsupported satellite table")
        }
    })
    df1 <- list(df1)
    names(df1) <- table
    ans <- c(df1, df_list)
    ans
}


# .extract_features_as_GRanges() -----------------------------------------------

.make_DataFrame_from_df_list <- GenomicFeatures:::.make_DataFrame_from_df_list
.assignMetadataList <- GenomicFeatures:::.assignMetadataList

.as_db_columns <- function(columns)
    sub("^(mod_id|tx_id|rx_id|spec_id|ref_id)$", "_\\1", columns)

.extract_modifications_as_GRanges <- function(epitxdb, mcolumns = character(0),
                                              filter = list(),
                                              use.names = FALSE)
{
    if (!isTRUEorFALSE(use.names)) {
        stop("'use.names' must be TRUE or FALSE")
    }
    table <- c("modification")
    # save the selected metadata columns for subsetting later on
    mcolumns0 <- mcolumns
    # add the transcript columns, since the always need to be extracted
    mcolumns <- unique(c(mcolumns,"tx_id","tx_name","tx_ensembl"))
    db_mcolumns <- db_mcolumns0 <- .as_db_columns(mcolumns)
    core_columns <- EPITXDB_table_columns(table)
    names(filter) <- .as_db_columns(names(filter))
    df_list <- .extract_modifications(epitxdb, table, db_mcolumns, filter,
                                      core_columns)
    DF <- .make_DataFrame_from_df_list(df_list)
    DF$seqnames <- .get_seqnames(DF)
    ans <- GenomicRanges::makeGRangesFromDataFrame(
        DF,
        seqinfo = .get_EpiTxDb_seqinfo(epitxdb),
        seqnames.field = "seqnames",
        start.field = "mod_start",
        end.field = "mod_end")
    if (use.names) {
        ans_names <- DF[ ,"mod_name"]
        names(ans) <- ans_names
    }
    mcols(ans) <- DF[db_mcolumns0]
    colnames(mcols(ans)) <- mcolumns
    # subset to the selected columns
    mcols(ans) <- mcols(ans)[,colnames(mcols(ans)) %in% mcolumns0,drop = FALSE]
    ans
}

# translate external to internal column names
translateCols <- function(columns, epitxdb){
    ## preserve any names
    oriColNames <- names(columns)
    ## and save the original column strings
    oriCols <- columns

    oriLen <- length(columns) ## not always the same as length(oriCols)
    ## get the available abbreviations as a translation vector (exp)
    names <- .makeColAbbreviations(epitxdb)
    exp <- sub("^_","", names(names))
    names(exp) <- names

    ## Then replace only those IDs that match the UC names
    m <- match(oriCols, names(exp))
    idx = which(!is.na(m))
    columns[idx] <- exp[m[idx]]

    if(length(columns) == oriLen && is.null(oriColNames)) {
        names(columns) <- oriCols
    } else if(length(columns) == oriLen && !is.null(oriColNames)) {
        names(columns) <- oriColNames
    } else if(length(columns) != oriLen) {
        stop("names were lost in translateCols() helper")
    }
    columns
}

.extractFromEpiTxDb <- function(epitxdb, mcolumns = NULL,
                                filter = NULL, use.names = FALSE){
    user_mcolumns <- mcolumns
    mcolumns <- translateCols(mcolumns, epitxdb)
    if (is.null(filter))
        filter <- list()
    names(filter) <- translateCols(names(filter), epitxdb)
    ans <- .extract_modifications_as_GRanges(epitxdb, mcolumns, filter,
                                             use.names)
    names(mcols(ans)) <- if (is.null(names(user_mcolumns))) user_mcolumns
    else names(user_mcolumns)
    .assignMetadataList(ans, epitxdb)
}

#' @rdname modifications
#' @export
setMethod("modifications", "EpiTxDb",
    function(x, column = c("mod_id","mod_type","mod_name"),
             filter = NULL, use.names = FALSE){
        .extractFromEpiTxDb(x, mcolumns = column, filter = filter,
                            use.names = use.names)
    }
)
