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
#' modifications(EpiTxDb.Hsapiens.hg38))
#' }
NULL

# helper functions for extracting feature data ---------------------------------

.extract_features <- function(epitxdb, table, mcolumns = character(0),
                              filter = list(), core_columns){
    schema_version <- EpiTxDb_schema_version(epitxdb)
    names(mcolumns) <- EPITXDB_column2table(mcolumns, from_table = table,
                                            schema_version = schema_version)
    proxy_column <- orderby <- c(modification = "_mod_id")

    ## 1st SELECT: extract stuff from the proxy table.
    columns1 <- union(core_columns, mcolumns[names(mcolumns) == table])
    df1 <- EpiTxDb_SELECT_from_INNER_JOIN(epitxdb, table, columns1,
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
        columns2 <- c(proxy_column, columns2)
        orderby <- c("_mod_id")
        EpiTxDb_SELECT_from_INNER_JOIN(epitxdb, satellite_table, columns2,
                                       filter = filter2, orderby = orderby)
    })
    df1 <- list(df1)
    names(df1) <- table
    c(df1, df_list)
}


# .extract_features_as_GRanges() -----------------------------------------------

.make_DataFrame_from_df_list <- GenomicFeatures:::.make_DataFrame_from_df_list
.assignMetadataList <- GenomicFeatures:::.assignMetadataList

.as_db_columns <- function(columns)
    sub("^(mod_id)$", "_\\1", columns)

.extract_features_as_GRanges <- function(epitxdb, table,
                                         mcolumns = character(0),
                                         filter = list(), use.names = FALSE) {
    if (!isTRUEorFALSE(use.names))
        stop("'use.names' must be TRUE or FALSE")
    db_mcolumns <- db_mcolumns0 <- .as_db_columns(mcolumns)
    proxy_columns <- EPITXDB_table_columns(table)
    if (use.names) {
        if ("name" %in% names(proxy_columns)) {
            proxy_name_column <- proxy_columns[["name"]]
            if (!(proxy_name_column %in% db_mcolumns0))
                db_mcolumns <- c(db_mcolumns0, proxy_name_column)
        } else {
            warning("no column in '", table, "' table to retrieve the feature ",
                    "names from")
            use.names <- FALSE
        }
    }
    names(filter) <- .as_db_columns(names(filter))
    core_columns <- proxy_columns[proxy_columns %in% EPITXDB_MOD_COLUMNS]
    df_list <- .extract_features(epitxdb, table, db_mcolumns,
                                 filter, core_columns)
    DF <- .make_DataFrame_from_df_list(df_list)
    DF$seqnames <- .get_seqnames(DF)
    ans <- GenomicRanges::makeGRangesFromDataFrame(
        DF,
        seqinfo = .get_EpiTxDb_seqinfo(epitxdb),
        seqnames.field = "seqnames",
        start.field = "mod_start",
        end.field = "mod_end")
    if (use.names) {
        ans_names <- DF[ , proxy_name_column]
        ans_names[is.na(ans_names)] <- ""  # replace NAs with empty strings
        names(ans) <- ans_names
    }
    mcols(ans) <- setNames(DF[db_mcolumns0], mcolumns)
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

.extractFromEpiTxDb <- function(epitxdb, table, mcolumns = NULL,
                                filter = NULL, use.names = FALSE){
    user_mcolumns <- mcolumns
    mcolumns <- translateCols(mcolumns, epitxdb)
    if (is.null(filter))
        filter <- list()
    names(filter) <- translateCols(names(filter), epitxdb)
    ans <- .extract_features_as_GRanges(epitxdb, table, mcolumns, filter,
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
        .extractFromEpiTxDb(x, "modification", mcolumns = column,
                            filter = filter, use.names = use.names)
    }
)
