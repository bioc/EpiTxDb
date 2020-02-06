# Low-level helpers imported from GenomicFeatures ------------------------------

.as_qualified <- GenomicFeatures:::.as_qualified
.tables_in_joins <- GenomicFeatures:::.tables_in_joins
.build_SQL_SELECT <- GenomicFeatures:::.build_SQL_SELECT
queryAnnotationDb <- GenomicFeatures:::queryAnnotationDb

# .EPITXDB_join_tables() -------------------------------------------------------

.EPITXDB_join_tables <- function(tables)
{
    tables <- unique(tables)
    if (length(tables) == 1L)
        return(tables)
    ## Order tables & remove duplicates.
    join_order <- c("modification", "reaction", "specifier", "transcript")
    tables <- intersect(join_order, tables)
    joins <- character(2L * length(tables) - 1L)
    ON_idx <- 2L * seq_len(length(tables) - 1L)
    ON <- sapply(2:length(tables), function(i) {
        Rtable <- tables[[i]]
        USING <- "_mod_id"
        Ltable <- tables[[1L]]
        Lcolumn <- .as_qualified(Ltable, USING)
        Rcolumn <- .as_qualified(Rtable, USING)
        paste(Lcolumn, Rcolumn, sep="=")
    })
    joins[ON_idx] <- ON
    joins[c(1L, ON_idx + 1L)] <- tables
    joins
}

EpiTxDb_SELECT_from_INNER_JOIN <- function(epitxdb, table, columns,
                                           filter = list(),
                                           orderby = character(0)){
    schema_version <- EpiTxDb_schema_version(epitxdb)
    tables <- EPITXDB_column2table(columns, from_table = table,
                                   schema_version = schema_version)
    where_columns <- names(filter)
    where_tables <- EPITXDB_column2table(where_columns, from_table = table,
                                         schema_version = schema_version)
    joins <- .EPITXDB_join_tables(c(table, tables, where_tables))
    orderby_tables <- EPITXDB_column2table(orderby, from_table = table,
                                           schema_version = schema_version)
    stopifnot(all(orderby_tables %in% .tables_in_joins(joins)))
    use_joins <- length(joins) > 1L
    if (use_joins) {
        columns <- .as_qualified(tables, columns)
        names(filter) <- .as_qualified(where_tables, where_columns)
        orderby <- .as_qualified(orderby_tables, orderby)
    }
    ## .build_SQL_SELECT() uses INNER joins.
    SQL <- .build_SQL_SELECT(columns, joins, distinct = use_joins,
                             filter = filter, orderby = orderby)
    queryAnnotationDb(epitxdb, SQL)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

EpiTxDb_SELECT_from_modification <- function(epitxdb, filter=list(),
                                             orderby="_mod_id")
{
    schema_version <- EpiTxDb_schema_version(epitxdb)
    columns <- EPITXDB_table_columns("modification",
                                     schema_version = schema_version)
    EpiTxDb_SELECT_from_INNER_JOIN(epitxdb, "modification", columns,
                                   filter = filter, orderby = orderby)
}

EpiTxDb_SELECT_from_reaction <- function(epitxdb, filter=list(),
                                         orderby=c("_mod_id", "mod_rank"))
{
    schema_version <- EpiTxDb_schema_version(epitxdb)
    columns <- EPITXDB_table_columns("reaction",
                                     schema_version = schema_version)
    EpiTxDb_SELECT_from_INNER_JOIN(epitxdb, "reaction", columns,
                                   filter = filter, orderby = orderby)
}

EpiTxDb_SELECT_from_specifier <- function(epitxdb, filter=list(),
                                          orderby=c("_mod_id"))
{
    schema_version <- EpiTxDb_schema_version(epitxdb)
    columns <- EPITXDB_table_columns("specifier",
                                     schema_version = schema_version)
    EpiTxDb_SELECT_from_INNER_JOIN(epitxdb, "specifier", columns,
                                   filter = filter, orderby = orderby)
}

EpiTxDb_SELECT_from_reference <- function(epitxdb, filter=list(),
                                          orderby=c("_mod_id"))
{
    schema_version <- EpiTxDb_schema_version(epitxdb)
    columns <- EPITXDB_table_columns("reference",
                                     schema_version = schema_version)
    EpiTxDb_SELECT_from_INNER_JOIN(epitxdb, "reference", columns,
                                   filter = filter, orderby = orderby)
}
