# Low-level helpers imported from GenomicFeatures ------------------------------

.as_qualified <- GenomicFeatures:::.as_qualified
.tables_in_joins <- GenomicFeatures:::.tables_in_joins
.build_SQL_SELECT <- GenomicFeatures:::.build_SQL_SELECT
queryAnnotationDb <- GenomicFeatures:::queryAnnotationDb

# .TXMODDB_join_tables() -------------------------------------------------------

.TXMODDB_join_tables <- function(tables)
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

TxModDb_SELECT_from_INNER_JOIN <- function(txmoddb, table, columns,
                                           filter = list(),
                                           orderby = character(0)){
    schema_version <- TxModDb_schema_version(txmoddb)
    tables <- TXMODDB_column2table(columns, from_table = table,
                                   schema_version = schema_version)
    where_columns <- names(filter)
    where_tables <- TXMODDB_column2table(where_columns, from_table = table,
                                         schema_version = schema_version)
    joins <- .TXMODDB_join_tables(c(table, tables, where_tables))
    orderby_tables <- TXMODDB_column2table(orderby, from_table = table,
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
    queryAnnotationDb(txmoddb, SQL)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

TxModDb_SELECT_from_modification <- function(txmoddb, filter=list(),
                                             orderby="_mod_id")
{
    schema_version <- TxModDb_schema_version(txmoddb)
    columns <- TXMODDB_table_columns("modification",
                                     schema_version = schema_version)
    TxModDb_SELECT_from_INNER_JOIN(txmoddb, "modification", columns,
                                   filter = filter, orderby = orderby)
}

TxModDb_SELECT_from_reaction <- function(txmoddb, filter=list(),
                                         orderby=c("_mod_id", "mod_rank"))
{
    schema_version <- TxModDb_schema_version(txmoddb)
    columns <- TXMODDB_table_columns("reaction",
                                     schema_version = schema_version)
    TxModDb_SELECT_from_INNER_JOIN(txmoddb, "reaction", columns,
                                   filter = filter, orderby = orderby)
}

TxModDb_SELECT_from_specifier <- function(txmoddb, filter=list(),
                                          orderby=c("_mod_id"))
{
    schema_version <- TxModDb_schema_version(txmoddb)
    columns <- TXMODDB_table_columns("specifier",
                                     schema_version = schema_version)
    TxModDb_SELECT_from_INNER_JOIN(txmoddb, "specifier", columns,
                                   filter = filter, orderby = orderby)
}

TxModDb_SELECT_from_transcript <- function(txmoddb, filter=list(),
                                           orderby=c("_mod_id", "entrezid"))
{
    schema_version <- TxModDb_schema_version(txmoddb)
    columns <- TXMODDB_table_columns("transcript",
                                     schema_version = schema_version)
    TxModDb_SELECT_from_INNER_JOIN(txmoddb, "transcript", columns,
                                   filter = filter, orderby = orderby)
}

