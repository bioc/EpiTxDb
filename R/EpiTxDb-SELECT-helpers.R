# Low-level helpers imported from GenomicFeatures ------------------------------

.as_qualified <- GenomicFeatures:::.as_qualified
.tables_in_joins <- GenomicFeatures:::.tables_in_joins
.build_SQL_SELECT <- GenomicFeatures:::.build_SQL_SELECT
queryAnnotationDb <- GenomicFeatures:::queryAnnotationDb

# .EPITXDB_join_tables() -------------------------------------------------------

.EPITXDB_add_table_bundle <- function(tables){
  if (length(tables) == 1L)
    return(tables)
  skip_tables <- c("modification",  "seqnames")
  skip_tables <- tables %in% skip_tables
  if(all(skip_tables)){
    return(tables)
  }
  tables <- c(tables[skip_tables],
              pc(paste0("modification2",tables[!skip_tables]),
                 tables[!skip_tables]))
  unlist(tables)
}

.EPITXDB_join_tables <- function(tables, joinColumn = "_mod_id")
{
    # browser()
    tables <- unique(tables)
    if (length(tables) == 1L)
        return(tables)
    joinColumn <- unique(joinColumn)
    tables <- .EPITXDB_add_table_bundle(tables)
    joinColumn <- unique(joinColumn)
    stopifnot(length(tables) == length(joinColumn) + (length(tables) - 1L) %/% 2)
    ## Order tables & remove duplicates.
    join_order <- c("modification", "seqnames",
                    "modification2reaction", "reaction",
                    "modification2specifier", "specifier",
                    "modification2reference", "reference")
    joinColumn_order <- c("_sn_id","_mod_id","_rx_id","_spec_id","_ref_id")
    tables <- intersect(join_order, tables)
    joinColumn <- intersect(joinColumn_order, joinColumn)
    jc_offset_tx <- sum("_sn_id" %in% joinColumn)
    jc_offset_sel <- c("2"=1,"3"=2,"4"=3,"5"=3,"6"=4,"7"=4,"8"=5)
    joins <- character(2L * length(tables) - 1L)
    ON_idx <- 2L * seq_len(length(tables) - 1L)
    ON <- sapply(2:length(tables), function(i) {
        Rtable <- tables[[i]]
        if(i%%2 == jc_offset_tx){
            USING <- joinColumn[[1L+jc_offset_tx]]
            Ltable <- tables[[1L]]
        } else {
            USING <- joinColumn[[jc_offset_sel[[as.character(i)]]]]
            Ltable <- tables[[i-1L]]
        }
        Lcolumn <- .as_qualified(Ltable, USING)
        Rcolumn <- .as_qualified(Rtable, USING)
        paste(Lcolumn, Rcolumn, sep="=")
    })
    joins[ON_idx] <- ON
    joins[c(1L, ON_idx + 1L)] <- tables
    joins
}

EpiTxDb_SELECT_from_LEFT_JOIN <- function(epitxdb, table, columns,
                                          filter = list(),
                                          orderby = character(0),
                                          joinColumn = character(0))
{
    schema_version <- EpiTxDb_schema_version(epitxdb)
    tables <- EPITXDB_column2table(columns, from_table = table,
                                   schema_version = schema_version)
    if(!("modification" %in% tables)){
        tables <- c("modification",tables)
    }
    where_columns <- names(filter)
    where_tables <- EPITXDB_column2table(where_columns, from_table = table,
                                         schema_version = schema_version)
    joinColumn <- EPITXDB_table2joinColumns(tables = c(table, tables, where_tables),
                                            schema_version = schema_version)
    joins <- .EPITXDB_join_tables(c(table, tables, where_tables), joinColumn)
    orderby_tables <- EPITXDB_column2table(orderby, from_table = table,
                                           schema_version = schema_version)
    stopifnot(all(orderby_tables %in% .tables_in_joins(joins)))
    use_joins <- length(joins) > 1L
    if (use_joins) {
        columns <- .as_qualified(tables[names(tables) != ""], columns)
        names(filter) <- .as_qualified(where_tables, where_columns)
        orderby <- .as_qualified(orderby_tables, orderby)
    }
    ## .build_SQL_SELECT() uses LEFT joins.
    SQL <- .build_SQL_SELECT(columns, joins, distinct = use_joins,
                             filter = filter, orderby = orderby)
    queryAnnotationDb(epitxdb, SQL)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

EpiTxDb_SELECT_from_modification <- function(epitxdb, filter = list(),
                                             orderby = c("_mod_id", "_sn_id"))
{
    schema_version <- EpiTxDb_schema_version(epitxdb)
    columns <- EPITXDB_table_columns("modification",
                                     schema_version = schema_version)
    EpiTxDb_SELECT_from_LEFT_JOIN(epitxdb, "modification", columns,
                                  filter = filter, orderby = orderby)
}
EpiTxDb_SELECT_from_seqnames <- function(epitxdb, filter = list(),
                                         orderby = "_sn_id")
{
    schema_version <- EpiTxDb_schema_version(epitxdb)
    columns <- EPITXDB_table_columns("seqnames",
                                     schema_version = schema_version)
    EpiTxDb_SELECT_from_LEFT_JOIN(epitxdb, "seqnames", columns,
                                  filter = filter, orderby = orderby)
}

EpiTxDb_SELECT_from_reaction <- function(epitxdb, filter = list(),
                                         orderby = c("_rx_id", "rx_rank"))
{
    schema_version <- EpiTxDb_schema_version(epitxdb)
    columns <- EPITXDB_table_columns("reaction",
                                     schema_version = schema_version)
    EpiTxDb_SELECT_from_LEFT_JOIN(epitxdb, "reaction", columns,
                                  filter = filter, orderby = orderby)
}

EpiTxDb_SELECT_from_specifier <- function(epitxdb, filter = list(),
                                          orderby = c("_spec_id"))
{
    schema_version <- EpiTxDb_schema_version(epitxdb)
    columns <- EPITXDB_table_columns("specifier",
                                     schema_version = schema_version)
    EpiTxDb_SELECT_from_LEFT_JOIN(epitxdb, "specifier", columns,
                                  filter = filter, orderby = orderby)
}

EpiTxDb_SELECT_from_reference <- function(epitxdb, filter = list(),
                                          orderby = c("_ref_id"))
{
    schema_version <- EpiTxDb_schema_version(epitxdb)
    columns <- EPITXDB_table_columns("reference",
                                     schema_version = schema_version)
    EpiTxDb_SELECT_from_LEFT_JOIN(epitxdb, "reference", columns,
                                  filter = filter, orderby = orderby)
}
