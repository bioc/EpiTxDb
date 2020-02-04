# =========================================================================
# TxModDb schema
# -------------------------------------------------------------------------
#
# Nothing in this file is exported.
#
# 3 tables:
#   - modification
#   - reaction
#   - specifier
#   - transcript
#   - metadata (not described here)

DB_TYPE_NAME <- "Db type"
DB_TYPE_VALUE <- "TxModDb"
DB_SCHEMA_VERSION <- "1.0"

# Return the *effective* schema version.
TxModDb_schema_version <- function(txmoddb)
{
    conn <- if (is(txmoddb, "TxModDb")) dbconn(txmoddb) else txmoddb
    version <- AnnotationDbi:::.getMetaValue(conn, "DBSCHEMAVERSION")
    numeric_version(version)
}


# Table columns ----------------------------------------------------------------

# 'modification' table

TXMODDB_MOD_COLDEFS <- c(
    `_mod_id` = "INTEGER PRIMARY KEY",
    mod_type = "TEXT NOT NULL",
    mod_name = "TEXT NULL",
    mod_start = "INTEGER NOT NULL",
    mod_end = "INTEGER NOT NULL",
    transcript_id = "TEXT NOT NULL",
    transcript_ensembltrans = "TEXT NULL",
    transcript_entrezid = "TEXT NULL"
)
TXMODDB_MOD_COLUMNS <- names(TXMODDB_MOD_COLDEFS)

# modification 'reaction' table

TXMODDB_REACTION_COLDEFS <- c(
    `_mod_id` = "INTEGER PRIMARY KEY",
    mod_rank = "INTEGER NOT NULL",
    reaction_genename = "TEXT NULL",
    reaction_ensembl = "TEXT NULL",
    reaction_ensembltrans = "TEXT NULL",
    reaction_entrezid = "TEXT NULL",
    reaction_enzyme = "TEXT NOT NULL"
)

TXMODDB_REACTION_COLUMNS <- names(TXMODDB_REACTION_COLDEFS)

# modification 'specifier' table

TXMODDB_SPECIFIER_COLDEFS <- c(
    `_mod_id` = "INTEGER NOT NULL",
    specifier_type = "TEXT NOT NULL",
    specifier_genename = "TEXT NOT NULL",
    specifier_entrezid = "TEXT NULL",
    specifier_ensembl = "TEXT NULL"
)

TXMODDB_SPECIFIER_COLUMNS <- names(TXMODDB_SPECIFIER_COLDEFS)

#

TXMODDB_COLUMNS <- list(
    modification = TXMODDB_MOD_COLUMNS,
    reaction = TXMODDB_REACTION_COLUMNS,
    specifier = TXMODDB_SPECIFIER_COLUMNS
)


# Build CREATE TABLE statements ------------------------------------------------

.build_SQL_CREATE_TABLE <- function(table, coldefs, constraints = NULL)
{
    SQL <- "CREATE TABLE %s (%s\n)"
    coldefs <- c(paste(names(coldefs), coldefs), constraints)
    coldefs <- paste("\n  ", coldefs, collapse = ",")
    sprintf(SQL, table, coldefs)
}

build_SQL_CREATE_modification_table <- function()
{
    unique_key <- "UNIQUE (_mod_id)"
    constraints <- c(unique_key)
    .build_SQL_CREATE_TABLE("modification", TXMODDB_MOD_COLDEFS, constraints)
}

build_SQL_CREATE_reaction_table <- function()
{
    unique_key <- "UNIQUE (_mod_id, mod_rank)"
    foreign_key <- "FOREIGN KEY (_mod_id) REFERENCES modification"
    constraints <- c(unique_key, foreign_key)
    .build_SQL_CREATE_TABLE("reaction", TXMODDB_REACTION_COLDEFS, constraints)
}

build_SQL_CREATE_specifier_table <- function()
{
    foreign_key <- "FOREIGN KEY (_mod_id) REFERENCES modification"
    constraints <- c(foreign_key)
    .build_SQL_CREATE_TABLE("specifier", TXMODDB_SPECIFIER_COLDEFS, constraints)
}


# helper functions -------------------------------------------------------------

TXMODDB_tables <- function() names(TXMODDB_COLUMNS)

TXMODDB_table_columns <- function(table, schema_version = NA){
    columns <- TXMODDB_COLUMNS[[table]]
    columns
}

TXMODDB_column2table <- function(columns, from_table = NA, schema_version = NA){
    if (length(columns) == 0L)
        return(character(0))
    tables <- sapply(columns,
                     function(column) {
                         for (table in TXMODDB_tables()) {
                             table_columns <-
                                 TXMODDB_table_columns(table,
                                                       schema_version = schema_version)
                             if (column %in% table_columns)
                                 return(table)
                         }
                         if (is.na(schema_version)) {
                             in_schema <- ""
                         } else {
                             in_schema <- c(" in db schema ", as.character(schema_version))
                         }
                         stop(column, ": no such column", in_schema)
                     })
    if (!is.na(from_table)) {
        table_columns <- TXMODDB_table_columns(from_table,
                                               schema_version = schema_version)
        tables[columns %in% table_columns] <- from_table
    }
    tables
}
