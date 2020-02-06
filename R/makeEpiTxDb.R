#' @include EpiTxDb-class.R
#' @include utils.R
NULL

#' @name makeEpiTxDb
#'
#' @title makeEpiTxDb
#'
#' @description
#' title
#'
#'
NULL

# check helper functions -------------------------------------------------------

.is_character_or_factor <- GenomicFeatures:::.is_character_or_factor
.check_foreign_key <- GenomicFeatures:::.check_foreign_key
translateIds <- GenomicFeatures:::translateIds
check_colnames <- GenomicFeatures:::check_colnames
has_col <- GenomicFeatures:::has_col
dbEasyQuery <- GenomicFeatures:::dbEasyQuery

.makeEpiTxDb_normarg_modifications <- function(modifications){
    .REQUIRED_COLS <- c("mod_id", "mod_type", "mod_start", "mod_end",
                        "transcript_id")
    .OPTIONAL_COLS <- c("mod_name", "transcript_name", "ensembltrans")
    # make sure 'transcript_id is set'
    if(!has_col(modifications, "transcript_id") &&
       !has_col(modifications, "ensembltrans") &&
       !has_col(modifications, "transcript_name")){
        stop("'modifications' must contain a column 'ensembltrans' or ",
             "'transcript_name', if no column 'transcript_id' is given.")
    }
    if(!has_col(modifications, "transcript_id") &&
       (has_col(modifications, "ensembltrans") ||
        has_col(modifications, "transcript_name"))){
        transcript_id <- factor(paste0(modifications$ensembltrans,"_",
                                       modifications$transcript_name))
        modifications$transcript_id <- as.integer(transcript_id)
    }
    check_colnames(modifications, .REQUIRED_COLS, .OPTIONAL_COLS,
                   "modifications")
    ## Check 'mod_id'.
    if (!is.integer(modifications$mod_id) || any(is.na(modifications$mod_id)))
        stop("'modifications$mod_id' must be an integer vector ",
                  "with no NAs")
    if (any(duplicated(modifications$mod_id)))
        stop("'modifications$mod_id' contains duplicated values")
    ## Check 'mod_start'.
    if (!is.numeric(modifications$mod_start)
        || any(is.na(modifications$mod_start)))
        stop("'modifications$mod_start' must be an integer vector ",
                  "with no NAs")
    if (!is.integer(modifications$mod_start))
        modifications$mod_start <- as.integer(modifications$mod_start)
    ## Check 'mod_end'.
    if (!is.numeric(modifications$mod_end)
        || any(is.na(modifications$mod_end)))
        stop("'transcripts$mod_end' must be an integer vector ",
                  "with no NAs")
    if (!is.integer(modifications$mod_end))
        modifications$mod_end <- as.integer(modifications$mod_end)
    ## Check 'mod_start <= mod_end'.
    if (any(modifications$mod_start > modifications$mod_end))
        stop("modification starts must be <= modification ends")
    ## Check 'mod_name'.
    if (!.is_character_or_factor(modifications$mod_name))
        stop("'modifications$mod_name' must be a character vector ",
                  "(or factor)")
    ## Check 'mod_type'.
    if (!.is_character_or_factor(modifications$mod_type))
        stop("'modifications$mod_type' must be a character vector ",
                  "(or factor)")
    if(!all(as.character(modifications$mod_type) %in% shortName(ModRNAString()))){
        stop("'modifications$mod_type' must be a valid modification shortName ",
             "(shortName(ModRNAString()))")
    }
    ## Check 'transcript_id'.
    if (!is.integer(modifications$transcript_id)
        || any(is.na(modifications$transcript_id)))
        stop("'modifications$transcript_id' must be a integer vector with no ",
             "NAs")
    ## Check 'transcript_name'.
    if (has_col(modifications, "transcript_name")){
        if(!.is_character_or_factor(modifications$transcript_name)
           || any(is.na(modifications$transcript_name))){
            stop("'modifications$transcript_name' must be a character vector ",
                 "(or factor) with no NAs")
        }
    } else {
        modifications$transcript_name <- character(1)
    }
    ## Check 'ensembltrans'.
    if (has_col(modifications, "ensembltrans")){
        if(!.is_character_or_factor(modifications$ensembltrans)
           || any(is.na(modifications$ensembltrans))){
            stop("'modifications$ensembltrans' must be a character vector",
                 " (or factor) with no NAs")
        }
    } else {
        modifications$ensembltrans <- character(1)
    }
    modifications
}

.makeEpiTxDb_normarg_reactions <- function(reactions, modifications_mod_id){
    if (is.null(reactions)) {
        reactions <- data.frame(mod_id = modifications_mod_id[FALSE],
                                mod_rank = integer(0),
                                reaction_genename = character(0),
                                reaction_ensembl = character(0),
                                reaction_ensembltrans = character(0),
                                reaction_entrezid = character(0),
                                reaction_enzyme = character(0),
                                check.names = FALSE, stringsAsFactors = FALSE)
        return(reactions)
    }
    .REQUIRED_COLS <- c("mod_id", "mod_rank")
    .OPTIONAL_COLS <- c("reaction_genename", "reaction_ensembl",
                        "reaction_ensembltrans", "reaction_entrezid",
                        "reaction_enzyme")
    check_colnames(reactions, .REQUIRED_COLS, .OPTIONAL_COLS, "reactions")
    ## Check 'mod_id'.
    .check_foreign_key(reactions$mod_id, "integer", "reactions$mod_id",
                       modifications_mod_id, "integer", "modifications$mod_id")
    ## Check 'mod_rank'.
    if (!is.numeric(reactions$mod_rank)
        || any(is.na(reactions$mod_rank)))
        stop("'reactions$mod_rank' must be an integer vector ",
             "with no NAs")
    if (!is.integer(reactions$mod_rank))
        reactions$mod_rank <- as.integer(reactions$mod_rank)
    if (any(reactions$mod_rank <= 0L))
        stop("'reactions$mod_rank' contains non-positive values")
    ## Check uniqueness of (mod_id, mod_rank) pairs.
    if (any(S4Vectors:::duplicatedIntegerPairs(reactions$mod_id,
                                               reactions$mod_rank)))
        stop("'reactions' must contain unique (mod_id, mod_rank) pairs")
    ## Check 'reaction_genename'.
    if (has_col(reactions, "reaction_genename")
        && !.is_character_or_factor(reactions$reaction_genename))
        stop("'reactions$reaction_genename' must be a character vector ",
             "(or factor)")
    ## Check 'reaction_ensembl'.
    if (has_col(reactions, "reaction_ensembl")
        && !.is_character_or_factor(reactions$reaction_ensembl))
        stop("'reactions$reaction_ensembl' must be a character vector ",
             "(or factor)")
    ## Check 'reaction_ensembltrans'.
    if (has_col(reactions, "reaction_ensembltrans")
        && !.is_character_or_factor(reactions$reaction_ensembltrans))
        stop("'reactions$reaction_ensembltrans' must be a character vector ",
             "(or factor)")
    ## Check 'reaction_entrezid'.
    if (has_col(reactions, "reaction_entrezid")
        && !.is_character_or_factor(reactions$reaction_entrezid))
        stop("'reactions$reaction_entrezid' must be a character vector ",
             "(or factor)")
    ## Check 'reaction_enzyme'.
    if (has_col(reactions, "reaction_enzyme")
        && !.is_character_or_factor(reactions$reaction_enzyme))
        stop("'reactions$reaction_enzyme' must be a character vector ",
             "(or factor)")
    reactions
}

.makeEpiTxDb_normarg_specifiers <- function(specifier, modifications_mod_id){
    if (is.null(specifier)) {
        specifier <- data.frame(mod_id = modifications_mod_id[FALSE],
                                specifier_type = character(0),
                                specifier_genename = character(0),
                                specifier_entrezid = character(0),
                                specifier_ensembl = character(0),
                                check.names = FALSE, stringsAsFactors = FALSE)
        return(specifier)
    }
    .REQUIRED_COLS <- c("mod_id", "specifier_type", "specifier_genename")
    .OPTIONAL_COLS <- c("specifier_entrezid", "specifier_ensembl")
    check_colnames(specifier, .REQUIRED_COLS, .OPTIONAL_COLS, "specifier")
    ## Check 'mod_id'.
    .check_foreign_key(specifier$mod_id, "integer", "specifier$mod_id",
                       modifications_mod_id, "integer", "modifications$mod_id")
    ## Check 'specifier_type'.
    if (has_col(specifier, "specifier_type")
        && !.is_character_or_factor(specifier$specifier_type))
        stop("'specifier$specifier_type' must be a character vector ",
             "(or factor)")
    ## Check 'specifier_genename'.
    if (has_col(specifier, "specifier_genename")
        && !.is_character_or_factor(specifier$specifier_genename))
        stop("'specifier$specifier_genename' must be a character vector ",
             "(or factor)")
    ## Check 'specifier_entrezid'.
    if (has_col(specifier, "specifier_entrezid")
        && !.is_character_or_factor(specifier$specifier_entrezid))
        stop("'specifier$specifier_entrezid' must be a character vector ",
             "(or factor)")
    ## Check 'specifier_ensembl'.
    if (has_col(specifier, "specifier_ensembl")
        && !.is_character_or_factor(specifier$specifier_ensembl))
        stop("'specifier$specifier_ensembl' must be a character vector ",
             "(or factor)")
    specifier
}

.makeEpiTxDb_normarg_references <- function(references, modifications_mod_id){
    if (is.null(references)) {
        references <- data.frame(mod_id = modifications_mod_id[FALSE],
                                 reference_type = character(0),
                                 reference = character(0),
                                 check.names = FALSE, stringsAsFactors = FALSE)
        return(references)
    }
    .REQUIRED_COLS <- c("mod_id")
    .OPTIONAL_COLS <- c("reference_type", "reference")
    check_colnames(references, .REQUIRED_COLS, .OPTIONAL_COLS, "reference")
    ## Check 'mod_id'.
    .check_foreign_key(references$mod_id, "integer", "references$mod_id",
                       modifications_mod_id, "integer", "modifications$mod_id")
    ## Check 'reference_type'.
    if (has_col(references, "reference_type")){
        if(!.is_character_or_factor(references$reference_type)){
            stop("'references$reference_type' must be a character vector ",
                 "(or factor)")
        }
    } else {
        references$reference_type <- character(1)
    }
    ## Check 'reference'.
    if (has_col(references, "reference")){
        if(!.is_character_or_factor(references$reference)){
            stop("'references$reference' must be a character vector ",
                 "(or factor)")
        }
    } else {
        references$reference <- character(1)
    }
    references
}

.makeEpiTxDb_normarg_metadata <- function(metadata)
{
    if (!is.null(metadata)) {
        if (!is.data.frame(metadata))
            stop("'metadata' must be NULL or a data.frame")
        if (!setequal(colnames(metadata), c("name", "value")))
            stop("'metadata' columns must be \"name\" and \"value\"")
    }
    metadata
}

# id helper functions ----------------------------------------------------------

.make_modifications_internal_mod_id <- function(modifications, reassign.ids)
{
    if (!reassign.ids)
        return(modifications$mod_id)
    makeFeatureIds(name = factor(modifications$mod_name,
                                 unique(modifications$mod_name)),
                   type = factor(modifications$mod_type,
                                 unique(modifications$mod_type)),
                   start = factor(modifications$mod_start,
                                  unique(modifications$mod_start)),
                   end = factor(modifications$mod_end,
                                unique(modifications$mod_end)))
}

# write table helper functions -------------------------------------------------

.write_modifications_table <- function(conn,
                                       mod_internal_id,
                                       mod_type,
                                       mod_name,
                                       mod_start,
                                       mod_end,
                                       transcript_id,
                                       ensembltrans,
                                       entrezid){
    if (is.null(mod_name))
        mod_name <- rep.int(NA_character_, length(mod_internal_id))
    data <- data.frame(mod_internal_id = mod_internal_id,
                       mod_type = mod_type,
                       mod_name = mod_name,
                       mod_start = mod_start,
                       mod_end = mod_end,
                       transcript_id = transcript_id,
                       ensembltrans = ensembltrans,
                       entrezid = entrezid,
                       check.names = FALSE, stringsAsFactors = FALSE)
    data <- unique(data)

    ## Create the table.
    SQL <- build_SQL_CREATE_modification_table()
    dbExecute(conn, SQL)

    ## Fill the table.
    insert_data_into_table(conn, "modification", data)
}

.write_reactions_table <- function(conn,
                                   internal_mod_id,
                                   mod_rank,
                                   reaction_genename,
                                   reaction_ensembl,
                                   reaction_ensembltrans,
                                   reaction_entrezid,
                                   reaction_enzyme){
    data <- data.frame(internal_mod_id = internal_mod_id,
                       mod_rank = mod_rank,
                       reaction_genename = reaction_genename,
                       reaction_ensembl = reaction_ensembl,
                       reaction_ensembltrans = reaction_ensembltrans,
                       reaction_entrezid = reaction_entrezid,
                       reaction_enzyme = reaction_enzyme,
                       check.names = FALSE, stringsAsFactors = FALSE)

    ## Create the 'reaction' table and related indices.
    SQL <- build_SQL_CREATE_reaction_table()
    dbExecute(conn, SQL)

    ## Fill the 'reaction' table.
    insert_data_into_table(conn, "reaction", data)
}

.write_specifiers_table <- function(conn,
                                    internal_mod_id,
                                    specifier_type,
                                    specifier_genename,
                                    specifier_entrezid,
                                    specifier_ensembl){
    data <- data.frame(internal_mod_id = internal_mod_id,
                       specifier_type = specifier_type,
                       specifier_genename = specifier_genename,
                       specifier_entrezid = specifier_entrezid,
                       specifier_ensembl = specifier_ensembl,
                       check.names = FALSE, stringsAsFactors = FALSE)

    ## Create the 'specifier' table and related indices.
    SQL <- build_SQL_CREATE_specifier_table()
    dbExecute(conn, SQL)

    ## Fill the 'specifier' table.
    insert_data_into_table(conn, "specifier", data)
}

.write_references_table <- function(conn,
                                    internal_mod_id,
                                    reference_type,
                                    reference){
    data <- data.frame(internal_mod_id = internal_mod_id,
                       reference_type = reference_type,
                       reference = reference,
                       check.names = FALSE, stringsAsFactors = FALSE)

    ## Create the 'reference' table and related indices.
    SQL <- build_SQL_CREATE_reference_table()
    dbExecute(conn, SQL)

    ## Fill the 'reference' table.
    insert_data_into_table(conn, "reference", data)
}

.write_metadata_table <- function(conn, metadata){
    nb_modifications <- dbEasyQuery(conn,
                                    "SELECT COUNT(*) FROM modification")[[1L]]
    thispkg_version <- packageDescription("EpiTxDb")$Version
    rsqlite_version <- packageDescription("RSQLite")$Version
    mat1 <- matrix(c(
        DB_TYPE_NAME, DB_TYPE_VALUE,
        "Supporting package", "EpiTxDb"),
        ncol = 2, byrow = TRUE
    )
    mat2 <- matrix(c(
        "Nb of transcripts", nb_modifications,
        "Db created by",     "EpiTxDb package from Bioconductor",
        "Creation time",     svn.time(),
        "EpiTxDb version at creation time", thispkg_version,
        "RSQLite version at creation time", rsqlite_version,
        "DBSCHEMAVERSION",   DB_SCHEMA_VERSION),
        ncol = 2, byrow = TRUE
    )
    colnames(mat1) <- colnames(mat2) <- c("name", "value")
    metadata <- rbind(data.frame(name = mat1[ , "name"], value = mat1[ , "value"],
                                 check.names = FALSE, stringsAsFactors = FALSE),
                      metadata,
                      data.frame(name = mat2[ , "name"], value = mat2[ , "value"],
                                 check.names = FALSE, stringsAsFactors = FALSE))
    dbWriteTable(conn, "metadata", metadata, row.names = FALSE)
}


# makeEpiTxDb------------------------------------------------------------------

#' @rdname makeEpiTxDb
#' @export
makeEpiTxDb <- function(modifications, reactions = NULL, specifiers = NULL,
                        references = NULL, metadata = NULL,
                        reassign.ids = FALSE){
    if (!isTRUEorFALSE(reassign.ids))
        stop("'reassign.ids' must be TRUE or FALSE")

    modifications <- .makeEpiTxDb_normarg_modifications(modifications)
    reactions <- .makeEpiTxDb_normarg_reactions(reactions, modifications$mod_id)
    specifiers <- .makeEpiTxDb_normarg_specifiers(specifiers, modifications$mod_id)
    references <- .makeEpiTxDb_normarg_references(references, modifications$mod_id)
    metadata <- .makeEpiTxDb_normarg_metadata(metadata)

    modifications_internal_mod_id <-
        .make_modifications_internal_mod_id(modifications, reassign.ids)
    reactions_internal_mod_id <- translateIds(modifications$mod_id,
                                              modifications_internal_mod_id,
                                              reactions$mod_id)
    specifiers_internal_mod_id <- translateIds(modifications$mod_id,
                                               modifications_internal_mod_id,
                                               specifiers$mod_id)
    references_internal_mod_id <- translateIds(modifications$mod_id,
                                               modifications_internal_mod_id,
                                               references$mod_id)

    ## Create the db in a temp file.
    conn <- dbConnect(SQLite(), dbname="")
    .write_modifications_table(conn,
                               modifications_internal_mod_id,
                               modifications$mod_type,
                               modifications$mod_name,
                               modifications$mod_start,
                               modifications$mod_end,
                               modifications$transcript_id,
                               modifications$transcript_name,
                               modifications$ensembltrans)
    .write_reactions_table(conn,
                           reactions_internal_mod_id,
                           reactions$mod_rank,
                           reactions$reaction_genename,
                           reactions$reaction_ensembl,
                           reactions$reaction_ensembltrans,
                           reactions$reaction_entrezid,
                           reactions$reaction_enzyme)
    .write_specifiers_table(conn,
                            specifiers_internal_mod_id,
                            specifiers$specifier_type,
                            specifiers$specifier_genename,
                            specifiers$specifier_entrezid,
                            specifiers$specifier_ensembl)
    .write_references_table(conn,
                            references_internal_mod_id,
                            references$reference_type,
                            references$reference)
    .write_metadata_table(conn, metadata)
    EpiTxDb(conn)
}
