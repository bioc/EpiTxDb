#' @include EpiTxDb-class.R
NULL

#' @name makeEpiTxDb
#'
#' @title Creating a \code{EpiTxDb} from user supplied annotations as
#' \code{data.frame}s
#'
#' @description
#' \code{makeEpiTxDb} is a low-level constructor for creating a
#' \code{\link[=EpiTxDb-class]{EpiTxDb}} object from user supplied annotations.
#'
#' This functions typically will not be used by regular users.
#'
#' @param modifications A \code{data.frame} containg the following columns:
#'   \itemize{
#'     \item{\code{mod_id}:} {a unique \code{integer} value per modification.}
#'     \item{\code{mod_type}:} {the modification type as a \code{character} or
#'     \code{factor} value. Must be a value from
#'     \code{shortName(ModRNAString())}.}
#'     \item{\code{mod_name}:} {a \code{character} or \code{factor} name for the
#'     specific modification}
#'     \item{\code{mod_start}:} {the start position for the modification as
#'     \code{integer} value. Usually \code{mod_start = mod_end} }
#'     \item{\code{mod_end}:} {the end position for the modification as
#'     \code{integer} value. Usually \code{mod_start = mod_end} }
#'     \item{\code{mod_strand}:} {the strand information for the modificaion as a
#'     \code{character} or \code{factor}. }
#'     \item{\code{sn_id}:} {a \code{integer} value per unique sequence }
#'     \item{\code{sn_name}:} {a \code{character} or \code{factor} as sequence
#'     name, e.g a chromosome or a transcript identifier like \code{chr1}.}
#'   }
#'   The first six are mandatory, whereas one of the last two has to be set.
#'   \code{sn_id} will be generated from \code{sn_name}, if \code{sn_id} is not
#'   set.
#' @param reactions An optional \code{data.frame} containg the following
#'   columns:
#'   \itemize{
#'     \item{\code{mod_id}:} {a \code{integer} value per modification and the
#'     link to the \code{modification} \code{data.frame}.}
#'     \item{\code{rx_genename}:} {a \code{character} or \code{factor} referencing
#'     a genename for the enzyme incorporating the modification. }
#'     \item{\code{rx_rank}:} {a \code{integer} for sorting enzyme reactions, if
#'     multiple enzymes are involved in the modification's
#'     incorporation/maintenance.}
#'     \item{\code{rx_ensembl}:} {a \code{character} or \code{factor} with an
#'     ensembl identifier for the genename of the enzyme.}
#'     \item{\code{rx_ensembltrans}:} {a \code{character} or \code{factor} with an
#'     ensembl identifier for the transcript being translated into the enzyme.}
#'     \item{\code{rx_entrezid}:} {a \code{character} or \code{factor} with an
#'     entrezid for the genename of the enzyme.}
#'   }
#'   (default: \code{reactions = NULL})
#' @param specifiers An optional \code{data.frame} containg the following
#'   columns:
#'     \itemize{
#'     \item{\code{mod_id}:} {a \code{integer} value per modification and the
#'     link to the \code{modification} \code{data.frame}.}
#'     \item{\code{spec_type}:} {a \code{character} or \code{factor}
#'     referencing a type of specifier, e.g. \code{snoRNA}. Not checked for
#'     validity.}
#'     \item{\code{spec_genename}:} {a \code{character} or \code{factor}
#'     referencing a genename for the specifier directing an enzyme to the
#'     specific location for the modification to be incorporated.}
#'     \item{\code{spec_ensembl}:} {a \code{character} or \code{factor} with an
#'     ensembl identifier for the genename of the specifier.}
#'     \item{\code{spec_ensembltrans}:} {a \code{character} or \code{factor} with an
#'     ensembl identifier for the transcript being translated into the specifier.}
#'     \item{\code{spec_entrezid}:} {a \code{character} or \code{factor} with an
#'     entrezid for the genename of the specifier.}
#'   }
#'   (default: \code{specifiers = NULL})
#' @param references An optional \code{data.frame} containg the following
#'   columns:
#'   \itemize{
#'     \item{\code{mod_id}:} {a \code{integer} value per modification and the
#'     link to the \code{modification} \code{data.frame}.}
#'     \item{\code{ref_type}:} {a \code{character} or \code{factor} with a
#'     reference type, e.g. \code{PMID}. Is not checked for validity. }
#'     \item{\code{ref}:} {a \code{character} or \code{factor} with a reference
#'     value, e.g. a specific pubmed id or an journal article. Is not checked
#'     for validity.}
#'   }
#'   (default: \code{references = NULL})
#' @param metadata An optional \code{data.frame} containg the following columns:
#'   \itemize{
#'     \item{\code{name}:} {a \code{character} value used as name}
#'     \item{\code{value}:} {a \code{character} value}
#'   }
#'   This dataframe will be returned
#'   by \code{metadata()} (default: \code{metadata = NULL})
#' @param reassign.ids \code{TRUE} or \code{FALSE} Controls how internal
#'   \code{mod_id}s should be assigned. If \code{reassign.ids} is \code{FALSE}
#'   (the default) and if the ids are supplied, then they are used as the
#'   internal ids, otherwise the internal ids are assigned in a way that is
#'   compatible with the order defined by ordering the modifications first by
#'   chromosome, then by strand, then by start, and finally by end.
#'
#' @seealso
#' \itemize{
#'   \item{\code{\link[=makeEpiTxDbFromGRanges]{makeEpiTxDbFromGRanges}} for
#'     creating a \code{EpiTxDb} object from a
#'   \code{\link[GenomicRanges:GRanges-class]{GRanges}} object and it's metadata
#'     columns}
#'   \item{\code{\link[=makeEpiTxDbFromRMBase]{makeEpiTxDbFromRMBase}} for
#'     creating a \code{EpiTxDb} object from RMBase online resources}
#'   \item{\code{\link[=makeEpiTxDbFromtRNAdb]{makeEpiTxDbFromtRNAdb}} for
#'     creating a \code{EpiTxDb} object from tRNAdb online resources}
#'   \item{\code{\link[Modstrings:alphabet]{shortName}} and
#'   \code{\link[Modstrings:ModRNAString]{ModRNAString}} for information on
#'   \code{ModRNAString} objects.}
#' }
#'
#' @return a \code{EpiTxDb} object.
#'
#' @export
#'
#' @examples
#' mod <- data.frame("mod_id" = 1L,
#'                   "mod_type" = "m1A",
#'                   "mod_name" = "m1A_1",
#'                   "mod_start" = 1L,
#'                   "mod_end" = 1L,
#'                   "mod_strand" = "+",
#'                   "sn_id" = 1L,
#'                   "sn_name" = "test")
#' rx <- data.frame(mod_id = 1L,
#'                  rx_genename = "test",
#'                  rx_rank = 1L,
#'                  rx_ensembl = "test",
#'                  rx_ensembltrans = "test",
#'                  rx_entrezid = "test")
#' spec <- data.frame(mod_id = 1L,
#'                    spec_type = "test",
#'                    spec_genename = "test",
#'                    spec_ensembl = "test",
#'                    spec_ensembltrans = "test",
#'                    spec_entrezid = "test")
#' ref <- data.frame(mod_id = 1L,
#'                   ref_type = "test",
#'                   ref = "test")
#' etdb <- makeEpiTxDb(mod,rx,spec,ref)
NULL

# helper functions -------------------------------------------------------------

insert_data_into_table <- function(conn, table, data){
    stopifnot(is.list(data))
    placeholders <- paste(rep.int("?", length(data)), collapse = ",")
    SQL <- sprintf("INSERT INTO %s VALUES (%s)", table, placeholders)
    ## dbExecute() emits annoying warnings if 'params' is a named list or if
    ## some of its list elements are factors.
    params <- unname(as.list(data))
    params <- lapply(params,
                     function(x) if (is.factor(x)) as.character(x) else x)
    dbExecute(conn, SQL, params = params)
}

duplicatedIntegerQuads <- S4Vectors:::duplicatedIntegerQuads
orderIntegerQuads <- S4Vectors:::orderIntegerQuads
matchIntegerQuads <- S4Vectors:::matchIntegerQuads

makeFeatureIds <- function(name = NULL, type, start, end,
                           same.id.for.dups = FALSE){
    a <- type
    b <- name
    c <- start
    d <- end
    if (!same.id.for.dups) {
        oo <- orderIntegerQuads(a, b, c, d)
        ans <- integer(length(oo))
        ans[oo] <- seq_len(length(oo))
        return(ans)
    }
    ## There should be a better way to do this...
    is_not_dup <- !duplicatedIntegerQuads(a, b, c ,d)
    ua <- a[is_not_dup]
    ub <- b[is_not_dup]
    uc <- c[is_not_dup]
    ud <- d[is_not_dup]
    oo <- orderIntegerQuads(ua, ub, uc, ud)
    ua <- ua[oo]
    ub <- ub[oo]
    uc <- uc[oo]
    ud <- ud[oo]
    matchIntegerQuads(a, b, c, d, ua, ub, uc, ud)
}

# check helper functions -------------------------------------------------------

.is_character_or_factor <- GenomicFeatures:::.is_character_or_factor
.check_foreign_key <- GenomicFeatures:::.check_foreign_key
translateIds <- GenomicFeatures:::translateIds
check_colnames <- GenomicFeatures:::check_colnames
has_col <- GenomicFeatures:::has_col
dbEasyQuery <- GenomicFeatures:::dbEasyQuery

.makeEpiTxDb_normarg_modifications <- function(modifications){
    .REQUIRED_COLS <- c("mod_id","mod_type","mod_start","mod_end","mod_strand",
                        "sn_id")
    .OPTIONAL_COLS <- c("mod_name", "sn_name")
    # make sure 'sn_id is set'
    if(!has_col(modifications, "sn_id") &&
       !has_col(modifications, "sn_name")){
        stop("'modifications' must contain a column 'sn_name', if no column ",
             "'sn_id' is given.")
    }
    if(!has_col(modifications, "sn_id") &&
       has_col(modifications, "sn_name")){
        modifications$sn_id <- as.integer(factor(modifications$sn_name,
                                                 unique(modifications$sn_name)))
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
    if (!is.integer(modifications$mod_start)
        || any(is.na(modifications$mod_start)))
        stop("'modifications$mod_start' must be an integer vector ",
             "with no NAs")
    ## Check 'mod_end'.
    if (!is.integer(modifications$mod_end)
        || any(is.na(modifications$mod_end)))
        stop("'modifications$mod_end' must be an integer vector ",
             "with no NAs")
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
    ## Check 'mod_strand'.
    if (!.is_character_or_factor(modifications$mod_strand) ||
        any(!(modifications$mod_strand %in% c("+","-","*"))))
        stop("'modifications$mod_strand' must be a character vector ",
             "(or factor) and contain only '+', '-' or '*'")
    ## Check 'sn_id'.
    if (!is.integer(modifications$sn_id)
        || any(is.na(modifications$sn_id)))
        stop("'modifications$sn_id' must be a integer vector with no ",
             "NAs")
    ## Check 'sn_name'.
    if (has_col(modifications, "sn_name")){
        if(!.is_character_or_factor(modifications$sn_name)
           || any(is.na(modifications$sn_name))){
            stop("'modifications$sn_name' must be a character vector ",
                 "(or factor) with no NAs")
        }
        if(length(unique(modifications$sn_id)) !=
           length(unique(modifications$sn_name))){
            stop("'modifications$sn_name' must match ",
                 "'modifications$sn_id' in meaning.")
        }
        f1 <- factor(modifications$sn_id,
                     unique(modifications$sn_id))
        p1 <- PartitioningByWidth(split(modifications$sn_id,f1))
        f2 <- factor(modifications$sn_name,
                     unique(modifications$sn_name))
        p2 <- PartitioningByWidth(split(modifications$sn_name,f2))
        if(!all(p1 == p2)){
            stop("'modifications$sn_name' must match ",
                 "'modifications$sn_id' in meaning.")
        }
    } else {
        modifications$sn_name <- character(1)
    }
    modifications
}

.makeEpiTxDb_normarg_reactions <- function(reactions, modifications_mod_id){
    if (is.null(reactions)) {
        reactions <- data.frame(mod_id = modifications_mod_id[FALSE],
                                rx_genename = character(0),
                                rx_rank = integer(0),
                                rx_ensembl = character(0),
                                rx_ensembltrans = character(0),
                                rx_entrezid = character(0),
                                check.names = FALSE, stringsAsFactors = FALSE)
        return(reactions)
    }
    .REQUIRED_COLS <- c("mod_id", "rx_genename", "rx_rank")
    .OPTIONAL_COLS <- c("rx_ensembl", "rx_ensembltrans", "rx_entrezid")
    check_colnames(reactions, .REQUIRED_COLS, .OPTIONAL_COLS, "reactions")
    ## Check 'mod_id'.
    .check_foreign_key(reactions$mod_id, "integer", "reactions$mod_id",
                       modifications_mod_id, "integer", "modifications$mod_id")
    ## Check 'rx_rank'.
    if (!is.integer(reactions$rx_rank)
        || any(is.na(reactions$rx_rank)))
        stop("'reactions$rx_rank' must be an integer vector ",
             "with no NAs")
    if (any(reactions$rx_rank <= 0L))
        stop("'reactions$rx_rank' contains non-positive values")
    ## Check uniqueness of (mod_id, rx_rank) pairs.
    if (any(S4Vectors::duplicatedIntegerPairs(reactions$mod_id,
                                              reactions$rx_rank)))
        stop("'reactions' must contain unique (mod_id, rx_rank) pairs")
    ## Check 'rx_genename'.
    if (has_col(reactions, "rx_genename")
        && !.is_character_or_factor(reactions$rx_genename))
        stop("'reactions$rx_genename' must be a character vector ",
             "(or factor)")
    ## Check 'rx_ensembl'.
    if (has_col(reactions, "rx_ensembl")){
        if(!.is_character_or_factor(reactions$rx_ensembl))
            stop("'reactions$rx_ensembl' must be a character vector ",
                 "(or factor)")
    } else {
        reactions$rx_ensembl <- character(1)
    }
    ## Check 'rx_ensembltrans'.
    if (has_col(reactions, "rx_ensembltrans")){
        if(!.is_character_or_factor(reactions$rx_ensembltrans))
            stop("'reactions$rx_ensembltrans' must be a character vector ",
                 "(or factor)")
    } else {
        reactions$rx_ensembltrans <- character(1)
    }
    ## Check 'rx_entrezid'.
    if (has_col(reactions, "rx_entrezid")){
        if(!.is_character_or_factor(reactions$rx_entrezid))
            stop("'reactions$rx_entrezid' must be a character vector ",
                 "(or factor)")
    } else {
        reactions$rx_entrezid <- character(1)
    }
    reactions
}

.makeEpiTxDb_normarg_specifiers <- function(specifier, modifications_mod_id){
    if (is.null(specifier)) {
        specifier <- data.frame(mod_id = modifications_mod_id[FALSE],
                                spec_type = character(0),
                                spec_genename = character(0),
                                spec_ensembl = character(0),
                                spec_ensembltrans = character(0),
                                spec_entrezid = character(0),
                                check.names = FALSE, stringsAsFactors = FALSE)
        return(specifier)
    }
    .REQUIRED_COLS <- c("mod_id", "spec_type", "spec_genename")
    .OPTIONAL_COLS <- c("spec_ensembl","spec_ensembltrans","spec_entrezid")
    check_colnames(specifier, .REQUIRED_COLS, .OPTIONAL_COLS, "specifier")
    ## Check 'mod_id'.
    .check_foreign_key(specifier$mod_id, "integer", "specifier$mod_id",
                       modifications_mod_id, "integer", "modifications$mod_id")
    ## Check 'spec_type'.
    if (has_col(specifier, "spec_type")
        && !.is_character_or_factor(specifier$spec_type))
        stop("'specifier$spec_type' must be a character vector ",
             "(or factor)")
    ## Check 'spec_genename'.
    if (has_col(specifier, "spec_genename")
        && !.is_character_or_factor(specifier$spec_genename))
        stop("'specifier$spec_genename' must be a character vector ",
             "(or factor)")
    ## Check 'spec_ensembl'.
    if (has_col(specifier, "spec_ensembl")){
        if(!.is_character_or_factor(specifier$spec_ensembl))
            stop("'specifier$spec_ensembl' must be a character vector ",
                 "(or factor)")
    } else {
        specifier$spec_ensembl <- character(1)
    }
    ## Check 'spec_ensembltrans'.
    if (has_col(specifier, "spec_ensembltrans")){
        if(!.is_character_or_factor(specifier$spec_ensembltrans))
            stop("'specifier$spec_ensembltrans' must be a character vector ",
                 "(or factor)")
    } else {
        specifier$spec_ensembltrans <- character(1)
    }
    ## Check 'spec_entrezid'.
    if (has_col(specifier, "spec_entrezid")){
        if(!.is_character_or_factor(specifier$spec_entrezid))
            stop("'specifier$spec_entrezid' must be a character vector ",
                 "(or factor)")
    } else {
        specifier$spec_entrezid <- character(1)
    }
    specifier
}

.makeEpiTxDb_normarg_references <- function(references, modifications_mod_id){
    if (is.null(references)) {
        references <- data.frame(mod_id = modifications_mod_id[FALSE],
                                 ref_type = character(0),
                                 ref = character(0),
                                 check.names = FALSE, stringsAsFactors = FALSE)
        return(references)
    }
    .REQUIRED_COLS <- c("mod_id")
    .OPTIONAL_COLS <- c("ref_type", "ref")
    check_colnames(references, .REQUIRED_COLS, .OPTIONAL_COLS, "reference")
    ## Check 'mod_id'.
    .check_foreign_key(references$mod_id, "integer", "references$mod_id",
                       modifications_mod_id, "integer", "modifications$mod_id")
    ## Check 'ref_type'.
    if (has_col(references, "ref_type")){
        if(!.is_character_or_factor(references$ref_type)){
            stop("'references$ref_type' must be a character vector ",
                 "(or factor)")
        }
    } else {
        references$ref_type <- character(1)
    }
    ## Check 'ref'.
    if (has_col(references, "ref")){
        if(!.is_character_or_factor(references$ref)){
            stop("'references$ref' must be a character vector ",
                 "(or factor)")
        }
    } else {
        references$ref <- character(1)
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

.make_ids <- function(df, cols){
    df[is.na(df)] <- ""
    id <- factor(do.call(paste,as.list(df[,colnames(df) %in% cols])))
    as.integer(id)
}

.make_reactions_ids <- function(reactions){
    cols <- c("rx_genename","rx_rank","rx_ensembl",
              "rx_ensembltrans","rx_entrezid")
    reactions$rx_id <- .make_ids(reactions, cols)
    reactions
}

.make_specifiers_ids <- function(specifiers){
    cols <- c("spec_type","spec_genename","spec_ensembl",
              "spec_ensembltrans","spec_entrezid")
    specifiers$spec_id <- .make_ids(specifiers, cols)
    specifiers
}

.make_references_ids <- function(references){
    cols <- c("ref_type","ref")
    references$ref_id <- .make_ids(references, cols)
    references
}

.shrink_df <- function(df, col){
    df[!duplicated(df[,colnames(df) %in% col]),]
}

# write table helper functions -------------------------------------------------

.write_seqnames_table <- function(conn,
                                  sn_id,
                                  sn_name){
    data <- data.frame(sn_id = sn_id,
                       sn_name = sn_name,
                       check.names = FALSE, stringsAsFactors = FALSE)
    data <- unique(data)

    ## Create the table.
    SQL <- build_SQL_CREATE_seqnames_table()
    dbExecute(conn, SQL)

    ## Fill the table.
    insert_data_into_table(conn, "seqnames", data)
}

.write_modification_table <- function(conn,
                                      mod_internal_id,
                                      mod_type,
                                      mod_name,
                                      mod_start,
                                      mod_end,
                                      mod_strand,
                                      sn_id){
    if (is.null(mod_name))
        mod_name <- rep.int(NA_character_, length(mod_internal_id))
    data <- data.frame(mod_internal_id = mod_internal_id,
                       mod_type = mod_type,
                       mod_name = mod_name,
                       mod_start = mod_start,
                       mod_end = mod_end,
                       mod_strand = mod_strand,
                       sn_id = sn_id,
                       check.names = FALSE, stringsAsFactors = FALSE)
    data <- unique(data)

    ## Create the table.
    SQL <- build_SQL_CREATE_modification_table()
    dbExecute(conn, SQL)

    ## Fill the table.
    insert_data_into_table(conn, "modification", data)
}

.write_reaction_table <- function(conn,
                                  rx_id,
                                  rx_genename,
                                  rx_rank,
                                  rx_ensembl,
                                  rx_ensembltrans,
                                  rx_entrezid){
    data <- data.frame(rx_id = rx_id,
                       rx_genename = rx_genename,
                       rx_rank = rx_rank,
                       rx_ensembl = rx_ensembl,
                       rx_ensembltrans = rx_ensembltrans,
                       rx_entrezid = rx_entrezid,
                       check.names = FALSE, stringsAsFactors = FALSE)

    ## Create the 'reaction' table and related indices.
    SQL <- build_SQL_CREATE_reaction_table()
    dbExecute(conn, SQL)

    ## Fill the 'reaction' table.
    insert_data_into_table(conn, "reaction", data)
}

.write_specifier_table <- function(conn,
                                   spec_id,
                                   spec_type,
                                   spec_genename,
                                   spec_ensembl,
                                   spec_ensembltrans,
                                   spec_entrezid){
    data <- data.frame(spec_id = spec_id,
                       spec_type = spec_type,
                       spec_genename = spec_genename,
                       spec_ensembl = spec_ensembl,
                       spec_ensembltrans = spec_ensembltrans,
                       spec_entrezid = spec_entrezid,
                       check.names = FALSE, stringsAsFactors = FALSE)

    ## Create the 'specifier' table and related indices.
    SQL <- build_SQL_CREATE_specifier_table()
    dbExecute(conn, SQL)

    ## Fill the 'specifier' table.
    insert_data_into_table(conn, "specifier", data)
}

.write_reference_table <- function(conn,
                                   ref_id,
                                   ref_type,
                                   ref){
    data <- data.frame(ref_id = ref_id,
                       ref_type = ref_type,
                       ref = ref,
                       check.names = FALSE, stringsAsFactors = FALSE)

    ## Create the 'reference' table and related indices.
    SQL <- build_SQL_CREATE_reference_table()
    dbExecute(conn, SQL)

    ## Fill the 'reference' table.
    insert_data_into_table(conn, "reference", data)
}

.write_modification2reaction_table <- function(conn,
                                                mod_id,
                                                rx_id){
    data <- data.frame(mod_id = mod_id,
                       rx_id = rx_id,
                       check.names = FALSE, stringsAsFactors = FALSE)

    ## Create the 'modification2reaction' table and related indices.
    SQL <- build_SQL_CREATE_modification2reaction_table()
    dbExecute(conn, SQL)

    ## Fill the 'modification2reaction' table.
    insert_data_into_table(conn, "modification2reaction", data)
}

.write_modification2specifier_table <- function(conn,
                                                mod_id,
                                                spec_id){
    data <- data.frame(mod_id = mod_id,
                       spec_id = spec_id,
                       check.names = FALSE, stringsAsFactors = FALSE)

    ## Create the 'modification2specifier' table and related indices.
    SQL <- build_SQL_CREATE_modification2specifier_table()
    dbExecute(conn, SQL)

    ## Fill the 'modification2specifier' table.
    insert_data_into_table(conn, "modification2specifier", data)
}

.write_modification2reference_table <- function(conn,
                                                mod_id,
                                                ref_id){
    data <- data.frame(mod_id = mod_id,
                       ref_id = ref_id,
                       check.names = FALSE, stringsAsFactors = FALSE)

    ## Create the 'modification2reference' table and related indices.
    SQL <- build_SQL_CREATE_modification2reference_table()
    dbExecute(conn, SQL)

    ## Fill the 'modification2reference' table.
    insert_data_into_table(conn, "modification2reference", data)
}

#' @importFrom utils packageDescription
.write_metadata_table <- function(conn, metadata){
    nb_modifications <- dbEasyQuery(conn,
                                    "SELECT COUNT(*) FROM modification")[[1L]]
    thispkg_version <- utils::packageDescription("EpiTxDb")$Version
    rsqlite_version <- utils::packageDescription("RSQLite")$Version
    mat1 <- matrix(c(
        DB_TYPE_NAME, DB_TYPE_VALUE,
        "Supporting package", "EpiTxDb"),
        ncol = 2, byrow = TRUE
    )
    mat2 <- matrix(c(
        "Nb of modifications", nb_modifications,
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
    message("Creating EpiTxDb object ... ", appendLF = FALSE)
    if (!isTRUEorFALSE(reassign.ids)){
        stop("'reassign.ids' must be TRUE or FALSE")
    }
    #
    modifications <- .makeEpiTxDb_normarg_modifications(modifications)
    reactions <- .makeEpiTxDb_normarg_reactions(reactions, modifications$mod_id)
    specifiers <- .makeEpiTxDb_normarg_specifiers(specifiers, modifications$mod_id)
    references <- .makeEpiTxDb_normarg_references(references, modifications$mod_id)
    metadata <- .makeEpiTxDb_normarg_metadata(metadata)
    #
    internal_mod_id <- .make_modifications_internal_mod_id(modifications,
                                                      reassign.ids)
    mod2rx_mod_id <- translateIds(modifications$mod_id,
                                  internal_mod_id,
                                  reactions$mod_id)
    mod2spec_mod_id <- translateIds(modifications$mod_id,
                                    internal_mod_id,
                                  specifiers$mod_id)
    mod2ref_mod_id <- translateIds(modifications$mod_id,
                                   internal_mod_id,
                                  references$mod_id)
    #
    reactions <- .make_reactions_ids(reactions)
    specifiers <- .make_specifiers_ids(specifiers)
    references <- .make_references_ids(references)
    mod2rx_rx_id <- reactions$rx_id
    mod2spec_spec_id <- specifiers$spec_id
    mod2ref_ref_id <- references$ref_id
    reactions <- .shrink_df(reactions, "rx_id")
    specifiers <- .shrink_df(specifiers, "spec_id")
    references <- .shrink_df(references, "ref_id")
    ## Create the db in a temp file.
    conn <- dbConnect(SQLite(), dbname="")
    .write_seqnames_table(conn,
                          modifications$sn_id,
                          modifications$sn_name)
    .write_modification_table(conn,
                              internal_mod_id,
                              modifications$mod_type,
                              modifications$mod_name,
                              modifications$mod_start,
                              modifications$mod_end,
                              modifications$mod_strand,
                              modifications$sn_id)
    .write_reaction_table(conn,
                          reactions$rx_id,
                          reactions$rx_genename,
                          reactions$rx_rank,
                          reactions$rx_ensembl,
                          reactions$rx_ensembltrans,
                          reactions$rx_entrezid)
    .write_specifier_table(conn,
                           specifiers$spec_id,
                           specifiers$spec_type,
                           specifiers$spec_genename,
                           specifiers$spec_ensembl,
                           specifiers$spec_ensembltrans,
                           specifiers$spec_entrezid)
    .write_reference_table(conn,
                           references$ref_id,
                           references$ref_type,
                           references$ref)
    # write link tables
    .write_modification2reaction_table(conn,
                                       mod2rx_mod_id,
                                       mod2rx_rx_id)

    .write_modification2specifier_table(conn,
                                        mod2spec_mod_id,
                                        mod2spec_spec_id)

    .write_modification2reference_table(conn,
                                        mod2ref_mod_id,
                                        mod2ref_ref_id)
    .write_metadata_table(conn, metadata)
    ans <- EpiTxDb(conn)
    message("done")
    ans
}
