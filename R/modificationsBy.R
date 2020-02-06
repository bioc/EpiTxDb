#' @include TxModDb-class.R
NULL

.format_feature_columns <- function(txmoddb, columns, by){
    abbr <- .makeColAbbreviations(txmoddb)
    user_columns <- abbr[abbr %in% columns]
    names(user_columns) <- gsub("^_","",names(user_columns))
    names(user_columns) <- gsub("^mod_type","mod",names(user_columns))
    user_columns
}

.split_into_GRL <- function(txmoddb, gr, columns, by, use.names){
    gr_mcols <- mcols(gr)
    colnames(gr_mcols)[match(columns,colnames(gr_mcols))] <- names(columns)
    mcols(gr) <- gr_mcols
    FUN <- function(x){
        x <- unique(x)
        x <- vapply(x,paste,character(1),collapse = "-")
        x[x == ""] <- "unknown"
        x
    }
    f <- switch(by,
                "transcript" = seqnames(gr),
                "type" = gr_mcols$mod,
                "reaction" = FUN(gr_mcols$reaction_enzyme),
                "specifier" = FUN(gr_mcols$specifier_type),
                stop("unsupported 'by' value."))
    grl <- S4Vectors::split(gr, f)
    .assignMetadataList(grl, txmoddb)
}


.extract_features_by <- function(txmoddb, by = "transcript"){
    ans <- switch(by,
                  "transcript" = modifications(txmoddb,
                                               c("MODID", "MODTYPE","MODNAME")),
                  "type" = modifications(txmoddb,
                                         c("MODTYPE", "MODID", "MODNAME")),
                  "reaction" = modifications(txmoddb,
                                             c("RXENSEMBL", "RXENSEMBLTRANS",
                                               "RXENTREZID", "RXENZYME",
                                               "RXGENENAME")),
                  "specifier" = modifications(txmoddb,
                                              c("SPECENSEMBL","SPECENTREZID",
                                                "SPECGENENAME", "SPECTYPE")),
                  stop("unsupported 'by' value."))
    mcolumns <- .format_feature_columns(txmoddb, colnames(mcols(ans)), by)
    .split_into_GRL(txmoddb, ans, mcolumns, by, use.names)
}

#' @rdname modifications
#' @export
setMethod("modificationsBy", "TxModDb",
    function(x, by = c("transcript","type","reaction","specifier")){
      by <- match.arg(by)
      .extract_features_by(x, by = by)
    }
)


#' @rdname modifications
#' @export
setMethod("modifiedSeqsByTranscript", c("TxModDb","DNAStringSet"),
    function(x, sequences){
      grl <- modificationsBy(x, "transcript")
      combineIntoModstrings(sequences, grl)
    }
)
