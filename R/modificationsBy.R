#' @include EpiTxDb-class.R
NULL

.format_feature_columns <- function(epitxdb, columns, by){
    abbr <- .makeColAbbreviations(epitxdb)
    user_columns <- abbr[abbr %in% columns]
    names(user_columns) <- gsub("^_","",names(user_columns))
    names(user_columns) <- gsub("^mod_type","mod",names(user_columns))
    user_columns
}

.split_into_GRL <- function(epitxdb, gr, columns, by){
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
                "seqnames" = seqnames(gr),
                "modtype" = gr_mcols$mod,
                "reaction" = FUN(gr_mcols$rx_genename),
                "specifier" = FUN(gr_mcols$spec_genename),
                "specifiertype" = FUN(gr_mcols$spec_type),
                stop("unsupported 'by' value."))
    grl <- S4Vectors::split(gr, f)
    .assignMetadataList(grl, epitxdb)
}


.extract_features_by <- function(epitxdb, by = "seqnames"){
    ans <- switch(by,
                  "seqnames" = modifications(epitxdb,
                                             c("MODID", "MODTYPE","MODNAME")),
                  "modtype" = modifications(epitxdb,
                                             c("MODTYPE", "MODID", "MODNAME")),
                  "reaction" = modifications(epitxdb,
                                             c("RXENSEMBL", "RXENSEMBLTRANS",
                                               "RXENTREZID", "RXGENENAME")),
                  "specifier" = modifications(epitxdb,
                                              c("SPECENSEMBL","SPECENTREZID",
                                                "SPECGENENAME", "SPECTYPE")),
                  "specifiertype" = modifications(epitxdb,
                                                   c("SPECENSEMBL","SPECENTREZID",
                                                     "SPECGENENAME", "SPECTYPE")),
                  stop("unsupported 'by' value."))
    mcolumns <- .format_feature_columns(epitxdb, colnames(mcols(ans)), by)
    .split_into_GRL(epitxdb, ans, mcolumns, by)
}

#' @rdname modifications
#' @export
setMethod("modificationsBy", "EpiTxDb",
    function(x, by = c("seqnames","modtype","reaction","specifier",
                       "specifiertype")){
      by <- match.arg(by)
      .extract_features_by(x, by = by)
    }
)


#' @rdname modifications
#' @export
setMethod("modifiedSeqsByTranscript", c("EpiTxDb","DNAStringSet"),
    function(x, sequences){
      grl <- modificationsBy(x, "seqnames")
      combineIntoModstrings(sequences, grl)
    }
)
