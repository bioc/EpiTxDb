
library(AnnotationHub)
library(BSgenome.Scerevisiae.UCSC.sacCer3)

context("makeEpiTxDbFromtRNAdb")
test_that("makeEpiTxDbFromtRNAdb:",{
    httptest::skip_if_disconnected(url = "http://trna.bioinf.uni-leipzig.de/")
    bs <- BSgenome.Scerevisiae.UCSC.sacCer3
    # get tx
    ah <- AnnotationHub()
    edb <- query(ah, c("EnsDb","Saccharomyces cerevisiae", "99"))[[1]]
    seqlevelsStyle(edb) <- "UCSC"
    tx <- exonsBy(edb,"tx")
    tx_id <- IRanges::CharacterList(Map(rep,names(tx),lengths(tx)))
    mcols(tx, level="within")[,"tx_id"] <- tx_id
    genome(tx) <- "sacCer3"
    tx <- tx[lengths(tx) != 0L]
    #
    seq <- getSeq(bs,tx)
    seq <- relist(unlist(unlist(seq)),
                  IRanges::PartitioningByWidth(sum(nchar(seq))))
    seq_rna <- as(seq,"RNAStringSet")
    expect_message(etdb <- makeEpiTxDbFromtRNAdb("Saccharomyces cerevisiae",
                                                 seq_rna))
    expect_s4_class(etdb,"EpiTxDb")
    expect_s3_class(metadata(etdb),"data.frame")
    #
    actual <- modifications(etdb)
    expect_s4_class(actual,"GRanges")
    expect_length(actual,2990L)
    dbDisconnect(etdb$conn)
})
