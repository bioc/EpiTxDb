
library(GenomicRanges)

context("makeEpiTxDbFromGRanges")
test_that("makeEpiTxDbFromGRanges:",{
    gr <- GRanges(seqnames = "test",
                  ranges = IRanges::IRanges(1,1),
                  strand = "+",
                  DataFrame(mod_id = 1L,
                            mod_type = "Am",
                            mod_name = "Am_1"))
    etdb <- makeEpiTxDbFromGRanges(gr)
    expect_s4_class(etdb,"EpiTxDb")
    expect_s3_class(metadata(etdb),"data.frame")
    dbDisconnect(etdb$conn)
})
