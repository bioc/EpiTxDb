
library(GenomicRanges)

context("Ranges helpers")
test_that("Ranges helpers:",{
    gr <- GRanges("chr1:1-5:+")
    actual <- positionSequence(gr)
    expect_type(actual,"integer")
    expect_equal(actual,seq.int(1L,5L))
    #
    grl <- GRangesList("1" = gr,"2" = gr,"3" = gr)
    actual2 <- positionSequence(grl)
    expect_named(actual2,c("1","2","3"))
    expect_s4_class(actual2,"IntegerList")
    expect_equal(actual2[[1L]],actual)
    expect_equal(actual2[[2L]],actual)
    expect_equal(actual2[[3L]],actual)
    #
    strand(grl[[2]]) <- "-"
    actual2 <- positionSequence(grl, TRUE, TRUE)
    expect_equal(actual2[[1L]],actual)
    expect_equal(actual2[[2L]],rev(actual))
    expect_equal(actual2[[3L]],actual)
    # rescale
    x <- IRanges("5-10")
    actual <- rescale(x, 100, 10)
    expect_s4_class(actual,"IRanges")
    expect_equal(start(actual), 41L)
    expect_equal(end(actual), 100L)
    #
    actual <- rescale(x, "31-60", "5-14")
    expect_s4_class(actual,"IRanges")
    expect_equal(start(actual), 31L)
    expect_equal(end(actual), 48L)
    #
    x <- IRangesList(x,x)
    actual <- rescale(x, 100, 10)
    expect_s4_class(actual,"IRangesList")
    expect_equal(start(actual[[2L]]), 41L)
    expect_equal(end(actual[[2L]]), 100L)
    #
    actual <- rescale(x, "31-60", "5-14")
    expect_s4_class(actual,"IRangesList")
    expect_equal(start(actual[[2L]]), 31L)
    expect_equal(end(actual[[2L]]), 48L)
})
