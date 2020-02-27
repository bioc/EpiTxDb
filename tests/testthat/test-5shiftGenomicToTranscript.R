context("shiftGenomicToTranscript")
test_that("shiftGenomicToTranscript:",{
    subject1 <- GRanges("chr1", IRanges(3, 6),
                        strand = "+")
    subject2 <- GRanges("chr1", IRanges(c(17,23), width=3),
                        strand = c("+","-"))
    subject3 <- GRanges("chr2", IRanges(c(51, 54), c(53, 59)),
                        strand = "-")
    subject <- GRangesList(a=subject1, b=subject2, c=subject3)
    tx1 <- GRanges("chr1", IRanges(1, 40),
                   strand="+")
    tx2 <- GRanges("chr1", IRanges(10, 30),
                   strand="+")
    tx3 <- GRanges("chr2", IRanges(50, 60),
                   strand="-")
    tx <- GRangesList(a=tx1, b=tx2, c=tx3)
    expect_warning(actual <- shiftGenomicToTranscript(unlist(subject),tx))
    expect_s4_class(actual,"GRanges")
    expect_equal(colnames(mcols(actual)),c("seq_start","seq_end","seq_strand",
                                           "seq_name"))
    expect_equal(start(actual),c(3L,17L,8L,8L,2L))
    expect_equal(end(actual),c(6L,19L,10L,10L,7L))
    actual2 <- shiftTranscriptToGenomic(actual,tx)
    expect_s4_class(actual2,"GRanges")
    expect_equal(ranges(actual2),
                 unname(ranges(unlist(subject)[c(1L,2L,2L,4L,5L)])))
    #
    expect_warning(actual <- shiftGenomicToTranscript(subject,tx))
    expect_s4_class(actual,"GRangesList")
    expect_equal(unlist(start(actual), use.names = FALSE),
                 c(3L, 17L, 8L, 8L, 2L))
    expect_equal(unlist(end(actual), use.names = FALSE),
                 c(6L, 19L, 10L, 10L, 7L))
    actual2 <- shiftTranscriptToGenomic(actual,tx)
    expect_s4_class(actual2,"GRangesList")
    expect_equal(unlist(start(actual2), use.names = FALSE),
                 c(3L, 26L, 17L, 51L, 54L))
    expect_equal(unlist(end(actual2), use.names = FALSE),
                 c(6L, 28L, 19L, 53L, 59L))
})
