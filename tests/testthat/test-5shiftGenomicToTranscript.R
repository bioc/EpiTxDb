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
    expect_warning(actual <- shiftGenomicToTranscript(subject,tx))
    expect_s4_class(actual,"GRangesList")
    #
    actual2 <- shiftTranscriptToGenomic(actual,tx)
    actual3 <- all(ranges(actual2) == ranges(subject[list(1,c(1,1),c(1,2))]))
    expect_true(all(actual3))
})
