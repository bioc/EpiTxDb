context("makeEpiTxDbFromRMBase")
test_that("makeEpiTxDbFromRMBase:",{
    #
    expect_type(listAvailableOrganismsFromRMBase(),"character")
    expect_error(listAvailableGenomesFromRMBase("human1"),
                 "'organism' must be a valid organism")
    expect_type(listAvailableGenomesFromRMBase("human"),"character")
    expect_error(listAvailableModFromRMBase("human1"),
                 "'organism' must be a valid organism")
    expect_error(listAvailableModFromRMBase("human"),
                 'argument "genome" is missing')
    expect_error(listAvailableModFromRMBase("human","hg191"),
                 "'genome' must be a valid genome for the")
    expect_type(listAvailableModFromRMBase("human","hg19"),"character")
    #
    # files <- system.file("extdata", "RMBase_testdata.txt.gz",
    #                      package = "EpiTxDb")
    # expect_warning(etdb <- makeEpiTxDbFromRMBaseFiles(files))
    # expect_s4_class(etdb,"EpiTxDb")
    # expect_s3_class(metadata(etdb),"data.frame")
    # #
    # actual <- modifications(etdb)
    # expect_s4_class(actual,"GRanges")
    # expect_length(actual,512L)
    # dbDisconnect(etdb$conn)
})