context("EpiTxDb class")
test_that("EpiTxDb:",{
    expect_error(EpiTxDb:::.format_modifications(),
                 'argument "modifications" is missing')
    expect_error(EpiTxDb:::.format_seqnames(),
                 'argument "modifications" is missing')
    expect_error(EpiTxDb:::.format_reactions(),
                 'argument "modifications" is missing')
    expect_error(EpiTxDb:::.format_specifiers(),
                 'argument "modifications" is missing')
    expect_error(EpiTxDb:::.format_references(),
                 'argument "modifications" is missing')
    #
    actual <- EpiTxDb:::.format_modifications(NULL)
    expect_s3_class(actual,"data.frame")
    expect_equal(colnames(actual),c("mod_id","mod_type","mod_name","mod_start",
                                    "mod_end","sn_id"))
    df <- data.frame(mod_id = 1, mod_type = 1, mod_name = 1, mod_start = 1,
                     mod_end = 1, sn_id = 1)
    actual <- EpiTxDb:::.format_modifications(df)
    expect_type(actual[,1],"double")
    expect_type(actual[,3],"double")
    actual <- EpiTxDb:::.format_modifications(df, set.col.class = TRUE)
    expect_type(actual[,1],"integer")
    expect_type(actual[,3],"character")
    #
    actual <- EpiTxDb:::.format_seqnames(NULL)
    expect_s3_class(actual,"data.frame")
    expect_equal(colnames(actual),c("sn_id","sn_name"))
    df <- data.frame(sn_id = 1, sn_name = 1)
    actual <- EpiTxDb:::.format_seqnames(df)
    expect_type(actual[,1],"double")
    expect_type(actual[,3],"double")
    actual <- EpiTxDb:::.format_seqnames(df, set.col.class = TRUE)
    expect_type(actual[,1],"integer")
    expect_type(actual[,3],"character")
    #
    actual <- EpiTxDb:::.format_reactions(NULL)
    expect_s3_class(actual,"data.frame")
    expect_equal(colnames(actual),c("rx_id","rx_genename","rx_rank",
                                    "rx_ensembl","rx_ensembltrans",
                                    "rx_entrezid"))
    df <- data.frame(rx_id = 1, rx_genename = 1, rx_rank = 1, rx_ensembl = 1,
                     rx_ensembltrans = 1, rx_entrezid = 1)
    actual <- EpiTxDb:::.format_reactions(df)
    expect_type(actual[,1],"double")
    expect_type(actual[,3],"double")
    actual <- EpiTxDb:::.format_reactions(df, set.col.class = TRUE)
    expect_type(actual[,1],"integer")
    expect_type(actual[,2],"character")
    #
    actual <- EpiTxDb:::.format_specifiers(NULL)
    expect_s3_class(actual,"data.frame")
    expect_equal(colnames(actual),c("spec_id","spec_type","spec_genename",
                                    "spec_ensembl","spec_ensembltrans",
                                    "spec_entrezid"))
    df <- data.frame(spec_id = 1, spec_type = 1, spec_genename = 1,
                     spec_ensembl = 1, spec_ensembltrans = 1, spec_entrezid = 1)
    actual <- EpiTxDb:::.format_specifiers(df)
    expect_type(actual[,1],"double")
    expect_type(actual[,3],"double")
    actual <- EpiTxDb:::.format_specifiers(df, set.col.class = TRUE)
    expect_type(actual[,1],"integer")
    expect_true(is.factor(actual[,2]))
    #
    actual <- EpiTxDb:::.format_references(NULL)
    expect_s3_class(actual,"data.frame")
    expect_equal(colnames(actual),c("ref_id","ref_type","ref"))
    df <- data.frame(ref_id = 1, ref_type = 1, ref = 1)
    actual <- EpiTxDb:::.format_references(df)
    expect_type(actual[,1],"double")
    expect_type(actual[,3],"double")
    actual <- EpiTxDb:::.format_references(df, set.col.class = TRUE)
    expect_type(actual[,1],"integer")
    expect_true(is.factor(actual[,2]))
    expect_type(actual[,3],"character")
    ############################################################################
    etdb_file <- system.file("extdata", "EpiTxDb.Hs.hg38.snoRNAdb.sqlite",
                             package="EpiTxDb")
    etdb <- loadDb(etdb_file)
    expect_s4_class(etdb,"EpiTxDb")
    # .get_EpiTxDb_seqinfo
    actual <- EpiTxDb:::.get_EpiTxDb_seqinfo(etdb)
    expect_s4_class(actual,"Seqinfo")
    expect_equal(seqinfo(etdb),actual)
    expect_type(seqlevels(etdb),"character")
    # as.list
    actual <- as.list(etdb)
    expect_type(actual,"list")
    expect_named(actual,c("modifications","seqnames","reactions",
                          "specifiers","references"))
    # compareEpiTxDbs
    expect_true(EpiTxDb:::compareEpiTxDbs(etdb,etdb))
})
