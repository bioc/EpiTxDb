context("EpiTxDb select")
test_that("EpiTxDb select:",{
    etdb_file <- system.file("extdata", "EpiTxDb.Hs.hg38.snoRNAdb.sqlite",
                             package="EpiTxDb")
    etdb <- loadDb(etdb_file)
    #
    actual <- EpiTxDb:::.getTableColMapping(etdb)
    expect_type(actual,"list")
    expect_type(actual[[1L]],"character")
    expect_named(actual, c("modification","modification2reaction",
                           "modification2reference","modification2specifier",
                           "reaction","reference","seqnames","specifier"))
    #
    actual <- EpiTxDb:::.makeColAbbreviations(etdb)
    expect_type(actual,"character")
    expect_named(actual)
    #
    expect_error(EpiTxDb:::.reverseColAbbreviations(etdb),
                 'argument "cnames" is missing')
    expect_equal(EpiTxDb:::.reverseColAbbreviations(etdb, "_mod_id"),
                 NA_character_)
    actual <- EpiTxDb:::.reverseColAbbreviations(etdb, c("MODID","MODNAME"))
    expect_equal(actual,c("_mod_id","mod_name"))
    #
    #.getTableNames
    expect_error(EpiTxDb:::.getTableNames(etdb),
                 'argument "cnames" is missing')
    actual <- EpiTxDb:::.getTableNames(etdb, "_mod_id")
    expect_type(actual,"list")
    expect_true(all(is.na(actual[[1L]])))
    actual <- EpiTxDb:::.getTableNames(etdb, c("MODID","MODNAME"))
    expect_type(actual,"list")
    expect_named(actual,c("_mod_id","mod_name"))
    expect_equal(actual[[2L]],"modification")
    #
    expect_error(EpiTxDb:::.getSimpleTableNames(etdb),
                 'argument "cnames" is missing')
    expect_equal(EpiTxDb:::.getSimpleTableNames(etdb, "_mod_id"),NA_character_)
    actual <- EpiTxDb:::.getSimpleTableNames(etdb,
                                             c("MODID","MODNAME","RXGENENAME"))
    expect_type(actual, "character")
    expect_equal(actual[c(1,5)], c("modification","reaction"))
    #
    expect_error(EpiTxDb:::.encodeSortedTableKey(),
                 'argument "sTNames" is missing')
    expect_equal(EpiTxDb:::.encodeSortedTableKey(c("seq","mod")),"modseq")
    expect_equal(EpiTxDb:::.encodeSortedTableKey(c("spec","seq","mod")),
                 "modseq")
    expect_equal(EpiTxDb:::.encodeSortedTableKey(c("spe","seq","mod")),
                 "modseqspe")
    #
    expect_error(EpiTxDb:::.makeTableKey(etdb),
                 'argument "cnames" is missing')
    actual <- EpiTxDb:::.makeTableKey(etdb,c("MODID","MODNAME"))
    expect_equal(actual,"mod")
    actual <- EpiTxDb:::.makeTableKey(etdb,c("MODID","MODNAME","RXGENENAME"))
    expect_equal(actual,"modrea")
    actual <- EpiTxDb:::.makeTableKey(etdb,c("SPECGENENAME","MODID","MODNAME"))
    expect_equal(actual,"modspe")
    #
    expect_error(EpiTxDb:::.tableJoinSelector(),
                 'argument "tName" is missing')
    expect_error(EpiTxDb:::.tableJoinSelector("test"),
                 'No query for this combination of tables')
    actual <- EpiTxDb:::.tableJoinSelector("modrea")
    expect_type(actual,"character")
    #
    expect_error(EpiTxDb:::.makeJoinSQL(etdb),
                 'argument "cnames" is missing')
    expect_equal(EpiTxDb:::.makeJoinSQL(etdb,c("MODID","MODNAME","RXGENENAME")),
                 actual)
    #
    expect_error(EpiTxDb:::.makeSelectList(etdb),
                 'argument "cnames" is missing')
    actual <- EpiTxDb:::.makeSelectList(etdb, c("MODID","MODNAME","RXGENENAME"))
    expect_type(actual,"character")
    expect_equal(actual,"m._mod_id, m.mod_name, r.rx_genename")
    actual <- EpiTxDb:::.makeSelectList(etdb, c("MODID","MODNAME","RXGENENAME"),
                                        abbrev = FALSE)
    expect_type(actual,"character")
    expect_equal(actual,"_mod_id, mod_name, rx_genename")
    #
    expect_error(EpiTxDb:::.makeKeyList(etdb),
                 'argument "keytype" is missing')
    expect_error(EpiTxDb:::.makeKeyList(etdb, "1"),
                 'argument "keytype" is missing')
    actual <- EpiTxDb:::.makeKeyList(etdb, "1", "MODID")
    expect_type(actual,"character")
    expect_equal(actual, "m._mod_id IN ( '1' )")
    #
    expect_error(EpiTxDb:::.select(etdb, "1", "RX_GENENAME", "MODID"),
                 'Invalid columns: RX_GENENAME')
    actual <- columns(etdb)
    expect_type(actual, "character")
    expect_equal(sort(unname(EpiTxDb:::.makeColAbbreviations(etdb))),actual)
    #
    expect_message(actual <- EpiTxDb:::.select(etdb, "1", "RXGENENAME", "MODID"),
                   "'select\\(\\)' returned 1:1 mapping")
    expect_s3_class(actual, "data.frame")
    expect_equal(colnames(actual),c("MODID","RXGENENAME"))
    expect_equal(expect_message(select(etdb, "1", "RXGENENAME", "MODID"),
                                "'select\\(\\)' returned 1:1 mapping"),actual)
    #
    actual <- keys(etdb, "RXGENENAME")
    expect_type(actual,"character")
    expect_equal(actual,c("dyskerin pseudouridine synthase 1","fibrillarin"))
    actual <- keys(etdb)
    expect_type(actual,"character")
    expect_equal(actual,as.character(seq(1,235)))
    #
    actual <- keytypes(etdb)
    expect_type(actual,"character")
    dbDisconnect(etdb$conn)
})
