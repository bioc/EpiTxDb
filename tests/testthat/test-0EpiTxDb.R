context("EpiTxDb class")
test_that("EpiTxDb:",{
    expect_error(EpiTxDb:::.format_modifications(),
                 'argument "modifications" is missing')
    expect_error(EpiTxDb:::.format_seqnames(),
                 'argument "seqnames" is missing')
    expect_error(EpiTxDb:::.format_reactions(),
                 'argument "reactions" is missing')
    expect_error(EpiTxDb:::.format_specifiers(),
                 'argument "specifiers" is missing')
    expect_error(EpiTxDb:::.format_references(),
                 'argument "references" is missing')
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
    expect_type(actual[,2],"double")
    actual <- EpiTxDb:::.format_seqnames(df, set.col.class = TRUE)
    expect_type(actual[,1],"integer")
    expect_type(actual[,2],"character")
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
    dbDisconnect(etdb$conn)
})

context("EpiTxDb SQL")
test_that("EpiTxDb SQL:",{
    expect_type(EpiTxDb:::build_SQL_CREATE_modification_table(),"character")
    expect_type(EpiTxDb:::build_SQL_CREATE_seqnames_table(),"character")
    expect_type(EpiTxDb:::build_SQL_CREATE_reaction_table(),"character")
    expect_type(EpiTxDb:::build_SQL_CREATE_specifier_table(),"character")
    expect_type(EpiTxDb:::build_SQL_CREATE_reference_table(),"character")
    expect_type(EpiTxDb:::build_SQL_CREATE_modification2reaction_table(),
                "character")
    expect_type(EpiTxDb:::build_SQL_CREATE_modification2specifier_table(),
                "character")
    expect_type(EpiTxDb:::build_SQL_CREATE_modification2reference_table(),
                "character")
    #
    expect_error(EpiTxDb:::EPITXDB_column2table(),
                 'argument "columns" is missing')
    expect_null(EpiTxDb:::EPITXDB_table_columns("test"))
    expect_error(EpiTxDb:::EPITXDB_column2table("mod_id"),
                 "mod_id: no such column")
    actual <- EpiTxDb:::EPITXDB_column2table("_mod_id")
    expect_equal(actual,c(`_mod_id` = "modification"))
    actual <- EpiTxDb:::EPITXDB_column2table("_sn_id")
    expect_equal(actual,c(`_sn_id` = "seqnames"))
    actual <- EpiTxDb:::EPITXDB_column2table("_rx_id")
    expect_equal(actual,c(`_rx_id` = "reaction"))
    actual <- EpiTxDb:::EPITXDB_column2table("_spec_id")
    expect_equal(actual,c(`_spec_id` = "specifier"))
    actual <- EpiTxDb:::EPITXDB_column2table("_ref_id")
    expect_equal(actual,c(`_ref_id` = "reference"))
    #
    expect_error(EpiTxDb:::EPITXDB_table2joinColumns(),
                 'argument "tables" is missing')
    expect_error(EpiTxDb:::EPITXDB_table2joinColumns("test"),
                 'test: no such table')
    actual <- EpiTxDb:::EPITXDB_table2joinColumns("modification")
    expect_equal(actual,c(modification = "_mod_id"))
    actual <- EpiTxDb:::EPITXDB_table2joinColumns("seqnames")
    expect_equal(actual,c(seqnames = "_sn_id"))
    actual <- EpiTxDb:::EPITXDB_table2joinColumns("reaction")
    expect_equal(actual,c(reaction = "_rx_id"))
    actual <- EpiTxDb:::EPITXDB_table2joinColumns("specifier")
    expect_equal(actual,c(specifier = "_spec_id"))
    actual <- EpiTxDb:::EPITXDB_table2joinColumns("reference")
    expect_equal(actual,c(reference = "_ref_id"))
    #
    expect_error(EpiTxDb:::.EPITXDB_add_table_bundle(),
                 'argument "tables" is missing')
    expect_equal(EpiTxDb:::.EPITXDB_add_table_bundle("modification"),
                 "modification")
    expect_equal(EpiTxDb:::.EPITXDB_add_table_bundle(c("modification",
                                                       "seqnames")),
                 c("modification", "seqnames"))
    expect_equal(EpiTxDb:::.EPITXDB_add_table_bundle(c("modification",
                                                       "reaction")),
                 c("modification", "modification2reaction", "reaction"))
    expect_equal(EpiTxDb:::.EPITXDB_add_table_bundle(c("modification",
                                                       "specifier")),
                 c("modification", "modification2specifier", "specifier"))
    expect_equal(EpiTxDb:::.EPITXDB_add_table_bundle(c("modification",
                                                       "reference")),
                 c("modification", "modification2reference", "reference"))
    #
    expect_error(EpiTxDb:::.EPITXDB_join_tables(),
                 'argument "tables" is missing')
    expect_error(EpiTxDb:::.EPITXDB_join_tables(c("specifier", "reference")),
                 'length\\(tables\\) == length\\(joinColumn\\)')
    expect_equal(EpiTxDb:::.EPITXDB_join_tables(c("modification", "reference"),
                                                c("_mod_id","_ref_id")),
                 c("modification",
                   "modification._mod_id=modification2reference._mod_id",
                   "modification2reference",
                   "modification2reference._ref_id=reference._ref_id",
                   "reference"))
    #
    etdb_file <- system.file("extdata", "EpiTxDb.Hs.hg38.snoRNAdb.sqlite",
                             package="EpiTxDb")
    etdb <- loadDb(etdb_file)
    actual <- EpiTxDb:::EpiTxDb_SELECT_from_modification(etdb)
    expect_s3_class(actual,"data.frame")
    expect_equal(colnames(actual),c("_mod_id","mod_type","mod_name","mod_start",
                                    "mod_end","_sn_id"))
    actual <- EpiTxDb:::EpiTxDb_SELECT_from_seqnames(etdb)
    expect_s3_class(actual,"data.frame")
    expect_equal(colnames(actual),c("_sn_id","sn_name"))
    actual <- EpiTxDb:::EpiTxDb_SELECT_from_reaction(etdb)
    expect_s3_class(actual,"data.frame")
    expect_equal(colnames(actual),c("_rx_id","rx_genename","rx_rank",
                                    "rx_ensembl","rx_ensembltrans",
                                    "rx_entrezid"))
    actual <- EpiTxDb:::EpiTxDb_SELECT_from_specifier(etdb)
    expect_s3_class(actual,"data.frame")
    expect_equal(colnames(actual),c("_spec_id","spec_type","spec_genename",
                                    "spec_ensembl","spec_ensembltrans",
                                    "spec_entrezid"))
    actual <- EpiTxDb:::EpiTxDb_SELECT_from_reference(etdb)
    expect_s3_class(actual,"data.frame")
    expect_equal(colnames(actual),c("_ref_id","ref_type","ref"))
    #
    expect_error(EpiTxDb:::EpiTxDb_SELECT_from_LEFT_JOIN(),
                 'argument "epitxdb" is missing')
    expect_error(EpiTxDb:::EpiTxDb_SELECT_from_LEFT_JOIN(etdb),
                 'argument "columns" is missing')
    table <- "modification"
    cols <- c("_mod_id","mod_name","spec_genename")
    actual <- EpiTxDb:::EpiTxDb_SELECT_from_LEFT_JOIN(etdb, table, cols)
    expect_s3_class(actual,"data.frame")
    expect_equal(colnames(actual),cols)
    table <- "specifier"
    cols <- c("rx_genename","spec_genename")
    actual <- EpiTxDb:::EpiTxDb_SELECT_from_LEFT_JOIN(etdb, table, cols)
    expect_s3_class(actual,"data.frame")
    expect_equal(colnames(actual),cols)
    dbDisconnect(etdb$conn)
})
