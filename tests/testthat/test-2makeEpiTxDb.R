context("makeEpiTxDb")
test_that("makeEpiTxDb:",{
    mod <- data.frame("mod_id" = 1L,
                      "mod_type" = "m1A",
                      "mod_name" = "m1A_1",
                      "mod_start" = 1L,
                      "mod_end" = 1L,
                      "sn_id" = 1L,
                      "sn_name" = "test")
    expect_equal(EpiTxDb:::.makeEpiTxDb_normarg_modifications(mod),mod)
    df <- mod
    df[,"mod_id"] <- 1
    expect_error(EpiTxDb:::.makeEpiTxDb_normarg_modifications(df),
                 "'modifications\\$mod_id' must be an integer")
    df <- mod
    df[,"mod_type"] <- "Z"
    expect_error(EpiTxDb:::.makeEpiTxDb_normarg_modifications(df),
                 "'modifications\\$mod_type' must be a valid")
    df <- mod
    df[,"mod_start"] <- 1
    expect_error(EpiTxDb:::.makeEpiTxDb_normarg_modifications(df),
                 "'modifications\\$mod_start' must be an integer")
    df <- mod
    df[,"mod_end"] <- 1
    expect_error(EpiTxDb:::.makeEpiTxDb_normarg_modifications(df),
                 "'modifications\\$mod_end' must be an integer")
    df <- mod
    df[,"mod_start"] <- 3L
    expect_error(EpiTxDb:::.makeEpiTxDb_normarg_modifications(df),
                 "modification starts must be <= modification ends")
    df <- mod
    df[,"sn_id"] <- 1
    expect_error(EpiTxDb:::.makeEpiTxDb_normarg_modifications(df),
                 "'modifications\\$sn_id' must be a integer")
    df <- mod
    df[,"mod_name"] <- 1
    expect_error(EpiTxDb:::.makeEpiTxDb_normarg_modifications(df),
                 "'modifications\\$mod_name' must be a character")
    df <- mod
    df[,"sn_name"] <- 1
    expect_error(EpiTxDb:::.makeEpiTxDb_normarg_modifications(df),
                 "'modifications\\$sn_name' must be a character")
    #
    rx <- data.frame(mod_id = 1L,
                     rx_genename = "test",
                     rx_rank = 1L,
                     rx_ensembl = "test",
                     rx_ensembltrans = "test",
                     rx_entrezid = "test")
    expect_equal(EpiTxDb:::.makeEpiTxDb_normarg_reactions(rx,1L),rx)
    df <- rx
    df[,"mod_id"] <- 1
    expect_error(EpiTxDb:::.makeEpiTxDb_normarg_reactions(df,1L),
                 "'reactions\\$mod_id' must be of type integer")
    df <- rx
    df[,"rx_genename"] <- 1
    expect_error(EpiTxDb:::.makeEpiTxDb_normarg_reactions(df,1L),
                 "'reactions\\$rx_genename' must be a character")
    df <- rx
    df[,"rx_rank"] <- 1
    expect_error(EpiTxDb:::.makeEpiTxDb_normarg_reactions(df,1L),
                 "'reactions\\$rx_rank' must be an integer vector")
    df <- rx
    df[,"rx_ensembl"] <- 1
    expect_error(EpiTxDb:::.makeEpiTxDb_normarg_reactions(df,1L),
                 "'reactions\\$rx_ensembl' must be a character")
    df <- rx
    df[,"rx_ensembltrans"] <- 1
    expect_error(EpiTxDb:::.makeEpiTxDb_normarg_reactions(df,1L),
                 "'reactions\\$rx_ensembltrans' must be a character")
    df <- rx
    df[,"rx_entrezid"] <- 1
    expect_error(EpiTxDb:::.makeEpiTxDb_normarg_reactions(df,1L),
                 "'reactions\\$rx_entrezid' must be a character")
    #
    spec <- data.frame(mod_id = 1L,
                       spec_type = "test",
                       spec_genename = "test",
                       spec_ensembl = "test",
                       spec_ensembltrans = "test",
                       spec_entrezid = "test")
    expect_equal(EpiTxDb:::.makeEpiTxDb_normarg_specifiers(spec,1L),spec)
    df <- spec
    df[,"mod_id"] <- 1
    expect_error(EpiTxDb:::.makeEpiTxDb_normarg_specifiers(df,1L),
                 "'specifier\\$mod_id' must be of type integer")
    df <- spec
    df[,"spec_type"] <- 1
    expect_error(EpiTxDb:::.makeEpiTxDb_normarg_specifiers(df,1L),
                 "'specifier\\$spec_type' must be a character")
    df <- spec
    df[,"spec_genename"] <- 1
    expect_error(EpiTxDb:::.makeEpiTxDb_normarg_specifiers(df,1L),
                 "'specifier\\$spec_genename' must be a character")
    df <- spec
    df[,"spec_ensembl"] <- 1
    expect_error(EpiTxDb:::.makeEpiTxDb_normarg_specifiers(df,1L),
                 "'specifier\\$spec_ensembl' must be a character")
    df <- spec
    df[,"spec_ensembltrans"] <- 1
    expect_error(EpiTxDb:::.makeEpiTxDb_normarg_specifiers(df,1L),
                 "'specifier\\$spec_ensembltrans' must be a character")
    df <- spec
    df[,"spec_entrezid"] <- 1
    expect_error(EpiTxDb:::.makeEpiTxDb_normarg_specifiers(df,1L),
                 "'specifier\\$spec_entrezid' must be a character")
    #
    ref <- data.frame(mod_id = 1L,
                      ref_type = "test",
                      ref = "test")
    expect_equal(EpiTxDb:::.makeEpiTxDb_normarg_references(ref,1L),ref)
    df <- ref
    df[,"mod_id"] <- 1
    expect_error(EpiTxDb:::.makeEpiTxDb_normarg_references(df,1L),
                 "'references\\$mod_id' must be of type integer")
    df <- ref
    df[,"ref_type"] <- 1
    expect_error(EpiTxDb:::.makeEpiTxDb_normarg_references(df,1L),
                 "'references\\$ref_type' must be a character")
    df <- ref
    df[,"ref"] <- 1
    expect_error(EpiTxDb:::.makeEpiTxDb_normarg_references(df,1L),
                 "'references\\$ref' must be a character")
    #
    expect_error(EpiTxDb:::.makeEpiTxDb_normarg_metadata(),
                 'argument "metadata" is missing')
    expect_equal(EpiTxDb:::.makeEpiTxDb_normarg_metadata(NULL),NULL)
    expect_error(EpiTxDb:::.makeEpiTxDb_normarg_metadata("1"),
                 "'metadata' must be NULL or a data.frame")
    expect_error(EpiTxDb:::.makeEpiTxDb_normarg_metadata(data.frame(a = 1,
                                                                    b = 2)),
                 "'metadata' columns must be \"name\" and \"value\"")
    actual <- EpiTxDb:::.makeEpiTxDb_normarg_metadata(data.frame(name = 1,
                                                                 value = 2))
    expect_s3_class(actual,"data.frame")
    expect_equal(colnames(actual),c("name","value"))
    #
    mod2 <- data.frame("mod_id" = 2L,
                       "mod_type" = "m1A",
                       "mod_name" = "m1A_1",
                       "mod_start" = 1L,
                       "mod_end" = 1L,
                       "sn_id" = 1L,
                       "sn_name" = "test")
    actual <- EpiTxDb:::.make_modifications_internal_mod_id(mod, FALSE)
    expect_equal(actual, 1L)
    actual <- EpiTxDb:::.make_modifications_internal_mod_id(mod2, FALSE)
    expect_equal(actual, 2L)
    actual <- EpiTxDb:::.make_modifications_internal_mod_id(mod, TRUE)
    expect_equal(actual, 1L)
    actual <- EpiTxDb:::.make_modifications_internal_mod_id(mod2, TRUE)
    expect_equal(actual, 1L)
    #
    rx2 <- rx
    rx2$rx_id <- 2L
    spec2 <- spec
    spec2$spec_id <- 2L
    ref2 <- ref
    ref2$ref_id <- 2L
    expect_equal(EpiTxDb:::.make_reactions_ids(rx2)$rx_id,1L)
    expect_equal(EpiTxDb:::.make_specifiers_ids(spec2)$spec_id,1L)
    expect_equal(EpiTxDb:::.make_references_ids(ref2)$ref_id,1L)
    expect_equal(EpiTxDb:::.shrink_df(rbind(rx2,rx2),"rx_id"),rx2)
    #
    etdb <- makeEpiTxDb(mod,rx,spec,ref)
    expect_s4_class(etdb,"EpiTxDb")
    expect_s3_class(metadata(etdb),"data.frame")
    dbDisconnect(etdb$conn)
})
