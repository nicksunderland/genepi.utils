test_that("get_gwas_mapping() works", {

  # # the default map
  # default <- get_gwas_mapping('default')
  # expect_equal(default, list(SNP="SNP",CHR="CHR",BP="BP",EA="EA", OA="OA", EAF="EAF", P="P", BETA="BETA", SE="SE", OR="OR", OR_SE="OR_SE", OR_LB="OR_LB", OR_UB="OR_UB", RSID="RSID"))
  #
  # # the gwama output map
  # gwama <- get_gwas_mapping('gwama')
  # expect_equal(gwama, list(RSID="rs_number",CHR="chromosome",BP="position",EA="reference_allele", OA="other_allele", BETA="beta", SE="se", EAF="eaf", P="p-value", BETA_95U="beta_95U",BETA_95L="beta_95L", Z="z", LOG10_P="_-log10_p-value", Q_STAT="q_statistic",I2="i2", N_STUDIES="n_studies", N="n_samples", EFFECTS="effects"))
  #
  # # not going to check every defined map
  #
  #
  # # a custom map
  # custom <- get_gwas_mapping(list(foo='foo', bar='bar'))
  # expect_equal(custom, list(foo='foo', bar='bar'))
  #
  # # invalid map arguments
  # expect_error(get_gwas_mapping("foo"))
  # expect_error(get_gwas_mapping(1))
})
