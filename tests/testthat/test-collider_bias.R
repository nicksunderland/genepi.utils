test_that("slopehunter", {

  # params
  ip   = 0.9
  pi0  = 0.6
  sxy1 = 1e-5

  # data
  gwas1 <- data.table::fread(system.file("extdata", "example2_gwas_sumstats.tsv", package="genepi.utils"))
  mapping<- list(SNP="MARKER", CHR="CHR", BP="POS", BETA="BETA", SE="SE", P="P", EAF="EAF", EA="A1", OA="A2")
  gwas1 <- standardise_gwas(gwas1, mapping)
  gwas2 <- standardise_gwas(gwas1, mapping)
  gwas2[, BETA := runif(.N)]

  # genepi.utils SH
  sh_geu <- slopehunter(gwas1, gwas2, ip=ip, pi0=pi0, sxy1=sxy1, bootstraps=0)

  # SlopeHunter SH
  inc <- SlopeHunter::format_data(gwas1,
                           type = "incidence",
                           snp_col = "SNP",
                           beta_col = "BETA",
                           se_col = "SE",
                           pval_col = "P",
                           eaf_col = "EAF",
                           effect_allele_col = "EA",
                           other_allele_col = "OA",
                           chr_col = "CHR",
                           pos_col = "BP")

  pro <- SlopeHunter::format_data(gwas2,
                                type = "prognosis",
                                snp_col = "SNP",
                                beta_col = "BETA",
                                se_col = "SE",
                                pval_col = "P",
                                eaf_col = "EAF",
                                effect_allele_col = "EA",
                                other_allele_col = "OA",
                                chr_col = "CHR",
                                pos_col = "BP")

  h <- SlopeHunter::harmonise_effects(inc,pro) |> as.data.frame()

  h <- h[h$remove==FALSE & h$palindromic==FALSE, ]

  sh <- SlopeHunter::hunt(h, xp_thresh=ip, init_pi=pi0, init_sigmaIP=sxy1, seed=2023,
                           xp_col="PVAL.incidence", yp_col="PVAL.prognosis", Bootstrapping = FALSE)


  # test that we get the same b value
  expect_equal(sh_geu$b, sh$b)

})


test_that("apply_collider_correction", {

  # params
  ip   = 0.9
  pi0  = 0.6
  sxy1 = 1e-5

  # data
  gwas1 <- data.table::fread(system.file("extdata", "example2_gwas_sumstats.tsv", package="genepi.utils"))
  mapping<- list(SNP="MARKER", CHR="CHR", BP="POS", BETA="BETA", SE="SE", P="P", EAF="EAF", EA="A1", OA="A2")
  gwas1 <- standardise_gwas(gwas1, mapping)
  gwas2 <- data.table::copy(gwas1)
  gwas2[, BETA := runif(.N)]

  # run SH package
  inc <- SlopeHunter::format_data(gwas1,
                                  type = "incidence",
                                  snp_col = "SNP",
                                  beta_col = "BETA",
                                  se_col = "SE",
                                  pval_col = "P",
                                  eaf_col = "EAF",
                                  effect_allele_col = "EA",
                                  other_allele_col = "OA",
                                  chr_col = "CHR",
                                  pos_col = "BP")

  pro <- SlopeHunter::format_data(gwas2,
                                  type = "prognosis",
                                  snp_col = "SNP",
                                  beta_col = "BETA",
                                  se_col = "SE",
                                  pval_col = "P",
                                  eaf_col = "EAF",
                                  effect_allele_col = "EA",
                                  other_allele_col = "OA",
                                  chr_col = "CHR",
                                  pos_col = "BP")

  h <- SlopeHunter::harmonise_effects(inc,pro) |> as.data.frame()

  h <- h[h$remove==FALSE & h$palindromic==FALSE, ]

  sh <- SlopeHunter::hunt(h, xp_thresh=ip, init_pi=pi0, init_sigmaIP=sxy1, seed=2023,
                          xp_col="PVAL.incidence", yp_col="PVAL.prognosis", Bootstrapping = FALSE)

  # apply SH package
  adjusted_sh <- SlopeHunter::SHadj(sh, h)

  # apply genepi.utils
  adjusted <- apply_collider_correction(gwas1, gwas2, b_correction_factor = sh$b, b_std_err = sh$bse, keep_palindromic = FALSE)

  # test the same
  expect_equal(adjusted$adjusted_beta, adjusted_sh$ybeta.adj)
  expect_equal(adjusted$adjusted_se,   adjusted_sh$yse.adj)
  expect_equal(adjusted$adjusted_p,    adjusted_sh$yp.adj)
})

