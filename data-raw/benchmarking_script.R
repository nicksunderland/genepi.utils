# load the package
load_all()
library(microbenchmark)
library(ggplot2)


gwas_sumstats_4.3million_rows <- "/Users/xx20081/Documents/local_data/hermes_progression/bioshift_triumph/raw/bioshift_triumph.allcause.gz"

# some GWAS data
dt <- data.table::fread( gwas_sumstats_4.3million_rows )

# benchmarking
mbm <- microbenchmark("chrpos_to_rsid: no alleles, no alt rsids" = {
                              future::plan(future::multisession, workers = 10)
                              progressr::with_progress({
                                genepi.utils::chrpos_to_rsid(dt,
                                                             chr_col="CHR",
                                                             pos_col="POS",
                                                             ea_col=NULL,
                                                             nea_col=NULL,
                                                             build="b37_dbsnp156",
                                                             alt_rsids=FALSE,
                                                             flip="report")
                              })
                        },
                      "chrpos_to_rsid: no alleles, with alt rsids" = {
                        future::plan(future::multisession, workers = 10)
                        progressr::with_progress({
                          genepi.utils::chrpos_to_rsid(dt,
                                                       chr_col="CHR",
                                                       pos_col="POS",
                                                       ea_col=NULL,
                                                       nea_col=NULL,
                                                       build="b37_dbsnp156",
                                                       alt_rsids=TRUE,
                                                       flip="report")
                        })
                      },
                      "chrpos_to_rsid: with alleles, no alt rsids" = {
                        future::plan(future::multisession, workers = 10)
                        progressr::with_progress({
                          genepi.utils::chrpos_to_rsid(dt,
                                                       chr_col="CHR",
                                                       pos_col="POS",
                                                       ea_col="EFFECT_ALLELE",
                                                       nea_col="OTHER_ALLELE",
                                                       build="b37_dbsnp156",
                                                       alt_rsids=FALSE,
                                                       flip="report")
                        })
                      },
                      "chrpos_to_rsid: with alleles, with alt rsids" = {
                        future::plan(future::multisession, workers = 10)
                        progressr::with_progress({
                          genepi.utils::chrpos_to_rsid(dt,
                                                       chr_col="CHR",
                                                       pos_col="POS",
                                                       ea_col="EFFECT_ALLELE",
                                                       nea_col="OTHER_ALLELE",
                                                       build="b37_dbsnp156",
                                                       alt_rsids=TRUE,
                                                       flip="report")
                        })
                      },
                      times = 10L)

# plot
autoplot(mbm)


