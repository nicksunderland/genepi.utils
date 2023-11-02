# load the package
load_all()


gwas_sumstats_4.3million_rows <- "/Users/xx20081/Documents/local_data/hermes_progression/bioshift_triumph/raw/bioshift_triumph.allcause.gz"

# some GWAS data
dt <- data.table::fread( gwas_sumstats_4.3million_rows )

# run
future::plan(future::multisession, workers = 12)
progressr::with_progress({
  result <- genepi.utils::chrpos_to_rsid(dt,chr_col="CHR",pos_col="POS",ea_col="EFFECT_ALLELE",nea_col="OTHER_ALLELE",
                                         build="b37_dbsnp156", alt_rsids=FALSE, flip="report")
})


set.seed(123)

rand_idx = sample(1:nrow(result), size=10, replace = FALSE)

result[rand_idx, c("RSID","CHR","POS","EFFECT_ALLELE","OTHER_ALLELE")]
