# load the package
load_all()


gwas_i <- data.table::fread("/Users/xx20081/Documents/local_data/hermes_incidence/clumped/hf_incidence_pheno1_eur.clumps.gz")
gwas_p <- data.table::fread("/Users/xx20081/Documents/local_data/results/progression_meta_analysis/all/meta.all.allcause_death.autosomes.gz")
data.table::setkey(gwas_i, "CHR", "BP")
data.table::setkey(gwas_p, "CHR", "BP")
gwas_i <- gwas_i[1:50000,]
# gwas_p <- gwas_p[1:50000,]


collider = ColliderBias(gwas_i, gwas_p, ip=c(0.1, 0.01,0.001))

# foo = slopehunter(gwas_i=gwas_i, gwas_p=gwas_p, bootstraps=10)
# foo1 = slopehunter(collider, ip=c(0.1, 0.01,0.001), bootstraps=10)

# foo2 = dudbridge(gwas_i=gwas_i, gwas_p=gwas_p, bootstraps=10)
# foo3 = dudbridge(collider, ip=c(0.1, 0.01,0.001), bootstraps=10)

# foo4 = ivw_mr(gwas_i=gwas_i, gwas_p=gwas_p, bootstraps=10)
# foo5 = ivw_mr(collider, ip=c(0.1, 0.01,0.001), bootstraps=10)

collider = analyse(collider, ip=c(0.1, 0.01,0.001), bootstraps = 2)


p <- plot_slopehunter_iters(collider)

p



gwas = genepi.utils::generate_random_gwas_data(100000)
highlight_snps = c("15:135277414[b37]G,T")
annotate_snps = c("15:135277414[b37]G,T")

p <- manhattan(gwas,
               highlight_snps = highlight_snps,
               highlight_win = 250,
               annotate_snps = annotate_snps,
               hit_table = TRUE,
               title = "HELLO",
               subtitle = "subdsadsda")
p


m <- miami(gwases           = list(gwas, data.table::copy(gwas)),
           highlight_snps   = list("top"=highlight_snps, "bottom"=highlight_snps),
           highlight_win    = list("top"=100, "bottom"=100),
           annotate_snps    = list("top"=annotate_snps,"bottom"=annotate_snps),
           title            = "Miami",
           subtitle         = list("top"="A","bottom"="B"),
           hit_table        = TRUE)
m
















gwas_sumstats_4.3million_rows <- "/Users/xx20081/Documents/local_data/hermes_progression/biostat_val/pre_qc/biostat_val.allcause_death.autosomes.gz"

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
