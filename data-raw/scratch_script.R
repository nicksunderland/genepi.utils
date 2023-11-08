# load the package
load_all()







gwas_in = data.table::fread("/Users/xx20081/Documents/local_data/hermes_progression/biostat_disc/raw/BIOSTAT_Discovery.allcause.gz")
set.seed(2023)
gwas_in = gwas_in[sample(1:nrow(gwas_in), size=100000), c("CHR","POS","OTHER_ALLELE","EFFECT_ALLELE","EAF","BETA","SE","P","EUR")]
gwas_in[, BETA := gwas_in[sample(nrow(gwas_in)), BETA]]
gwas_in[, P    := gwas_in[sample(nrow(gwas_in)), P]]
gwas_in[, P    := gwas_in[sample(nrow(gwas_in)), P]]
gwas_in[, SNP  := gwas_in[, paste0(CHR,":",POS,"[b37]",OTHER_ALLELE,",",EFFECT_ALLELE)]]
data.table::setnames(gwas_in, c("CHR","POS","OTHER_ALLELE","EFFECT_ALLELE","EAF","BETA","SE","P","EUR"), c("CHR","BP","OA","EA","EAF","BETA","SE","P","EUR_EAF"))
data.table::setkey(gwas_in, CHR, BP)
gwas_path = "/Users/xx20081/git/genepi.utils/inst/extdata/example3_gwas_sumstats.tsv"
data.table::fwrite(gwas_in, gwas_path, sep="\t")
gwas_in = data.table::fread(gwas_path)

highlight_snps = gwas_in[SNP=="4:32205845[b37]C,T", ][["SNP"]]
annotate_snps = gwas_in[SNP=="4:32205845[b37]C,T", ][["SNP"]]

gwas=gwas_in

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
