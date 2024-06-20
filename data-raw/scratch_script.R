library(devtools)
load_all()

# ------------- #
#     BMI       #
# ------------- #
bmi_path <- "/Users/xx20081/Desktop/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz"
bmi_map <- list(
  "rsid"=    list("name"= "rsid",    "alias"= "SNP",  "type"= "character"),
  "chr"=     list("name"= "chr",     "alias"= "CHR",   "type"= "character"),
  "bp"=      list("name"= "bp",      "alias"= "POS",    "type"= "integer"),
  "se"=      list("name"= "se",      "alias"= "SE",    "type"= "numeric"),
  "p"=       list("name"= "p",       "alias"= "P",     "type"= "numeric"),
  "n"=       list("name"= "n",       "alias"= "N",     "type"= "integer"),
  "ea"=      list("name"= "ea",      "alias"= "Tested_Allele",    "type"= "character"),
  "oa"=      list("name"= "oa",      "alias"= "Other_Allele",    "type"= "character"),
  "beta"=    list("name"= "beta",    "alias"= "BETA",  "type"= "numeric"),
  "eaf"=     list("name"= "eaf",     "alias"= "Freq_Tested_Allele",   "type"= "numeric")
)
bmi_mapping <- ColumnMap(lapply(bmi_map, function(x) do.call(genepi.utils::Column, x)))

# ------------- #
#      AF       #
# ------------- #
af_path <- "/Users/xx20081/Desktop/GCST90204201_buildGRCh37.tsv"
af_map <- list(
  "rsid"=    list("name"= "rsid",    "alias"= "variant_id",  "type"= "character"),
  "chr"=     list("name"= "chr",     "alias"= "chromosome",   "type"= "character"),
  "bp"=      list("name"= "bp",      "alias"= "base_pair_location",    "type"= "integer"),
  "se"=      list("name"= "se",      "alias"= "standard_error",    "type"= "numeric"),
  "p"=       list("name"= "p",       "alias"= "p_value",     "type"= "numeric"),
  "ea"=      list("name"= "ea",      "alias"= "effect_allele",    "type"= "character"),
  "oa"=      list("name"= "oa",      "alias"= "other_allele",    "type"= "character"),
  "beta"=    list("name"= "beta",    "alias"= "beta",  "type"= "numeric"),
  "eaf"=     list("name"= "eaf",     "alias"= "effect_allele_frequency",   "type"= "numeric")
)
af_mapping <- ColumnMap(lapply(af_map, function(x) do.call(genepi.utils::Column, x)))

# ------------- #
#      CAD      #
# ------------- #
cad_path <- "/Users/xx20081/Desktop/GCST90132315_buildGRCh37.tsv"
cad_map <- list(
  "rsid"=    list("name"= "rsid",    "alias"= "markername",  "type"= "character"),
  "chr"=     list("name"= "chr",     "alias"= "chromosome",   "type"= "character"),
  "bp"=      list("name"= "bp",      "alias"= "base_pair_location",    "type"= "integer"),
  "se"=      list("name"= "se",      "alias"= "standard_error",    "type"= "numeric"),
  "p"=       list("name"= "p",       "alias"= "p_value",     "type"= "numeric"),
  "ea"=      list("name"= "ea",      "alias"= "effect_allele",    "type"= "character"),
  "oa"=      list("name"= "oa",      "alias"= "other_allele",    "type"= "character"),
  "beta"=    list("name"= "beta",    "alias"= "beta",  "type"= "numeric"),
  "eaf"=     list("name"= "eaf",     "alias"= "effect_allele_frequency",   "type"= "numeric"),
  "n"=       list("name"= "n",       "alias"= "n",    "type"= "integer"),
  "ncase"=   list("name"= "ncase",   "alias"= "effective_cases",         "type"= "integer")
)
cad_mapping <- ColumnMap(lapply(cad_map, function(x) do.call(genepi.utils::Column, x)))

# ------------- #
#      HF       #
# ------------- #
hf_path <- "/Users/xx20081/Documents/local_data/hermes_incidence/raw/Pheno1_EUR/FORMAT-METAL_Pheno1_EUR.tsv.gz"
hf_map <- list(
  "rsid"=    list("name"= "rsid",    "alias"= "rsID",  "type"= "character"),
  "chr"=     list("name"= "chr",     "alias"= "chr",   "type"= "character"),
  "bp"=      list("name"= "bp",      "alias"= "pos_b37",    "type"= "integer"),
  "se"=      list("name"= "se",      "alias"= "se",    "type"= "numeric"),
  "p"=       list("name"= "p",       "alias"= "pval",     "type"= "numeric"),
  "n"=       list("name"= "n",       "alias"= "N_total",     "type"= "integer"),
  "ncase"=   list("name"= "ncase",   "alias"= "N_case", "type"= "integer"),
  "ea"=      list("name"= "ea",      "alias"= "A1",    "type"= "character"),
  "oa"=      list("name"= "oa",      "alias"= "A2",    "type"= "character"),
  "beta"=    list("name"= "beta",    "alias"= "A1_beta",  "type"= "numeric"),
  "eaf"=     list("name"= "eaf",     "alias"= "A1_freq",   "type"= "numeric")
)
hf_mapping <- ColumnMap(lapply(hf_map, function(x) do.call(genepi.utils::Column, x)))


# ------------- #
#      HR       #
# ------------- #
hr_path <- "/Users/xx20081/Desktop/ZhuZ_30940143_ukbb.bolt_460K_selfRepWhite.rhrmean.assoc.gz"
foo <- data.table::fread(hr_path)
hr_map <- list(
  "rsid"=    list("name"= "rsid",    "alias"= "SNP",  "type"= "character"),
  "chr"=     list("name"= "chr",     "alias"= "CHR",   "type"= "character"),
  "bp"=      list("name"= "bp",      "alias"= "BP",    "type"= "integer"),
  "se"=      list("name"= "se",      "alias"= "SE",    "type"= "numeric"),
  "p"=       list("name"= "p",       "alias"= "P",     "type"= "numeric"),
  "ea"=      list("name"= "ea",      "alias"= "A1",    "type"= "character"),
  "oa"=      list("name"= "oa",      "alias"= "A0",    "type"= "character"),
  "beta"=    list("name"= "beta",    "alias"= "BETA",  "type"= "numeric"),
  "eaf"=     list("name"= "eaf",     "alias"= "MAF",   "type"= "numeric")
)
hr_mapping <- ColumnMap(lapply(hr_map, function(x) do.call(genepi.utils::Column, x)))
gwas_hr  <- GWAS(dat = hr_path,
                  map = hr_mapping,
                  fill = TRUE,
                  drop = TRUE,
                  fill_rsid = "b37_dbsnp156",
                  missing_rsid = "fill_CHR:BP",
                  parallel_cores = 12,
                  id = "hr",
                  trait = "hr",
                  n = 458969)




# ------------- #
#     GWASs     #
# ------------- #
gwas_bmi <- GWAS(bmi_path, bmi_mapping, fill_rsid = F, missing_rsid = "fill_CHR:BP", id = "bmi", trait = "bmi")
gwas_af  <- GWAS(af_path,  af_mapping,  fill_rsid = F, missing_rsid = "fill_CHR:BP", id = "af",  trait = "af")

gwas_cad  <- GWAS(dat = cad_path,
                 map = cad_mapping,
                 fill = TRUE,
                 drop = TRUE,
                 fill_rsid = "b37_dbsnp156",
                 missing_rsid = "fill_CHR:BP",
                 parallel_cores = 12,
                 id = "cad",
                 trait = "cad")

caddt<- as.data.table(gwas_cad)

gwas_hf  <- GWAS(hf_path,  hf_mapping,  fill_rsid = F, missing_rsid = "fill_CHR:BP", id = "hf",  trait = "hf")

# ------------- #
#     MRobj     #
# ------------- #
mrobj <- MR(exposure = list(gwas_bmi, gwas_af), outcome = gwas_hf)

# saveRDS(mrobj, "/Users/xx20081/Desktop/mrobj.RDS")
#
mrobj <- readRDS("/Users/xx20081/Desktop/mrobj.RDS")
#
mr_dt_bmi <- as.data.table(mrobj, 1)
mr_dt_af <- as.data.table(mrobj, 2)

# ------------- #
#     Clump     #
# ------------- #
mrobj <- clump_mr(mrobj, p1 = 5e-8, p2 = 1, r2 = 0.001, kb = 10000,
                  plink2    = genepi.utils::which_plink2(),
                  plink_ref = "/Users/xx20081/Documents/local_data/genome_reference/ukb_reference_genome/uk10k",
                  parallel_cores = parallel::detectCores())

# ------------- #
#       MR      #
# ------------- #
res <- run_mr(mrobj, corr = FALSE, methods = c("mr_ivw"))

























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







# testing ld_matrix
load_all()
library(ggplot2)
library(data.table)


# rs910166
gwas = readRDS("/Users/xx20081/Desktop/gwas.RDS")

ld <- ld_matrix(gwas,
                method = "r",
                plink2    = genepi.utils::which_plink2(),
                plink_ref = "/Users/xx20081/Desktop/EUR_1kGv3_chr6_38996574_39075519")

r2 <- data.table(RSID_allele = rownames(ld$ld_mat), R2 = ld$ld_mat[, which(grepl("rs910166",colnames(ld$ld_mat)))])
r2[, RSID := sub("_[ACTG]+_[ACTG]+$","",RSID_allele)]

gwas[r2, R2:=i.R2, on="RSID"]



ggplot(gwas, aes(x=BP, y=nlog10P, color=abs(R2))) +
  geom_point() +
  xlim(c(38995543, 39075923)) +
  binned_scale(aesthetics = "color",
               scale_name = "stepsn",
               palette = function(x) c("darkblue", "lightblue", "green", "orange", "red"),
               breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
               limits = c(0, 1),
               show.limits = TRUE,
               guide = "colorsteps"
  )



fwrite(gwas[BP>38995543 & BP<39075923, ], "/Users/xx20081/Desktop/gwas.csv", sep=",")






library(devtools)
load_all()

map <- list(
  "rsid"=    list("name"= "rsid",    "alias"= "rsid",  "type"= "character"),
  "chr"=     list("name"= "chr",     "alias"= "chr",   "type"= "character"),
  "bp"=      list("name"= "bp",      "alias"= "bp",    "type"= "integer"),
  "se"=      list("name"= "se",      "alias"= "se",    "type"= "numeric"),
  "p"=       list("name"= "p",       "alias"= "p",     "type"= "numeric"),
  "n"=       list("name"= "n",       "alias"= "n",     "type"= "integer"),
  "ncase"=   list("name"= "ncase",   "alias"= "ncase", "type"= "integer"),
  "ea"=      list("name"= "ea",      "alias"= "ea",    "type"= "character"),
  "oa"=      list("name"= "oa",      "alias"= "oa",    "type"= "character"),
  "beta"=    list("name"= "beta",    "alias"= "beta",  "type"= "numeric"),
  "eaf"=     list("name"= "eaf",     "alias"= "eaf",   "type"= "numeric")
)
mapping <- ColumnMap(lapply(map, function(x) do.call(genepi.utils::Column, x)))

# bmi = "/Users/xx20081/Desktop/bmi_clumped.tsv"
incidence_gwas = "/Users/xx20081/Desktop/hermes_pheno1_eur_harmonised.gz"
progression_gwas = "/Users/xx20081/Desktop/meta.all.composite_2.autosomes.gz"
#"/Users/xx20081/Desktop/meta.all.allcause_death.autosomes.gz"


exposure <- genepi.utils::GWAS(incidence_gwas,
                               map          = mapping,
                               fill         = TRUE,
                               drop         = TRUE,
                               fill_rsid    = FALSE,
                               missing_rsid = "fill_CHR:BP",
                               id           = "hermes_incidence",
                               trait        = basename(incidence_gwas))

outcome  <- genepi.utils::GWAS(progression_gwas,
                               map          = mapping,
                               fill         = TRUE,
                               drop         = TRUE,
                               fill_rsid    = FALSE,
                               missing_rsid = "fill_CHR:BP",
                               id           = "hermes_progression",
                               trait        = basename(progression_gwas))

mr       <- genepi.utils::MR(exposure, outcome)

mr       <- genepi.utils::clump_mr(mr,
                                   p1 = 0.99,
                                   p2 = 1,
                                   r2 = 0.001,
                                   kb = 250,
                                   plink2    = genepi.utils::which_plink2(),
                                   plink_ref = genepi.utils::which_1000G_reference(build="GRCh37"),
                                   parallel_cores = parallel::detectCores())

t <- run_mr(mr, methods = "mr_weighted_mode")



incidence <- data.table::fread("/Users/xx20081/Desktop/hermes_pheno1_eur_harmonised.gz")
d <- data.table::fread("/Users/xx20081/Desktop/meta.all.allcause_death.autosomes.gz")
c1 <- data.table::fread("/Users/xx20081/Desktop/meta.all.composite_1.autosomes.gz")
c2 <- data.table::fread("/Users/xx20081/Desktop/meta.all.composite_2.autosomes.gz")
bmi <- data.table::fread("/Users/xx20081/Documents/local_data/giant_2018/bmi.giant-ukbb.meta-analysis.combined.23May2018.gz")
foo = list(incidence, d, c1, c2)

for (g in foo) {

  print(g$trait[[1]])
  print(max(g$n, na.rm=T))
  print(max(g$ncase, na.rm=T))
  print("---")
}

saveRDS(mr, "/Users/xx20081/Desktop/foo.RDS")

index_snp <- data.table::data.table(index_snp = mr@index_snp)
data.table::fwrite(index_snp, "/Users/xx20081/Desktop/foo.gz", sep = "\t")

mr = readRDS("/Users/xx20081/Desktop/foo.RDS")

res      <- genepi.utils::collider_bias(
  x           = mr,
  bias_method = "dudbridge",
  r2          = NA,
  p1          = NA,
  kb          = NA,
  plink2      = NA,
  plink_ref   = NA,
  ip          = 0.01,
  pi0         = 0.6,
  sxy1        = 1e-5,
  bootstraps  = 100,
  weighted    = TRUE,
  method      = "Simex",
  B           = 1000,
  seed        = 2023)






exposure    <- GWAS(incidence, map=map, fill=T, drop=T, fill_rsid = FALSE, missing_rsid = "fill_CHR:BP", id="BMI", trait="BMI")
progression <- GWAS(progression, map=map, fill=T, drop=T, fill_rsid = FALSE, missing_rsid = "fill_CHR:BP", id="prog", trait="prog")
mr          <- MR(exposure, progression)
bias        <- collider_bias(mr,r2=0.1,p1=5e-8,kb=10000, ip=0.01, pi=0.6, sxy1=1e-5, bootstraps=100,weighted=TRUE,method="Simex",B=1000,
                             bias_method = c("cwls"),
                             plink_ref = "/Users/xx20081/Documents/local_data/genome_reference/ukb_reference_genome/uk10k")





incidence   <- GWAS(incidence, map=map, fill=T, drop=T, fill_rsid = FALSE, missing_rsid = "fill_CHR:BP", id="hermes_pheno1", trait="hf_incidence")


foo1b <- mr(mr, method="mr_weighted_mode")

foo1c <- collider_bias(mr, bias_method = "mr_weighted_mode")



ivw <- mr_ivw(mr)

mr_results_to_data_table(ivw)[]

foo <- as.data.table(mr)
# rsid    chr        bp     ea      oa      bx   bxse      px        by     byse       py proxy_snp index_snp group ld_info
# <char> <char>     <int> <char>  <char>   <num>  <num>   <num>     <num>    <num>    <num>    <char>    <lgcl> <int>  <lgcl>
# 8383515:        rs9999995_A_G      4 185171608      G       A  0.0030 0.0091 0.74060  0.018037 0.040833 0.663437      <NA>      TRUE    NA   FALSE







# see if using compressed files
if(file.exists(paste0(plink_ref,".pvar.zst"))) {
  compressed_plink = TRUE
} else {
  compressed_plink = FALSE
}

# build command line command and run
cmd <- paste("/usr/local/plink/plink2",
             "--pfile", "/Users/xx20081/Documents/local_data/genome_reference/ukb_reference_genome/uk10k",
             "--snp", "rs117057532",
             "--window", "100",
             "--make-just-pvar",
             "--keep-allele-order",
             "--set-missing-var-ids @:#_$r_$a",
             "--out", "/Users/xx20081/Desktop/proxy_test.txt")
system(cmd)

foo <- get_pfile_variants(snp="rs117057532", win_kb = 100)
foo1 <- get_pfile_variants(chr="12", from_bp = min(foo$bp), to_bp = max(foo$bp))

foo <- get_proxies(snps="6:78655257_AT_A", win_r2 = 0.2, pfile = "/Users/xx20081/Documents/local_data/genome_reference/ukb_reference_genome/uk10k")




# VEP annotation








