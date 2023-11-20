# GTex data vb38

# download GTEx_Analysis_v8_eQTL.tar file from https://www.gtexportal.org/home/downloads/adult-gtex#qtl
# untar the file

# directory path
getx_v8_dir <- "/Users/xx20081/Downloads/GTEx_Analysis_v8_eQTL"

# output path
output_dir <- "/Users/xx20081/Documents/local_data/gtex_v8"

#######################################
# Process
#######################################

# list the significant variants files
gtex_sig_files <- list.files("/Users/xx20081/Downloads/GTEx_Analysis_v8_eQTL", pattern = "v8.signif_variant_gene_pairs.txt.gz", full.names = T)

# loop through the files and add bind into a large data.table
gtex_data <- lapply(gtex_sig_files, function(x) {

  cat("Processing", x, "\n")

  # read the file
  d <- data.table::fread(x)

  # split the id into individual columns
  d[, c("CHR","BP","OA","EA","build") := data.table::tstrsplit(variant_id, "_", fixed=TRUE)]
  d[, CHR := sub("chr","",CHR, ignore.case = TRUE)]
  d[, BP := as.numeric(BP)]
  d[, SNP := paste0(CHR,":",BP,"[b38]",OA,",",EA)] # v8 is build 38

  # standardise and annotate with b38 rsids
  mapping <- list(SNP="SNP",GENE_ID="gene_id",TSS_DIST="tss_distance",CHR="CHR",BP="BP",OA="OA",EA="EA",EAF="maf",BETA="slope",SE="slope_se",P="pval_nominal",MA_SAMPLES="ma_samples",MA_COUNT="ma_count")
  d <- genepi.utils::standardise_gwas(d, mapping, drop=TRUE, build="GRCh38", populate_rsid=FALSE)

  # return
  return(d)
}) |>
  `names<-`(sub(".v8.signif_variant_gene_pairs.txt.gz","",basename(gtex_sig_files))) |>
  data.table::rbindlist(idcol="tissue", fill=TRUE)

# annotate with b38 rsids
future::plan(future::multisession, workers = 12)
progressr::with_progress({
  gtex_data <- genepi.utils::chrpos_to_rsid(gtex_data, chr_col="CHR", pos_col="BP", ea_col="EA", nea_col="OA", build="b38_dbsnp156")
})

# store the rsIDs and positions then lift over to b37
gtex_data[, BP_b38 := BP]
data.table::setnames(gtex_data, c("SNP","RSID"), c("SNP_b38","RSID_b38"))
gtex_data <- genepi.utils::lift(gtex_data, from="Hg38", to="Hg19", snp_col = "SNP", chr_col = "CHR", pos_col = "BP", ea_col = "EA", oa_col = "OA")

# reannotate with b37
future::plan(future::multisession, workers = 12)
progressr::with_progress({
  gtex_data <- genepi.utils::chrpos_to_rsid(gtex_data, chr_col="CHR", pos_col="BP", ea_col="EA", nea_col="OA", build="b37_dbsnp156")
})
data.table::setnames(gtex_data, c("BP","RSID"), c("BP_b37","RSID_b37"))

# write whole file for genome wide use
data.table::fwrite(gtex_data, file.path(output_dir, "gtex_v8_full.tsv.gz"), sep="\t")

# write chr files for easier use with small genome segments
dt_by_chr_list <- split(gtex_data, by="CHR")
for(i in 1:length(dt_by_chr_list)) {

  file_name <- file.path(output_dir, paste0("gtex_v8_chr", names(dt_by_chr_list)[i], ".tsv.gz"))

  data.table::fwrite(dt_by_chr_list[[i]], file_name, sep="\t")

}

# end
