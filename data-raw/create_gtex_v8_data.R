#######################################
# GTEx V8 files
# Nicholas Sunderland
# nicholas.sunderland@bristol.ac.uk
# November 2023
#######################################

#######################################
# Input
#######################################
# download GTEx_Analysis_v8_eQTL.tar file from https://www.gtexportal.org/home/downloads/adult-gtex#qtl
# untar the file and set the directory below.

#######################################
# Set paths
#######################################
getx_v8_dir <- "/Users/xx20081/Downloads/GTEx_Analysis_v8_eQTL" # the untarred file directory
output_dir  <- "/Users/xx20081/Documents/local_data/gtex_v8"    # where to output the data to

#######################################
# Expected outputs
#######################################
# files `gtex_v8_chr{1-22|X|Y}.tsv.gz`
# file `gtex_v8_full.tsv.gz`
#
# STRUCTURE:
# Classes ‘data.table’ and 'data.frame':	4607665 obs. of  19 variables:
# $ tissue    : chr  "Adipose_Subcutaneous" "Nerve_Tibial" "Cells_Cultured_fibroblasts" "Brain_Frontal_Cortex_BA9" ...
# $ GENE_NAME : chr  "RP11-34P13.13" "RP11-34P13.18" "WASH7P" "NOC2L" ...
# $ GENE_ID   : chr  "ENSG00000241860.6" "ENSG00000279457.4" "ENSG00000227232.5" "ENSG00000188976.10" ...
# $ TSS_DIST  : int  -159185 -178570 -12548 -942162 -983731 -111667 -611503 -156140 -156132 -610820 ...
# $ SNP_b37   : chr  "1:14677[b37]G,A" "1:16841[b37]G,T" "1:17005[b37]A,G" "1:17147[b37]G,A" ...
# $ SNP_b38   : chr  "1:14677[b38]G,A" "1:16841[b38]G,T" "1:17005[b38]A,G" "1:17147[b38]G,A" ...
# $ RSID_b37  : chr  "rs201327123" "rs62636368" "rs201833382" "rs867691030" ...
# $ RSID_b38  : chr  "rs201327123" "rs62636368" "rs201833382" "rs867691030" ...
# $ CHR       : chr  "1" "1" "1" "1" ...
# $ BP_b37    : int  14677 16841 17005 17147 17407 17556 17559 17722 17730 20254 ...
# $ BP_b38    : int  14677 16841 17005 17147 17407 17556 17559 17722 17730 20254 ...
# $ EA        : chr  "A" "T" "G" "A" ...
# $ OA        : chr  "G" "G" "A" "G" ...
# $ EAF       : num  0.0516 0.0357 0.0176 0.0194 0.0127 ...
# $ BETA      : num  0.701 -0.546 -0.868 -0.622 1.247 ...
# $ SE        : num  0.129 0.15 0.209 0.147 0.283 ...
# $ P         : num  7.92e-08 2.95e-04 3.85e-05 4.05e-05 1.71e-05 ...
# $ MA_SAMPLES: int  60 38 17 6 6 13 9 19 18 4 ...
# $ MA_COUNT  : int  60 38 17 6 6 13 9 19 18 4 ...
#
#
# file `gtex_v8_genes.tsv.gz`
#
# STRUCTURE:
# Classes ‘data.table’ and 'data.frame':	39832 obs. of  6 variables:
# $ gene_id   : chr  "ENSG00000227232.5" "ENSG00000268903.1" "ENSG00000269981.1" "ENSG00000241860.6" ...
# $ gene_chr  : chr  "1" "1" "1" "1" ...
# $ gene_name : chr  "WASH7P" "RP11-34P13.15" "RP11-34P13.16" "RP11-34P13.13" ...
# $ gene_start: int  14410 135141 137682 141474 185217 257864 366053 629062 629640 631074 ...
# $ gene_end  : int  29553 135895 137965 173862 195411 297502 501617 629433 630683 632616 ...
# $ strand    : chr  "-" "-" "-" "-" ...


#######################################
# Processing
#######################################

# list the significant variants files
gtex_sig_files <- list.files(getx_v8_dir, pattern = "v8.signif_variant_gene_pairs.txt.gz", full.names = T)

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




# list the genes files
gtex_genes_files <- list.files("/Users/xx20081/Downloads/GTEx_Analysis_v8_eQTL", pattern = "v8.egenes.txt.gz", full.names = T)

# loop through the files and add bind into a large data.table
gtex_gene_data <- lapply(gtex_genes_files, function(x) {

  cat("Processing", x, "\n")

  # read the file
  d <- data.table::fread(x, select = c("gene_id", "gene_chr","gene_name", "gene_start", "gene_end", "strand"))

  # return
  return(d)
}) |>
  `names<-`(sub(".v8.egenes.txt.gz","",basename(gtex_genes_files))) |>
  data.table::rbindlist(idcol="tissue", fill=TRUE)

# unique genes data.table
gtex_gene_data <- unique(gtex_gene_data, by="gene_id")[, tissue := NULL]

# join the gene name
gtex_data[gtex_gene_data, GENE_NAME := i.gene_name, on=c("GENE_ID"="gene_id")]

# lift genes over to b37
gtex_gene_data[, c("CHR","SNP","gene_start_b38","gene_end_b38") := list(sub("chr","",gene_chr), NA, gene_start, gene_end)]
gtex_gene_data <- genepi.utils::lift(gtex_gene_data, from="Hg38", to="Hg19", snp_col = "SNP", chr_col = "CHR", pos_col = "gene_start", ea_col = NULL, oa_col = NULL)
gtex_gene_data <- genepi.utils::lift(gtex_gene_data, from="Hg38", to="Hg19", snp_col = "SNP", chr_col = "CHR", pos_col = "gene_end", ea_col = NULL, oa_col = NULL)
data.table::setnames(gtex_gene_data, c("gene_start","gene_end"), c("gene_start_b37","gene_end_b37"))
gtex_gene_data[, c("SNP","CHR") := NULL]

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
gtex_data[, SNP_b37 := sub("\\[b38\\]","[b37]",SNP_b38)]

# set names
data.table::setcolorder(gtex_data, c("tissue","GENE_NAME","GENE_ID","TSS_DIST","SNP_b37","SNP_b38","RSID_b37","RSID_b38","CHR","BP_b37",
                                     "BP_b38","EA","OA","EAF","BETA","SE","P","MA_SAMPLES","MA_COUNT"))

# write whole file for genome wide use
data.table::fwrite(gtex_data, file.path(output_dir, "gtex_v8_full.tsv.gz"), sep="\t")

# write chr files for easier use with small genome segments
dt_by_chr_list <- split(gtex_data, by="CHR")
for(i in 1:length(dt_by_chr_list)) {

  file_name <- file.path(output_dir, paste0("gtex_v8_chr", names(dt_by_chr_list)[i], ".tsv.gz"))

  data.table::fwrite(dt_by_chr_list[[i]], file_name, sep="\t")

}

# write whole genes file
gtex_gene_data[, gene_chr := sub("chr","",gene_chr)]
data.table::fwrite(gtex_gene_data, file.path(output_dir, "gtex_v8_genes.tsv.gz"), sep="\t")


# end
