#' @title Calculate LD matrix
#' @description
#' Based on the ieugwasr function (see reference)
#' @references [ieugwasr::ld_matrix_local()](https://github.com/MRCIEU/ieugwasr/blob/33e4629f4dacd635c68e690bb5648de529c333cc/R/ld_matrix.R#L92C1-L92C16)
#' @param variants data.frame like object, or file path, with at least column `RSID`
#' @inheritParams clump
#' @return an LD matrix
#' @export
#'
ld_matrix <- function(variants,
                      plink2    = genepi.utils::which_plink2(),
                      plink_ref = genepi.utils::which_1000G_reference(build="GRCh37"),
                      logging   = TRUE) {

  # to data.table format
  variants <- import_table(variants)

  # write out the variants file
  variants_file <- tempfile()
  data.table::fwrite(variants[, "RSID"], variants_file, sep="\t", col.names=FALSE)

  # checks
  stopifnot("At least column(s) RSID must be present in the `variants`" = all(c("RSID") %in% colnames(variants)))

  # output allele file
  plink_output_alleles <- "/Users/xx20081/Downloads/alleles" # tempfile()

  # build command line command
  cmd <- paste(ifelse(is.null(plink2), "plink2", plink2),
               "--pfile", plink_ref,
               "--extract", variants_file,
               "--make-just-pvar",
               "--out", plink_output_alleles)

  # run clumping
  system(cmd)

  # the pvar file
  allele_info <- data.table::fread(paste0(plink_output_alleles,".pvar"), skip="#CHROM	POS	ID	REF	ALT", select=c("ID","#CHROM","POS","REF","ALT"))
  data.table::setnames(allele_info, names(allele_info), c("RSID","CHR","BP","REF","ALT"))

  # output files
  plink_output_ldmat <- tempfile()

  # build command line command
  cmd <- paste(ifelse(is.null(plink2), "plink2", plink2),
               "--pfile", plink_ref,
               "--extract", variants_file,
               "--r-phased square",
               "--out", plink_output_ldmat)

  # run clumping
  system(cmd)

  # the ld mat file
  ld_mat     <- data.table::fread(paste0(plink_output_ldmat,".phased.vcor1"), sep="\t", header=FALSE)
  ld_mat_ids <- data.table::fread(paste0(plink_output_ldmat,".phased.vcor1.vars"), sep="\t", header=FALSE) |> `colnames<-`("RSID")
  ld_mat     <- as.matrix(ld_mat, rownames=ld_mat_ids$RSID) |> `colnames<-`(ld_mat_ids$RSID)

  # logging info
  if(logging) {
    attr(ld_mat, "log") <- paste(readLines(paste0(plink_output_ldmat,".log")), collapse="\n")
  }

  # add allele info
  attr(ld_mat, "allele_info") <- as.data.frame(allele_info)

  # return the matrix
  return(ld_mat)
}


#' @title Re-harmonise against LD matrix reference alleles
#' @description
#' Based on the TwoSampleMR function (see reference)
#' @references [TwoSampleMR::harmonise_ld_dat()](https://github.com/MRCIEU/TwoSampleMR/blob/cbd03e6ac58a81922248eb1265cb3ee36b9a76ce/R/other_formats.R#L76)
#' @param harm a data.table, output from `harmonise()`
#' @param ld_mat a matrix, output from `ld_matrix()`
#' @inheritParams harmonise
#' @return a two-element `list(harm=harm, ld_mat=ld_mat)`
#' @export
#'
harmonise_ld_dat <- function(harm, ld_mat, gwas1_trait="", gwas2_trait="") {

  # silence R CMD
  ALT_LD = REF_LD = keep_ld = EA_exp_store = EA_out_store = NULL

  # column names
  RSID_exposure <- paste0("RSID_",gwas1_trait)
  EA_exposure   <- paste0("EA_",gwas1_trait)
  OA_exposure   <- paste0("OA_",gwas1_trait)
  BETA_exposure <- paste0("BETA_",gwas1_trait)
  EAF_exposure  <- paste0("EAF_",gwas1_trait)
  EA_outcome    <- paste0("EA_",gwas2_trait)
  OA_outcome    <- paste0("OA_",gwas2_trait)
  BETA_outcome  <- paste0("BETA_",gwas2_trait)
  EAF_outcome   <- paste0("EAF_",gwas2_trait)

  # checks
  stopifnot("RSIDs in `ld_mat` not present in harmonised data" = all(colnames(ld_mat) %in% harm$RSID_exposure))

  # LD allele data
  allele_info <- attr(ld_mat, "allele_info") |> data.table::as.data.table()
  attr(ld_mat, "allele_info") <- NULL
  attr(ld_mat, "log") <- NULL
  data.table::setnames(allele_info, names(allele_info), paste0(names(allele_info),"_LD"))

  # size of harmonised set
  n_harm <- nrow(harm)

  # trim harmonised data by variants with LD info
  harm <- harm[allele_info, on=c(RSID_exposure="RSID_LD"), nomatch=NULL]

  # report
  if(n_harm < nrow(harm)) {
    message("[-] ", n_harm-nrow(harm), " harmonised rows removed as RSIDs not found in `ld_mat` allele information")
  }

  # flag appropriate way round alleles
  harm[, keep_ld := (get("EA_exposure")==ALT_LD & get("OA_exposure")==REF_LD) | # correct
                    (get("EA_exposure")==REF_LD & get("OA_exposure")==ALT_LD)]  # incorrect, but can be flipped

  # report
  if(any(!harm$keep_ld)) {
    message("[x] no SNPs could be aligned to the LD reference panel")
    return(NULL)
  }

  # flipping exposure
  harm[, (BETA_exposure) := ifelse(get("EA_exposure") != ALT_LD, get("BETA_exposure")*-1, get("BETA_exposure"))]
  harm[, (EAF_exposure)  := ifelse(get("EA_outcome") != ALT_LD, 1-get("EAF_exposure"),    get("EAF_exposure"))]
  harm[, EA_exp_store := get("EA_exposure")]
  harm[, (EA_exposure) := ifelse(get("EA_exposure") != ALT_LD, get("OA_exposure"),  get("EA_exposure"))]
  harm[, (OA_exposure) := ifelse(EA_exp_store != ALT_LD,       get("EA_exp_store"), get("OA_exposure"))]
  harm[, EA_exp_store := NULL]

  # flipping outcome
  harm[, (BETA_outcome)  := ifelse(get("EA_outcome") != ALT_LD, get("BETA_outcome")*-1,  get("BETA_outcome"))]
  harm[, (EAF_outcome)  := ifelse(get("EA_outcome") != ALT_LD, 1-get("EAF_outcome"),     get("EAF_outcome"))]
  harm[, EA_out_store := get("EA_outcome")]
  harm[, (EA_outcome) := ifelse(get("EA_outcome") != ALT_LD, get("OA_outcome"),  get("EA_outcome"))]
  harm[, (OA_outcome) := ifelse(EA_out_store != ALT_LD,      get("EA_out_store"), get("OA_outcome"))]
  harm[, EA_out_store := NULL]

  # remove
  if(any(!harm$keep_ld)) {
    message("[-] ", sum(!harm$keep_ld), " variants removed due to harmonisation issues")
    ld_mat <- ld_mat[harm$keep_ld, harm$keep_ld]
    harm <- harm[keep_ld==TRUE, ]
  }

  # return
  return(list(harm=harm, ld_mat=ld_mat))

}
