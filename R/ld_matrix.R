#' @title Calculate LD matrix
#' @description
#' Based on the ieugwasr function (see reference)
#' @references [ieugwasr::ld_matrix_local()](https://github.com/MRCIEU/ieugwasr/blob/33e4629f4dacd635c68e690bb5648de529c333cc/R/ld_matrix.R#L92C1-L92C16)
#' @param variants data.frame like object, or file path, with at least column `RSID`; if columns `EA`,`OA`,`BETA`,`EAF` are provided then the variants
#' will be return harmonised to the reference panel (effect allele, data = major allele, reference)
#' @inheritParams clump
#' @return an LD matrix if only variants provided, else if alleles provided a list(dat=harmonised data, ld_mat=ld_matrix)
#' @export
#'
ld_matrix <- function(dat,
                      plink2    = genepi.utils::which_plink2(),
                      plink_ref = genepi.utils::which_1000G_reference(build="GRCh37")) {


  dat = readRDS("/Users/xx20081/Desktop/gwas.RDS")
  plink2 = "/usr/local/plink/plink2"
  plink_ref = "/Users/xx20081/git/TargetExplorerData/inst/extdata/references/ALL_1kGv3/ALL_1kGv3"

  # checks
  stopifnot("dat must be a data.frame like object" = inherits(dat, "data.frame"))
  stopifnot("At least column(s) RSID must be present in `dat`" = all(c("RSID") %in% colnames(dat)))

  # remove variants with problematic rsid coding
  correct_rsid <- grepl("^rs[0-9]+$", trimws(dat$RSID))
  if(sum(!correct_rsid) > 0) {
    warning("[-] ", sum(!correct_rsid), " data rows removed due to incorrect RSIDs coding")
    dat <- dat[correct_rsid, ]
  }

  # extract the variants and allele information from the reference file
  alleles <- get_ref_allele_info(dat$RSID, plink2, plink_ref)

  # get the LD matrix
  ld_mat  <- get_ld_matrix_values(dat$RSID, plink2, plink_ref)

  # the extracted alleles may contain multi-allelic variants, so index from the matrix to match all the variants
  ld_mat <- ld_mat[alleles$RSID, alleles$RSID]

  # now rename with the allele appended RSIDs
  colnames(ld_mat) <- rownames(ld_mat) <- alleles$RSID_allele

  # if EA and OA columns provided, then assume we are harmonising a gwas dataset
  if(!all(c("EA","OA","BETA","EAF") %in% names(dat))) {

    message("ld_matrix(): no variant alleles provided, returning LD matrix")
    return(ld_mat)

  } else {

    message("ld_matrix(): variant alleles provided, returning harmonised data and LD matrix object")

    # join the LD allele columns - both ways
    dat[alleles, c("LD_REF","LD_ALT","LD_RSID_allele") := list(i.REF, i.ALT, i.RSID_allele), on=c("RSID"="RSID","EA"="REF","OA"="ALT")]
    dat[alleles, c("LD_REF","LD_ALT","LD_RSID_allele") := list(i.REF, i.ALT, i.RSID_allele), on=c("RSID"="RSID","OA"="REF","EA"="ALT")]

    # report & remove
    if(sum(is.na(dat$LD_REF)) > 0) {
      message("[-] ", sum(is.na(dat$LD_REF)), " data rows removed during harmonisation, as RSIDs not found in LD reference")
      print(dat[is.na(LD_REF), ])
    }
    dat <- dat[!is.na(LD_REF), ]

    # flag appropriate alleles
    dat[, keep_ld := (EA==LD_REF & OA==LD_ALT) | # correct
          (EA==LD_ALT & OA==LD_REF)]  # incorrect, but can be flipped

    # report & remove
    if(sum(!dat$keep_ld) > 0) {
      message("[-] ", sum(!dat$keep_ld), " data rows removed during harmonisation, as alleles cannot be harmonised to LD reference")
      print(dat[keep_ld==FALSE, ])
    }
    dat <- dat[keep_ld==TRUE, ]

    # report
    if(nrow(dat)==0) {
      warning("[x] no SNPs could be aligned to the LD reference panel")
      return(NULL)
    }

    # flipping
    dat[, BETA := ifelse(EA != LD_REF, BETA*-1, BETA)]
    dat[, EAF  := ifelse(EA != LD_REF, 1-EAF,   EAF )]
    dat[, EA_exp_store := EA]
    dat[, EA   := ifelse(EA           != LD_REF, OA,           EA)]
    dat[, OA   := ifelse(EA_exp_store != LD_REF, EA_exp_store, OA)]
    dat[, EA_exp_store := NULL]

    # flag appropriate alleles
    dat[, keep_ld := (EA==LD_REF & OA==LD_ALT)]

    # check
    if(any(!dat$keep_ld)) {
      warning("error harmonising LD matrix")
      return(NULL)
    }

    # order the LD matrix by the data
    ld_mat <- ld_mat[dat$LD_RSID_allele, dat$LD_RSID_allele]

    # clean up the dataset
    dat[, c("LD_REF","LD_ALT","keep_ld") := NULL]

    # return the harmonised data and LD matrix
    return(list(dat=dat, ld_mat=ld_mat))
  }

}


#' Title
#'
#' @param variants
#' @param plink2
#' @param plink_ref
#'
#' @return
#' @noRd
#'
get_ref_allele_info <- function(variants,
                                plink2    = genepi.utils::which_plink2(),
                                plink_ref = genepi.utils::which_1000G_reference(build="GRCh37")) {

  # checks
  stopifnot("variants must be a character vector" = is.character(variants))
  stopifnot("variants must be of the form `rs[0-9]+`" = all(grepl("^rs[0-9]+$", variants)))

  # write out the variants file
  variants_file <- tempfile()
  data.table::fwrite(data.table::data.table(RSID=variants), variants_file, sep="\t", col.names=FALSE)

  # output allele file
  alleles_file <- tempfile()

  # see if using compressed files
  if(file.exists(paste0(plink_ref,".pvar.zst"))) {
    compressed_plink = TRUE
  } else {
    compressed_plink = FALSE
  }

  # build command line command and run
  cmd <- paste(ifelse(is.null(plink2), "plink2", plink2),
               "--pfile", ifelse(compressed_plink, paste(plink_ref, "vzs"), plink_ref),
               "--extract", variants_file,
               "--make-just-pvar",
               "--keep-allele-order",
               "--set-missing-var-ids @:#_$r_$a",
               "--out", alleles_file)
  system(cmd)

  # read in the extracted alleles file
  alleles <- data.table::fread(paste0(alleles_file,".pvar"), skip="#CHROM	POS	ID	REF	ALT", select=c("ID","#CHROM","POS","REF","ALT"))
  data.table::setnames(alleles, names(alleles), c("RSID","CHR","BP","REF","ALT"))

  # some reference files give the ALT allele as comma separated alternate allele options - expand to extra rows
  alleles <- alleles[, list(RSID,CHR,BP,REF,ALT=unlist(strsplit(ALT, ",", fixed=TRUE))), by=seq_len(nrow(alleles))]
  alleles[, seq_len := NULL]
  alleles[, RSID_allele := paste0(RSID,"_",ALT,"_",REF)]

  # return alleles data.table
  return(alleles)
}



#' Title
#'
#' @param variants
#' @param plink2
#' @param plink_ref
#'
#' @return
#' @noRd
#'
get_ld_matrix_values <- function(variants,
                                 plink2    = genepi.utils::which_plink2(),
                                 plink_ref = genepi.utils::which_1000G_reference(build="GRCh37")) {

  # checks
  stopifnot("variants must be a character vector" = is.character(variants))
  stopifnot("variants must be of the form `rs[0-9]+`" = all(grepl("^rs[0-9]+$", variants)))

  # write out the variants file
  variants_file <- tempfile()
  data.table::fwrite(data.table::data.table(RSID=variants), variants_file, sep="\t", col.names=FALSE)

  # output allele file
  alleles_file <- tempfile()

  # see if using compressed files
  if(file.exists(paste0(plink_ref,".pvar.zst"))) {
    compressed_plink = TRUE
  } else {
    compressed_plink = FALSE
  }

  # create and run the LD matrix plink2 function
  plink_ldmat <- tempfile()
  cmd <- paste(ifelse(is.null(plink2), "plink2", plink2),
               "--pfile", ifelse(compressed_plink, paste(plink_ref, "vzs"), plink_ref),
               "--extract", variants_file,
               "--keep-allele-order ",
               "--r-phased square",
               "--out", plink_ldmat)
  system(cmd)

  # read in the LD matrix data
  ld_mat <- data.table::fread(paste0(plink_ldmat,".phased.vcor1"), sep="\t", header=FALSE) |> as.matrix()

  # read in the LD matrix variant rsids
  ld_mat_ids <- data.table::fread(paste0(plink_ldmat,".phased.vcor1.vars"), sep="\t", header=FALSE) |> `colnames<-`("RSID")

  # name the matrix
  colnames(ld_mat) <- rownames(ld_mat) <- ld_mat_ids$RSID

  # return the matrix
  return(ld_mat)
}
