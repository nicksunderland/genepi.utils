#' @title Calculate LD matrix
#' @description
#' Based on the ieugwasr function (see reference)
#' @references [ieugwasr::ld_matrix_local()](https://github.com/MRCIEU/ieugwasr/blob/33e4629f4dacd635c68e690bb5648de529c333cc/R/ld_matrix.R#L92C1-L92C16)
#' @param dat data.frame like object, or file path, with at least column `RSID`; if columns `EA`,`OA`,`BETA`,`EAF` are provided then the variants
#' will be return harmonised to the reference panel (effect allele, data = major allele, reference)
#' @param method a string, either `r` or `r2`
#' @param colmap a list, mapping to columns list(RSID=?,EA=?,OA=?,BETA=?,EAF=?) where ? can be a character vector in the case of harmonised datasets.
#' Warning - it is assumed that harmonised datasets are indeed harmonised, if not, any unharmonised variants will be
#' inappropriately removed.
#' @inheritParams clump
#' @param ukbb_ref path to a UKBB reference file
#' @return an LD matrix if only variants provided, else if alleles provided a list(dat=harmonised data, ld_mat=ld_matrix)
#' @export
#'
ld_matrix <- function(dat,
                      colmap    = NULL,
                      method    = "r",
                      plink2    = genepi.utils::which_plink2(),
                      plink_ref = genepi.utils::which_1000G_reference(build="GRCh37"),
                      ukbb_ref  = NULL) {

#
#   # testing
#   dat = readRDS("/Users/xx20081/Desktop/harm.RDS") |> data.table::as.data.table()
#   plink2 = "/usr/local/plink/plink2"
#   # plink_ref = "/Users/xx20081/git/TargetExplorerData/inst/extdata/references/ALL_1kGv3/ALL_1kGv3"
#   plink_ref <- NULL
#   colmap=list(
#     RSID = c("SNP_exp", "SNP_out"),
#     EA = c("effect_allele.exposure","effect_allele.outcome"),
#     OA = c("other_allele.exposure","other_allele.outcome"),
#     BETA = c("beta.exposure","beta.outcome"),
#     EAF = c("eaf.exposure","eaf.outcome")
#   )
#   method    = "r"
#   ukbb_ref <- "/Users/xx20081/git/TargetExplorerData/inst/extdata/references/UKB_LD/cache/UKB_LD_chr6_38000001_41000001"

  # R CMD checks
  i.REF = i.ALT = i.RSID_allele = LD_REF = keep_ld = LD_ALT = EA_exp_store = EA_exp_store = NULL

  # column mapping
  if(is.null(colmap)) {
    colmap <- list(RSID="RSID",EA="EA",OA="OA",BETA="BETA",EAF="EAF")
  } else {
    stopifnot(all(sapply(colmap, function(col) all(col %in% names(dat)))))
  }

  # checks
  method <- match.arg(method, choices=c("r","r2"))
  stopifnot("dat must be a data.frame like object" = inherits(dat, "data.frame"))
  dat <- data.table::as.data.table(dat)
  stopifnot("At least column(s) RSID must be present in `dat`" = all(colmap[["RSID"]] %in% names(dat)))

  # remove variants with problematic rsid coding
  correct_rsid <- rowSums(dat[, lapply(.SD, function(sd) !grepl("^rs[0-9]+$", trimws(sd))), .SDcols=colmap[['RSID']]]) == 0
  if(sum(!correct_rsid) > 0) {
    warning("[-] ", sum(!correct_rsid), " data rows removed due to incorrect RSIDs coding")
    dat <- dat[correct_rsid, ]
  }

  # extract the variants and allele information from the reference file
  if(!is.null(plink_ref)) {

    alleles <- get_ref_allele_info(dat[[colmap[['RSID']][[1]]]], plink2, plink_ref)
    ld_mat  <- get_ld_matrix_values(dat[[colmap[['RSID']][[1]]]], method, plink2, plink_ref)

    # the extracted alleles may contain multi-allelic variants in plink, so index from the matrix to match all the variants
    ld_mat <- ld_mat[alleles$RSID, alleles$RSID]

  } else if(!is.null(ukbb_ref)) {

    # UKBB allele and LD mats are already aligned
    ukbb_ld <- get_ukbb_ld(dat[[colmap[['RSID']][[1]]]], ukbb_ref)
    alleles <- ukbb_ld[["alleles"]]
    ld_mat  <- ukbb_ld[["ld_mat"]]

  }

  # filter out NaN
  na_rows <- which(colSums(is.nan(ld_mat))==nrow(ld_mat))
  if(length(na_rows) > 0) {
    na_alleles <- alleles$RSID[na_rows]
    alleles <- alleles[-na_rows, ]
    warning("[-] ", length(na_rows), " data rows removed due to NaN LD matrix values [", paste0(na_alleles, collapse = ", "), "]")
    ld_mat <- ld_mat[-na_rows, -na_rows]
  }

  # now rename with the allele appended RSIDs
  colnames(ld_mat) <- rownames(ld_mat) <- alleles$RSID_allele

  # if EA and OA columns provided, then assume we are harmonising a gwas dataset
  if(!all(sapply(colmap[c("EA","OA","BETA","EAF")], function(col) all(col %in% names(dat))))) {

    message("ld_matrix(): no variant alleles provided, returning LD matrix")
    return(ld_mat)

  } else {

    message("ld_matrix(): variant alleles provided, returning harmonised data and LD matrix object")

    # join the LD allele columns - both ways
    ld_cols        <- c("LD_REF","LD_ALT","LD_RSID_allele")
    join_cols      <- stats::setNames(c("RSID","REF","ALT"),c(colmap[["RSID"]][[1]],colmap[["EA"]][[1]],colmap[["OA"]][[1]]))
    join_flip_cols <- stats::setNames(c("RSID","REF","ALT"),c(colmap[["RSID"]][[1]],colmap[["OA"]][[1]],colmap[["EA"]][[1]]))

    dat[alleles, (ld_cols) := list(i.REF, i.ALT, i.RSID_allele), on=join_cols]
    dat[alleles, (ld_cols) := list(i.REF, i.ALT, i.RSID_allele), on=join_flip_cols]

    # report & remove
    if(sum(is.na(dat$LD_REF)) > 0) {
      message("[-] ", sum(is.na(dat$LD_REF)), " data rows removed during harmonisation, as RSIDs not found in LD reference")
      print(dat[is.na(LD_REF), ])
    }
    dat <- dat[!is.na(LD_REF), ]

    # flag appropriate alleles
    dat[, keep_ld := (get(colmap[['EA']][[1]])==LD_REF & get(colmap[['OA']][[1]])==LD_ALT) | # correct
                     (get(colmap[['EA']][[1]])==LD_ALT & get(colmap[['OA']][[1]])==LD_REF)]  # incorrect, but can be flipped

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
    dat[get(colmap[['EA']][[1]]) != LD_REF, colmap[['BETA']] := lapply(.SD, function(x) -1*x), .SDcols=colmap[['BETA']]]
    dat[get(colmap[['EA']][[1]]) != LD_REF, colmap[['EAF']]  := lapply(.SD, function(x)  1-x), .SDcols=colmap[['EAF']]]
    dat[                                  , EA_exp_store     := get(colmap[['EA']][[1]])]
    dat[get(colmap[['EA']][[1]]) != LD_REF, colmap[['EA']]   := lapply(.SD, function(x) get(colmap[['OA']][[1]])), .SDcols=colmap[['EA']]]
    dat[EA_exp_store != LD_REF            , colmap[['OA']]   := lapply(.SD, function(x) EA_exp_store), .SDcols=colmap[['OA']]]
    dat[                                  , EA_exp_store     := NULL]

    # flag appropriate alleles
    dat[, keep_ld := (get(colmap[['EA']][[1]])==LD_REF & get(colmap[['OA']][[1]])==LD_ALT)]

    # check
    if(any(!dat$keep_ld)) {
      warning("error harmonising LD matrix")
      return(NULL)
    }

    # order the LD matrix by the data
    ld_mat <- ld_mat[dat$LD_RSID_allele, dat$LD_RSID_allele]

    # round to 9 dp as might not be symmetric due to precision (some functions later complain)
    ld_mat <- round(ld_mat, 9)

    # clean up the dataset
    dat[, keep_ld := NULL]

    # return the harmonised data and LD matrix
    return(list(dat=dat, ld_mat=ld_mat))
  }

}


#' Title
#'
#' @param variants .
#' @param plink2 .
#' @param plink_ref .
#'
#' @return .
#' @noRd
#'
get_ref_allele_info <- function(variants,
                                plink2    = genepi.utils::which_plink2(),
                                plink_ref = genepi.utils::which_1000G_reference(build="GRCh37")) {

  # R CMD checks
  RSID_allele = RSID = REF = ALT = NULL

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
  alleles <- data.table::fread(paste0(alleles_file,".pvar"), skip="#CHROM	POS	ID	REF	ALT", select=c("ID","#CHROM","POS","REF","ALT"))#,"INFO"))
  data.table::setnames(alleles, names(alleles), c("RSID","CHR","BP","REF","ALT")) #,"INFO"))
  # alleles[, EUR_AF := as.numeric(sub(".*EUR_AF=(0\\.?[0-9]*).*", "\\1", INFO))]
  # alleles[, ALL_AF := as.numeric(sub(".*;AF=(0\\.?[0-9]*).*", "\\1", INFO))]

  # clean up
  unlink(paste0(alleles_file,".pvar"))

  # some reference files give the ALT allele as comma separated alternate allele options - expand to extra rows
  alleles <- alleles[, list(RSID,CHR,BP,REF,ALT=unlist(strsplit(ALT, ",", fixed=TRUE))), by=seq_len(nrow(alleles))]
  alleles[, seq_len := NULL]
  alleles[, RSID_allele := paste0(RSID,"_",REF,"_",ALT)]

  # return alleles data.table
  return(alleles)
}



#' Title
#'
#' @param variants .
#' @param plink2 .
#' @param plink_ref .
#'
#' @return .
#' @noRd
#'
get_ld_matrix_values <- function(variants,
                                 method    = "r",
                                 plink2    = genepi.utils::which_plink2(),
                                 plink_ref = genepi.utils::which_1000G_reference(build="GRCh37")) {

  # checks
  method <- match.arg(method, choices=c("r","r2"))
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
               ifelse(method=="r", "--r-phased square bin", "--r2-phased square bin"),
               "--out", plink_ldmat)
  system(cmd)

  # read in the LD matrix variant rsids
  ld_mat_ids <- data.table::fread(paste0(plink_ldmat,ifelse(method=="r",".phased.vcor1.bin.vars",".phased.vcor2.bin.vars")), sep="\t", header=FALSE) |> `colnames<-`("RSID")

  # read in the LD matrix data
  ld_dat <- readBin(paste0(plink_ldmat,ifelse(method=="r",".phased.vcor1.bin",".phased.vcor2.bin")), what="double", n=nrow(ld_mat_ids)^2)
  ld_mat <- matrix(ld_dat, nrow = nrow(ld_mat_ids))

  # clean up
  unlink(paste0(plink_ldmat,ifelse(method=="r",".phased.vcor1.bin",".phased.vcor2.bin")))
  unlink(paste0(plink_ldmat,ifelse(method=="r",".phased.vcor1.bin.vars",".phased.vcor2.bin.vars")))

  # name the matrix
  colnames(ld_mat) <- rownames(ld_mat) <- ld_mat_ids$RSID

  # return the matrix
  return(ld_mat)
}




# ukbb_ref_dt <- download_ukbb_ld("6",39016574,39055519,
#                                 "/Users/xx20081/git/TargetExplorerData/inst/extdata/references/UKB_LD/cache")
#
# variants      <- readRDS("/Users/xx20081/Desktop/gwas.RDS")$RSID
# allele_ld_obj <- get_ukbb_ld(variants, ukbb_ref_dt$root_file)


get_ukbb_ld <- function(variants, ukbb_ref) {

  # TODO: what if multiple reference files, need to join LD mats?

  # R CMD checks
  RSID_allele = RSID = REF = ALT = NULL

  # get the allele info
  allele_fst <- fst::fst(paste0(ukbb_ref,".fst"))
  ref_rsids  <- allele_fst[["RSID"]]

  # ensure variants coded like they are in the allele info
  variants <- sub("^(rs[0-9]+|[0-9]+:[0-9]+)_.*","\\1",variants)

  # index of variants that we have
  idx <- match(variants, ref_rsids)

  # full allele info for the variants
  alleles <- allele_fst[idx[!is.na(idx)], ] |> data.table::as.data.table()
  alleles[, RSID_allele := paste0(RSID,"_",REF,"_",ALT)]

  # full LD matrix for the variants we have
  ld_mat_fst <- fst::fst(paste0(ukbb_ref,"_ld_mat.fst"))
  ld_mat <- ld_mat_fst[idx[!is.na(idx)], idx[!is.na(idx)]] |> as.matrix()
  colnames(ld_mat) <- rownames(ld_mat) <- alleles$RSID

  # return
  return(list(alleles = alleles, ld_mat = ld_mat))
}


#' Title
#'
#' @param chr .
#' @param bp_start .
#' @param bp_end .
#' @param ukbb_ld_cache .
#'
#' @return a data.table of paths to the LD reference files
#' @export
#'
download_ukbb_ld <- function(chr, bp_start, bp_end, ukbb_ld_cache) {

  # R CMD checks
  available_alleles = allele_file = available_ld_matrix = ld_mat_file = RSID = NULL

  tryCatch({

    # round to nearest 500kb to try to reduced number of downloads
    bp_start <- round(bp_start - 5e5, digits = -5)
    bp_end   <- round(bp_end + 5e5,   digits = -5)

    # the LD files
    starts <- seq(1, 252000001, by = 1000000)
    ends   <- starts + 3000000

    # indices of the min and max files covering the region
    i_min <- max(which(starts <= bp_start))
    i_max <- max(c(i_min, min(which(ends >= bp_end))))

    # the file names and keys
    aws_keys <- paste0("chr",chr,"_",starts[i_min:i_max],"_",ends[i_min:i_max])
    req_files <- data.table::data.table("aws_key"     = aws_keys,
                                        "root_file"   = file.path(ukbb_ld_cache, paste0("UKB_LD_", aws_keys)),
                                        "allele_file" = file.path(ukbb_ld_cache, paste0("UKB_LD_", aws_keys,".fst")),
                                        "ld_mat_file" = file.path(ukbb_ld_cache, paste0("UKB_LD_", aws_keys,"_ld_mat.fst")))

    # if all files exist in cache, return those
    if(all(file.exists(c(req_files$allele_file, req_files$ld_mat_file)))) {

      return(req_files)

    }

    # connect to AWS
    Sys.setenv(
      AWS_ACCESS_KEY_ID = "AKIAVDGP42VXLAOJK5FF",
      AWS_SECRET_ACCESS_KEY = "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",
      AWS_REGION = "us-east-1"
    )
    s3 <- paws::s3()

    # download the files
    for(i in 1:nrow(req_files)) {

      # both files already exist
      if(file.exists(req_files$allele_file[[i]]) && file.exists(req_files$ld_mat_file[[i]])) {
        next
      }

      # get allele info file
      s3$download_file(Bucket   = "broad-alkesgroup-ukbb-ld",
                       Key      = paste0("UKBB_LD/",req_files$aws_key[[i]],".gz"),
                       Filename = paste0(req_files$root_file[[i]],".gz"))

      # convert allele file to .fst
      d <- data.table::fread(paste0(req_files$root_file[[i]],".gz"),
                             select    = list(character="rsid", character="chromosome", integer="position", character="allele1", character="allele2"),
                             col.names = c("RSID","CHR","BP","ALT","REF"))
      d[, RSID := sub("^(rs[0-9]+|[0-9]+:[0-9]+)_.*","\\1",RSID)]
      fst::write_fst(d, req_files$allele_file[[i]], compress = 100)
      rm(d)
      unlink(paste0(req_files$root_file[[i]],".gz"))

      # get LD matrix file
      s3$download_file(Bucket   = "broad-alkesgroup-ukbb-ld",
                       Key      = paste0("UKBB_LD/",req_files$aws_key[[i]],".npz"),
                       Filename = paste0(req_files$root_file[[i]],".npz"))

      # convert allele file to .fst
      ld_mat <- read_npz(paste0(req_files$root_file[[i]],".npz")) |> as.data.frame()
      fst::write_fst(ld_mat, req_files$ld_mat_file[[i]], compress = 100)
      rm(ld_mat)
      unlink(paste0(req_files$root_file[[i]],".npz"))

    }

    # check and return
    if(all(file.exists(c(req_files$allele_file, req_files$ld_mat_file)))) {

      return(req_files)

    } else {

      warning("Error downloading UKBB LD references")

    }

  },
  warning=function(e) {

    req_files[, available_alleles := file.exists(allele_file)]
    req_files[, available_ld_matrix := file.exists(ld_mat_file)]
    print(req_files)
    return(NULL)

  }) # end tryCatch

}


#' @title Read compressed NumPy file
#' @param path path to .npz file
#' @return the compressed python object as an R list
#' @noRd
#'
read_npz <- function(path){

  # python numpy to read numpy compressed files
  if(!reticulate::py_module_available("numpy")) {
    reticulate::py_install("numpy")
  }
  np <- reticulate::import("numpy")
  npz <- np$load(path)

  #create the matrix from the lower
  dims   <- npz[["shape"]]
  ld_matrix <- matrix(0, nrow=dims[[1]], ncol=dims[[2]])
  ld_matrix[cbind(npz[["row"]]+1, npz[["col"]]+1)] <- npz[["data"]]

  # add the upper triangle
  ld_matrix <- ld_matrix + t(ld_matrix)

  # return
  return(ld_matrix)
}

