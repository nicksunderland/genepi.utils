#' @title Clump a GWAS
#' @description
#' Clump variants in a GWAS using PLINK2 and an appropriate reference panel.
#' For example, the 1000 genomes phase 3 data can be downloaded from the PLINK
#' website (https://www.cog-genomics.org/plink/2.0/resources#phase3_1kg). To remove
#' duplicates you can run:
#' \cr \cr
#' plink2 \cr
#'   --pfile all_phase3 \cr
#'   --rm-dup force-first \cr
#'   --make-pgen \cr
#'   --out all_phase3_nodup \cr
#' \cr
#' The path to the reference (without the plink extensions) should be passed as the
#' `plink_ref` argument. The path to the plink2 executable should be passed as the
#' `plink2` argument.
#' @param gwas a data.frame like object with at least columns RSID, EA, OA, and P
#' @param p1 a numeric, the p-value threshold for inclusion as a clump
#' @param p2 a numeric, the p-value threshold for incorporation into a clump
#' @param r2 a numeric, the r2 value
#' @param kb a integer, the window for clumping
#' @param plink2 a string, path to the plink executable
#' @param plink_ref a string, path to the pfile genome reference
#' @param logging a logical, whether to set the plink logging information as attributes (`log`, `missing_id`, `missing_allele`) on the returned data.table
#' @return a data.table with additional columns `index` (logical, whether the variant is an index SNP) and `clump` (integer, the clump the variant belongs to)
#' @export
#'
clump <- function(gwas,
                  p1 = 1.0,
                  p2 = 1.0,
                  r2 = 0.1,
                  kb = 250,
                  plink2    = genepi.utils::which_plink2(),
                  plink_ref = genepi.utils::which_1000G_reference(build="GRCh37"),
                  logging = TRUE) {

  SNP = RSID = SP2 = ID = i.clump = clump_member = SNP_store = NULL

  # to data.table format
  gwas <- import_table(gwas)

  # checks
  stopifnot("At least columns RSID, EA, OA, and P must be present in the `gwas`" = all(c("RSID","EA","OA","P") %in% colnames(gwas)))

  # init files
  plink_input  <- tempfile()
  plink_output <- tempfile()

  # write the data out; plink wants A1/A2 coding
  data.table::setnames(gwas, c("EA","OA"), c("A1","A2"))
  gwas[, c("SNP", "SNP_store") := list(RSID, SNP)]
  data.table::fwrite(gwas, plink_input, sep="\t", na=".", quote=FALSE, nThread=parallel::detectCores()) # plink doesn't like empty strings as NA

  # revert coding now it's written out
  data.table::setnames(gwas, c("A1","A2"), c("EA","OA"))
  gwas[, c("SNP", "SNP_store") := list(SNP_store, NULL)]

  # build command line command
  cmd <- paste(ifelse(is.null(plink2), "plink2", plink2),
               "--pfile", plink_ref,
               "--clump", plink_input,
               "--out", plink_output,
               # "--clump-annotate A1,OR",
               "--clump-p1", p1,
               "--clump-p2", p2,
               "--clump-r2", r2,
               "--clump-kb", kb)

  # run clumping
  system(cmd)

  # the clumped plink output
  clumped <- data.table::fread(paste0(plink_output,".clumps"), nThread=parallel::detectCores())

  # add index
  clumped[, clump := seq_len(.N)]

  # pivot longer
  clumped_long <- clumped[, list(clump_member = unlist(data.table::tstrsplit(SP2, ","))), by=c("ID","clump")]

  # clean e.g. rs1763611(G) --> rs1763611
  clumped_long[, clump_member := sub("(rs[0-9]+).*", "\\1", clump_member)]

  # set the key
  data.table::setkey(clumped_long, ID)
  data.table::setkey(gwas, RSID)

  # flag the index SNPs with TRUE and the clump number
  gwas[clumped_long, c("index","clump") := list(TRUE, clump)]

  # set key to the clump member rsID to join again
  data.table::setkey(clumped_long, clump_member)

  # flag the clump members as not the index SNP and with the clump number
  gwas[clumped_long, c("index", "clump") := list(FALSE, i.clump)]

  # code clump as a factor
  gwas[, clump := factor(clump, levels=sort(unique(gwas$clump)))]

  # the (potential) plink clumping log files that might be produced
  # TODO: check if there are others
  log_files <- list(
    "log"            = paste0(plink_output,".log"),
    "missing_id"     = paste0(plink_output,".clumps.missing_id"),
    "missing_allele" = paste0(plink_output,".clumps.missing_allele")
  )

  # assess each log file
  for(i in seq_along(log_files)) {

    if(file.exists(log_files[[i]])) {

      # if logging flag - read and append logs
      if(logging) {
        # data.table logs
        if(names(log_files)[[i]] %in% c("missing_id", "missing_allele")) {

          log_info <- data.table::fread(log_files[[i]], sep="\t", header=FALSE, nThread=parallel::detectCores())

        # string / text logs
        } else {

          log_info <- paste(readLines(log_files[[i]]), collapse="\n")

        }

        # append info as data.table attribute
        data.table::setattr(gwas, names(log_files)[[i]], log_info)
      }

      # clean up log files
      unlink(log_files[[i]])
    }
  }

  # clean up the rest
  unlink(plink_input)
  unlink(plink_output)

  # return
  return(gwas)
}
