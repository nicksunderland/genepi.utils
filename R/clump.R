#' @title Clump a GWAS
#' @description
#' Clump variants in a GWAS using PLINK and an appropriate reference panel.
#' For example, the 1000 genomes phase 3 data can be downloaded from the PLINK
#' website (https://www.cog-genomics.org/plink/2.0/resources#phase3_1kg). To remove
#' duplicates you should run 'plink2 --pfile all_phase3 --rm-dup force-first --make-pgen --out all_phase3_nodup'
#' The path to the reference (without the plink extensions) should be passed as the
#' `plink_ref` argument. The path to the plink executable should be passed as the
#' `plink` argument.
#' @param gwas a data.frame like object with at least columns RSID, EA, OA, and P
#' @param p1 a numeric, the p-value threshold for inclusion as a clump
#' @param p2 a numeric, the p-value threshold for incorporation into a clump
#' @param r2 a numeric, the r2 value
#' @param kb a integer, the window for clumping
#' @param plink a string, path to the plink executable
#' @param plink_ref a string, path to the pfile genome reference
#' @return a data.table with additional columns `index` (logical, whether the variant is an index SNP) and `clump` (integer, the clump the variant belongs to)
#' @export
#'
clump <- function(gwas,
                  p1 = 1.0,
                  p2 = 1.0,
                  r2 = 0.1,
                  kb = 250,
                  plink     = genepi.utils::which_plink(),
                  plink_ref = genepi.utils::which_1000G_reference()) {

  SNP = RSID = SP2 = ID = i.clump = clump_member = NULL

  # checks
  stopifnot("`gwas` must be a data.frame like object" = inherits(gwas, "data.frame"))
  stopifnot("At least columns RSID, EA, OA, and P must be present in the `gwas`" = all(c("RSID","EA","OA","P") %in% colnames(gwas)))

  # to data.table format
  gwas <- data.table::as.data.table(gwas)

  # init files
  plink_input  <- tempfile()
  plink_output <- tempfile()

  # write the data out; plink wants A1/A2 coding
  data.table::setnames(gwas, c("EA","OA"), c("A1","A2"))
  gwas_out <- gwas[, SNP := RSID]
  data.table::fwrite(gwas_out, plink_input, sep="\t", nThread=parallel::detectCores())
  rm(gwas_out)

  # build command line command
  cmd <- paste(plink,
               "--pfile", plink_ref,
               "--clump", plink_input,
               "--out", plink_output,
               "--clump-p1", p1,
               "--clump-p2", p2,
               "--clump-r2", r2,
               "--clump-kb", kb)

  # run clumping
  system(cmd)

  # the clumped plink output
  clumped <- data.table::fread(paste0(plink_output,".clumps"), nThread=parallel::detectCores())

  # remove the temp files
  unlink(plink_input)
  unlink(plink_output)

  # add index
  clumped[, clump := seq_len(.N)]

  # pivot longer
  clumped_long <- clumped[, list(clump_member = unlist(data.table::tstrsplit(SP2, ","))), by=c("ID","clump")]

  # clean e.g. rs1763611(G) --> rs1763611
  clumped_long[, clump_member := sub("(rs[0-9]+).*", "\\1", clump_member)]

  # set the key
  data.table::setkey(clumped_long, ID)
  data.table::setkey(gwas, SNP)

  # flag the index SNPs with TRUE and the clump number
  gwas[clumped_long, c("index","clump") := list(TRUE, clump)]

  # set key to the clump member rsID to join again
  data.table::setkey(clumped_long, clump_member)

  # flag the clump members as not the index SNP and with the clump number
  gwas[clumped_long, c("index", "clump") := list(FALSE, i.clump)]

  # return
  return(gwas)
}
