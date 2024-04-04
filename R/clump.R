# Silence R CMD check
globalVariables(c("snp", "rsid", "SP2", "ID", "i.clump", "clump_member", "snp_store", "A1", "A2", "P"),
                package = "genepi.utils")

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
#' @param gwas a data.frame like object with at least columns rsid, ea, oa, and p
#' @param p1 a numeric, the p-value threshold for inclusion as a clump
#' @param p2 a numeric, the p-value threshold for incorporation into a clump
#' @param r2 a numeric, the r2 value
#' @param kb a integer, the window for clumping
#' @param plink2 a string, path to the plink executable
#' @param plink_ref a string, path to the pfile genome reference
#' @param logging a logical, whether to set the plink logging information as attributes (`log`, `missing_id`, `missing_allele`) on the returned data.table
#' @param parallel_cores an integer, how many cores / threads to use
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
                  logging = TRUE,
                  parallel_cores = parallel::detectCores()) {

  # to data.table format
  gwas <- import_table(gwas)

  # checks
  stopifnot("At least columns rsid, ea, oa, and p must be present in the `gwas`" = all(c("rsid","ea","oa","p") %in% colnames(gwas)))

  # init files
  plink_input  <- tempfile()
  plink_output <- tempfile()

  # write the data out; plink wants A1/A2 coding
  data.table::setnames(gwas, c("ea","oa","p"), c("A1","A2","P"))
  gwas[, SNP := rsid]
  data.table::fwrite(gwas[, list(SNP,A1,A2,P)] , plink_input, sep="\t", na=".", quote=FALSE, nThread=parallel::detectCores()) # plink doesn't like empty strings as NA

  # revert coding now it's written out
  data.table::setnames(gwas, c("A1","A2","SNP","P"), c("ea","oa","snp","p"))

  # see if using compressed files
  if(file.exists(paste0(plink_ref,".pvar.zst"))) {
    compressed_plink = TRUE
  } else {
    compressed_plink = FALSE
  }

  # build command line command
  cmd <- paste(ifelse(is.null(plink2), "plink2", plink2),
               "--pfile", ifelse(compressed_plink, paste(plink_ref, "vzs"), plink_ref),
               "--clump", plink_input,
               "--out", plink_output,
               "--clump-p1", p1,
               "--clump-p2", p2,
               "--clump-r2", r2,
               "--clump-kb", kb,
               "--threads", parallel_cores)

  # run clumping
  system(cmd)

  # check for clump output
  clump_path <- paste0(plink_output,".clumps")

  if(file.exists(clump_path)) {

    # the clumped plink output
    clumped <- data.table::fread(clump_path, nThread=parallel::detectCores())

    # add index
    clumped[, clump := seq_len(.N)]

    # pivot longer
    clumped_long <- clumped[, list(clump_member = unlist(data.table::tstrsplit(SP2, ","))), by=c("ID","clump")]

    # clean e.g. rs1763611(G) --> rs1763611
    clumped_long[, clump_member := sub("(rs[0-9]+).*", "\\1", clump_member)]

    # set key to the clump member rsID to join the clump members
    data.table::setkey(gwas, rsid)
    data.table::setkey(clumped_long, clump_member)

    # flag the clump members as not the index SNP and with the clump number
    gwas[clumped_long, c("index", "clump") := list(FALSE, i.clump)]

    # set the key to join the next SNPs, do this after so singleton clumps arent overwritten
    data.table::setkey(clumped_long, ID)

    # flag the index SNPs with TRUE and the clump number
    gwas[clumped_long, c("index","clump") := list(TRUE, i.clump)]

    # code clump as a factor
    gwas[, clump := factor(clump, levels=sort(unique(gwas$clump)))]

  } else {

    warning("No clumps were found with the provided parameters, or plink failed")
    gwas[, c("index", "clump") := list(FALSE, NA_integer_)]

  }

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
