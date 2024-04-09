#' @title Extract variants from plink binary
#' @param snp character, an rsid
#' @param win_kb numeric, window size around snp in kb
#' @param chr character, the chromosome (use instead of snp and win_kb, not in addition)
#' @param from_bp numeric, the start base position (use instead of snp and win_kb, not in addition)
#' @param to_bp numeric, the end base position (use instead of snp and win_kb, not in addition)
#' @param plink2 character / path, the plink2 executable
#' @param pfile character / path, the plink pfile set
#'
#' @return a data.table
#' @export
#'
get_pfile_variants <- new_generic('get_pfile_variants', c('snp','win_kb', 'chr', 'from_bp', 'to_bp'), function(snp, win_kb, chr, from_bp, to_bp, plink2 = genepi.utils::which_plink2(), pfile = genepi.utils::which_1000G_reference(build="GRCh37")) { S7_dispatch() })
#' @title get_pfile_variants snp and window version
#' @noRd
method(get_pfile_variants, list(class_character,
                                class_numeric,
                                class_missing,
                                class_missing,
                                class_missing)) <- function(snp, win_kb, chr, from_bp, to_bp,
                                                            plink2 = genepi.utils::which_plink2(),
                                                            pfile  = genepi.utils::which_1000G_reference(build="GRCh37")) {

  stopifnot("`snp` must be length 1" = length(snp) == 1)
  stopifnot("`win_kb` must be length 1" = length(win_kb) == 1)

  # see if using compressed files
  if(file.exists(paste0(pfile,".pvar.zst"))) {
    compressed_plink = TRUE
  } else {
    compressed_plink = FALSE
  }

  # output allele file
  variant_file <- tempfile()

  # build command line command and run
  cmd <- paste(ifelse(is.null(plink2), "plink2", plink2),
               "--pfile", ifelse(compressed_plink, paste(pfile, "vzs"), pfile),
               "--snp", snp,
               "--window", round(win_kb),
               "--make-just-pvar",
               "--keep-allele-order",
               "--set-missing-var-ids @:#_$r_$a",
               "--out", variant_file)
  system(cmd)

  if (!file.exists(paste0(variant_file,".pvar"))) {
    warning("No variants found, or plink failed")
    return(NULL)
  }

  # read in the extracted variants file
  variants <- data.table::fread(paste0(variant_file,".pvar"), skip="#CHROM	POS	ID	REF	ALT", select=c("ID","#CHROM","POS","REF","ALT"))
  data.table::setnames(variants, names(variants), c("rsid","chr","bp","ref","alt"))

  # clean up
  unlink(paste0(variant_file,".pvar"))

  return(variants)
}
#' @title get_pfile_variants chr and bp version
#' @noRd
method(get_pfile_variants, list(class_missing,
                                class_missing,
                                class_character,
                                class_numeric,
                                class_numeric)) <- function(snp, win_kb, chr, from_bp, to_bp,
                                                            plink2 = genepi.utils::which_plink2(),
                                                            pfile  = genepi.utils::which_1000G_reference(build="GRCh37")) {

  stopifnot("`chr` must be length 1" = length(chr) == 1)
  stopifnot("`from_bp` must be length 1" = length(from_bp) == 1)
  stopifnot("`to_bp` must be length 1" = length(to_bp) == 1)

  # see if using compressed files
  if(file.exists(paste0(pfile,".pvar.zst"))) {
    compressed_plink = TRUE
  } else {
    compressed_plink = FALSE
  }

  # output allele file
  variant_file <- tempfile()

  # build command line command and run
  cmd <- paste(ifelse(is.null(plink2), "plink2", plink2),
               "--pfile", ifelse(compressed_plink, paste(pfile, "vzs"), pfile),
               "--chr", chr,
               "--from-bp", round(from_bp),
               "--to-bp", round(to_bp),
               "--make-just-pvar",
               "--keep-allele-order",
               "--set-missing-var-ids @:#_$r_$a",
               "--out", variant_file)
  system(cmd)

  if (!file.exists(paste0(variant_file,".pvar"))) {
    warning("No variants found, or plink failed")
    return(NULL)
  }

  # read in the extracted variants file
  variants <- data.table::fread(paste0(variant_file,".pvar"), skip="#CHROM	POS	ID	REF	ALT", select=c("ID","#CHROM","POS","REF","ALT"))
  data.table::setnames(variants, names(variants), c("rsid","chr","bp","ref","alt"))

  # clean up
  unlink(paste0(variant_file,".pvar"))

  return(variants)
}


#' @title Get proxies for variants from plink binary
#' @param snps character vector, rsids
#' @param stat character, the R stat to calculate, one of "r2-unphased", "r2-phased", "r-unphased", "r-phased"
#' @param win_kb numeric, the window to look in around the variants
#' @param win_r2 numeric, the lower r2 limit to include in output, (for --r-phased and --r-unphased, this means |r|≥sqrt(0.2))
#' @param win_ninter numeric, controls the maximum number of other variants allowed between variant-pairs in the report. Inf = off.
#' @param plink2 character / path, the plink2 executable
#' @param pfile character / path, the plink pfile set
#'
#' @return a data.table
#' @export
#'
get_proxies <- new_generic('get_proxies', c('snps'), function(snps,
                                                              stat       = "r2-unphased",
                                                              win_kb     = 1000,
                                                              win_r2     = 0.2,
                                                              win_ninter = Inf,
                                                              plink2     = genepi.utils::which_plink2(),
                                                              pfile      = genepi.utils::which_1000G_reference(build="GRCh37")) { S7_dispatch() })
method(get_proxies, class_character) <- function(snps,
                                                 stat       = "r2-unphased",
                                                 win_kb     = 1000,
                                                 win_r2     = 0.2,
                                                 win_ninter = Inf,
                                                 plink2     = genepi.utils::which_plink2(),
                                                 pfile      = genepi.utils::which_1000G_reference(build="GRCh37")) {

  stats <- c(UNPHASED_R2 = "r2-unphased", PHASED_R2 = "r2-phased", UNPHASED_R = "r-unphased", PHASED_R = "r-phased")
  stat <- match.arg(stat, choices = unname(stats))
  stopifnot("`win_kb` must be length 1" = length(win_kb) == 1)
  stopifnot("`win_r2` must be length 1" = length(win_r2) == 1)
  stopifnot("`win_ninter` must be length 1" = length(win_ninter) == 1)

  # see if using compressed files
  if(file.exists(paste0(pfile,".pvar.zst"))) {
    compressed_plink = TRUE
  } else {
    compressed_plink = FALSE
  }

  # output file
  proxy_file <- tempfile()

  # build command line command and run
  cmd <- paste(ifelse(is.null(plink2), "plink2", plink2),
               "--pfile", ifelse(compressed_plink, paste(pfile, "vzs"), pfile),
               paste0("--", stat), "cols=+maj,+nonmaj",
               "--ld-snps", paste0(snps, collapse = ", "),
               "--ld-window-kb", win_kb,
               "--ld-window-r2", win_r2,
               ifelse(is.finite(win_ninter), paste0("--ld-window", win_ninter), ""),
               "--out", proxy_file)
  system(cmd)

  if (!file.exists(paste0(proxy_file,".vcor"))) {
    warning("No variants found, or plink failed")
    return(NULL)
  }

  # read in the extracted variants file
  cols <- c("#CHROM_A", "POS_A", "ID_A", "MAJ_A", "NONMAJ_A", "CHROM_B", "POS_B", "ID_B", "MAJ_B", "NONMAJ_B", names(which(stats == stat)))
  proxies <- data.table::fread(paste0(proxy_file,".vcor"), skip = paste(cols, collapse = "\t"), select = cols)
  proxies <- proxies[, list(rsid = ID_A,
                            chr  = `#CHROM_A`,
                            bp   = POS_A,
                            ref  = MAJ_A,
                            alt  = NONMAJ_A,
                            proxy_rsid = ID_B,
                            proxy_chr  = CHROM_B,
                            proxy_bp   = POS_B,
                            proxy_ref  = MAJ_B,
                            proxy_alt  = NONMAJ_B,
                            rstat      = get(names(which(stats == stat))))]

  # clean up
  unlink(paste0(proxy_file,".vcor"))

  return(proxies)
}


