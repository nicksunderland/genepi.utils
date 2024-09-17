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
#' @param x a character vector of rsids or a GWAS object
#' @param stat character, the R stat to calculate, one of "r2-unphased", "r2-phased", "r-unphased", "r-phased"
#' @param win_kb numeric, the window to look in around the variants
#' @param win_r2 numeric, the lower r2 limit to include in output, (for --r-phased and --r-unphased, this means |r|≥sqrt(0.2))
#' @param win_ninter numeric, controls the maximum number of other variants allowed between variant-pairs in the report. Inf = off.
#' @param proxy_eaf numeric, the minimal effect allele frequency for proxy variants. NULL = eaf filtering off.
#' @param plink2 character / path, the plink2 executable
#' @param pfile character / path, the plink pfile set
#' @param ... other arguments (see below)
#' @return a data.table of variants and their proxies (if `x` is a `character` vector) or a `GWAS` object if
#' `x` is a `GWAS` object.
#' @export
#'
get_proxies <- new_generic('get_proxies', "x", function(x,
                                                        stat       = "r2-unphased",
                                                        win_kb     = 125,
                                                        win_r2     = 0.8,
                                                        win_ninter = Inf,
                                                        proxy_eaf  = NULL,
                                                        plink2     = genepi.utils::which_plink2(),
                                                        pfile      = genepi.utils::which_1000G_reference(build="GRCh37"),
                                                        ...) { S7_dispatch() })

#' @name get_proxies
method(get_proxies, class_character) <- function(x,
                                                 stat       = "r2-unphased",
                                                 win_kb     = 125,
                                                 win_r2     = 0.8,
                                                 win_ninter = Inf,
                                                 proxy_eaf  = NULL,
                                                 plink2     = genepi.utils::which_plink2(),
                                                 pfile      = genepi.utils::which_1000G_reference(build="GRCh37")) {

  stats <- c(UNPHASED_R2 = "r2-unphased", PHASED_R2 = "r2-phased", UNPHASED_R = "r-unphased", PHASED_R = "r-phased")
  stat <- match.arg(stat, choices = unname(stats))
  stopifnot("`win_kb` must be length 1" = length(win_kb) == 1)
  stopifnot("`win_r2` must be length 1" = length(win_r2) == 1)
  stopifnot("`win_ninter` must be length 1" = length(win_ninter) == 1)

  # empty return
  proxies <- data.table::data.table(rsid           = character(),
                                    chr            = character(),
                                    bp             = integer(),
                                    ref            = character(),
                                    alt            = character(),
                                    freq_alt       = numeric(),
                                    proxy_rsid     = character(),
                                    proxy_chr      = character(),
                                    proxy_bp       = integer(),
                                    proxy_ref      = character(),
                                    proxy_alt      = character(),
                                    proxy_freq_alt = numeric(),
                                    rstat          = numeric())

  # see if using compressed files
  if(file.exists(paste0(pfile,".pvar.zst"))) {
    compressed_plink = TRUE
  } else {
    compressed_plink = FALSE
  }

  # first see which variants exist in the pfile (as the next bit fails if absent)
  snp_file <- tempfile()
  exists_file <- tempfile()

  # write the snps
  data.table::fwrite(as.list(x), snp_file, sep="\t", col.names=FALSE)

  # build command line command and run
  cmd <- paste(ifelse(is.null(plink2), "plink2", plink2),
               "--pfile", ifelse(compressed_plink, paste(pfile, "vzs"), pfile),
               "--extract", snp_file,
               "--make-just-pvar",
               "--out", exists_file)
  system(cmd)

  # exit if no SNPs in reference
  if (!file.exists(paste0(exists_file,".pvar"))) {
    warning("No variants found in the reference file")
    return(proxies)
  }

  # read the found SNPs
  snps_exist <- data.table::fread(paste0(exists_file,".pvar"), skip="#CHROM	POS	ID	REF	ALT", select=c("ID"))

  # clean up
  unlink(paste0(exists_file,".pvar"))

  # report
  message(paste0("[i] ", nrow(snps_exist), "/", length(x), " (", sprintf("%.1f", 100*(nrow(snps_exist)/length(x))), "%) input SNPs found in plink reference file"))

  # rewrite the found snps
  snps_found <- tempfile()
  data.table::fwrite(snps_exist, snps_found, sep="\t", col.names = FALSE)

  # output file
  proxy_file <- tempfile()

  # build command line command and run
  cmd <- paste(ifelse(is.null(plink2), "plink2", plink2),
               "--pfile", ifelse(compressed_plink, paste(pfile, "vzs"), pfile),
               paste0("--", stat), "cols=+maj,+nonmaj,+freq",
               "--ld-snp-list", snps_found,
               "--ld-window-kb", win_kb,
               "--ld-window-r2", win_r2,
               ifelse(is.finite(win_ninter), paste("--ld-window", win_ninter), ""),
               "--out", proxy_file)
  system(cmd)

  if (!file.exists(paste0(proxy_file,".vcor"))) {
    stop(paste0("Plink `", paste0("--", stat), "` failed"))
  }

  # read in the extracted variants file
  cols <- c("#CHROM_A", "POS_A", "ID_A", "MAJ_A", "NONMAJ_A", "NONMAJ_FREQ_A", "CHROM_B", "POS_B", "ID_B", "MAJ_B", "NONMAJ_B", "NONMAJ_FREQ_B", names(which(stats == stat)))
  proxies <- data.table::fread(paste0(proxy_file,".vcor"), skip = paste(cols, collapse = "\t"), select = cols, encoding="Latin-1")
  proxies <- proxies[, list(rsid           = as.character(ID_A),
                            chr            = as.character(`#CHROM_A`),
                            bp             = as.integer(sub(".*?([0-9]+).*", "\\1", POS_A)), # plink sometimes write random "<FF>317417841" strings for bp
                            ref            = as.character(MAJ_A),
                            alt            = as.character(NONMAJ_A),
                            freq_alt       = as.numeric(NONMAJ_FREQ_A),
                            proxy_rsid     = as.character(ID_B),
                            proxy_chr      = as.character(CHROM_B),
                            proxy_bp       = as.integer(sub(".*?([0-9]+).*", "\\1", POS_B)),
                            proxy_ref      = as.character(MAJ_B),
                            proxy_alt      = as.character(NONMAJ_B),
                            proxy_freq_alt = as.numeric(NONMAJ_FREQ_B),
                            rstat          = as.numeric(get(names(which(stats == stat)))))]

  # clean up
  unlink(paste0(proxy_file,".vcor"))
  unlink(snps_found)

  # filter on eaf
  if (!is.null(proxy_eaf)) {
    nproxies <- nrow(proxies)
    proxies <- proxies[proxy_freq_alt>proxy_eaf]
    if (nrow(proxies)==0){
      warning(paste0("No variants found at `proxy_eaf<", proxy_eaf, "` (", nproxies, " removed due to frequency filtering)"))
    }
  }

  return(proxies)
}


