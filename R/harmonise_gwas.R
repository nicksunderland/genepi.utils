# Silence R CMD check
globalVariables(c(),
                package = "genepi.utils")





# need to test if index beta are flipping





#' @title Harmonise GWAS
#' @param gwas
#' @param ref
#' @param action
#' @return harmonised GWAS
#' @export
#'
harmonise_gwas <- new_generic('harmonise', c('gwas','ref'), function(gwas, ref, join = "chr:bp", action = 2, ...) { S7_dispatch() })

method(harmonise_gwas, list(GWAS,GWAS)) <- function(gwas, ref, join = "chr:bp", action = 2) {
  cli::cli_h1("genepi.utils::harmonise_gwas")
  r <- as.data.table(ref)
  g <- as.data.table(gwas)
  h <- harmonise_internal(g, r, join, action)
}

method(harmonise_gwas, list(GWAS,new_S3_class('data.table'))) <- function(gwas, ref, join = "chr:bp", action = 2, rmap = NULL) {
  cli::cli_h1("genepi.utils::harmonise_gwas")
  g <- as.data.table(gwas)
  r <- validate_ref(ref, rmap, join)
  h <- harmonise_internal(g, r, join, action)
}

method(harmonise_gwas, list(GWAS,class_character)) <- function(gwas, ref, join = "chr:bp", action = 2, rmap = NULL) {
  cli::cli_h1("genepi.utils::harmonise_gwas")
  stopifnot("ref must be a valid file path" = file.exists(ref))
  r <- data.table::fread(ref)
  r <- validate_ref(r, rmap, join)
  cli::cli_progress_done()
  g <- as.data.table(gwas)
  h <- harmonise_internal(g, r, join, action)
}

method(harmonise_gwas, list(new_S3_class('data.table'),new_S3_class('data.table'))) <- function(gwas, ref, join = "chr:bp", action = 2, gmap = NULL, rmap = NULL) {
  cli::cli_h1("genepi.utils::harmonise_gwas")
  r <- validate_ref(ref, rmap, join)
  g <- validate_gwas(gwas, gmap, join)
  h <- harmonise_internal(g, r, join, action)
}

method(harmonise_gwas, list(new_S3_class('data.table'),class_character)) <- function(gwas, ref, join = "chr:bp", action = 2, gmap = NULL, rmap = NULL) {
  cli::cli_h1("genepi.utils::harmonise_gwas")
  stopifnot("ref must be a valid file path" = file.exists(ref))
  cli::cli_progress_step("reading file")
  r <- data.table::fread(ref)
  r <- validate_ref(r, rmap, join)
  cli::cli_progress_done()
  g <- validate_gwas(gwas, gmap, join)
  h <- harmonise_internal(g, r, join, action)
}

method(harmonise_gwas, list(class_character,class_character)) <- function(gwas, ref, join = "chr:bp", action = 2, gmap = NULL, rmap = NULL) {
  cli::cli_h1("genepi.utils::harmonise_gwas")
  stopifnot("ref must be a valid file path" = file.exists(ref))
  stopifnot("gwas must be a valid file path" = file.exists(gwas))
  cli::cli_progress_step("reading files")
  r <- data.table::fread(ref)
  r <- validate_ref(r, rmap, join)
  g <- data.table::fread(gwas)
  g <- validate_gwas(g, gmap, join)
  cli::cli_progress_done()
  h <- harmonise_internal(g, r, join, action)
}

harmonise_internal <- function(x, y, join, action) {

  # standardise column types for join
  join <- match.arg(join, choices = c("chr:bp", "rsid"))
  join_key <- unlist(strsplit(join, ":"))
  cli::cli_progress_step("standardising join key column type(s)")
  fn_list <- list(
    "rsid" = as.character,
    "chr"  = as.character,
    "bp"   = as.integer
  )
  x[, (join_key) := Map(function(f, col) f(col), fn_list[join_key], .SD), .SDcols = join_key]
  y[, (join_key) := Map(function(f, col) f(col), fn_list[join_key], .SD), .SDcols = join_key]
  cli::cli_progress_done()

  # remove rows without join keys
  pre_nx <- nrow(x)
  pre_ny <- nrow(y)
  x <- x[complete.cases(x[, mget(join_key)])]
  y <- y[complete.cases(y[, mget(join_key)])]
  if (pre_nx != nrow(x)) cli::cli_alert_warning("{pre_nx - nrow(x)} variants in dataset 1 have no key data - removed")
  if (pre_ny != nrow(y)) cli::cli_alert_warning("{pre_ny - nrow(y)} variants in dataset 2 have no key data - removed")

  # set join key
  cli::cli_progress_step("setting join key: `{join}`")
  data.table::setkeyv(x, join_key)
  data.table::setkeyv(y, join_key)
  cli::cli_progress_done()

  # join; exclude all variants not in both datasets
  cli::cli_alert_info("{nrow(x)} variants in dataset 1")
  cli::cli_alert_info("{nrow(y)} variants in dataset 2")
  cli::cli_progress_step("merging datasets")
  h <- data.table::merge.data.table(x, y, by=join_key, all=FALSE, suffixes=c("","_ref"), allow.cartesian=TRUE)
  h[, paste0(join_key, "_ref") := mget(join_key)]
  other_y_cols <- names(h)[names(h) %in% names(y) & !names(h) %in% names(x)]
  data.table::setnames(h, other_y_cols, paste0(other_y_cols, "_ref"))
  h[, paste0(join_key, "_ref") := mget(join_key)]
  cli::cli_progress_done()

  # add columns required for processing
  h[, `:=`(keep = TRUE, palindromic = NA, flipped = NA, strand_flip = NA)]

  # report common snps, return if none
  if (nrow(h) == 0) {
    cli::cli_alert_danger("no variants remaining after join")
    return(h)
  } else {
    cli::cli_alert_info("{nrow(h)} common variants")
  }

  # recode indels
  cli::cli_progress_step("recoding indels")
  h <- recode_indels(h)
  cli::cli_progress_done()
  n_indel_fails <- sum(!h$keep)
  if (n_indel_fails > 0) cli::cli_alert_warning("{sum(!h$keep)} variants failing after indel recoding")

  # SNPs with alleles in 2 need to swap (e.g. A/C Vs. C/A)
  cli::cli_progress_step("harmonising alleles assuming positive strand")
  h[, to_swap := possible_swap(ea, oa, ea_ref, oa_ref)]

  # For A's alleles that need to swap, do swap
  h[to_swap==TRUE, tmp  := ea]
  h[to_swap==TRUE, ea   := oa]
  h[to_swap==TRUE, oa   := tmp]
  h[to_swap==TRUE, beta := -beta]
  h[to_swap==TRUE, eaf  := 1 - eaf]
  h[to_swap==TRUE, flipped := TRUE]

  # Palindromics
  h[, palindromic := is.palindromic(ea, oa)]

  # For 'NON-palindromic and alleles still DON'T match' (e.g. A/C Vs. T/G), do try flipping strand
  cli::cli_progress_step("flipping strand for mismatching alleles")
  h[        , OK := alleles_ok(ea, oa, ea_ref, oa_ref)]
  h[        , ii := !palindromic & !OK] #to_swap &
  h[ii==TRUE, ea := flip_alleles(ea)]
  h[ii==TRUE, oa := flip_alleles(oa)]
  h[ii==TRUE, strand_flip := TRUE]

  # If still don't match (e.g. A/C Vs. G/T, which - after flipping in the previous step - becomes A/C Vs. C/A) then try swapping
  cli::cli_progress_step("harmonising alleles, trying reverse strand")
  h[                        , OK      := alleles_ok(ea, oa, ea_ref, oa_ref)]
  h[                        , ii      := !palindromic & !OK]
  h[                        , to_swap := possible_swap(ea, oa, ea_ref, oa_ref)]
  h[ii==TRUE & to_swap==TRUE, tmp     := ea]
  h[ii==TRUE & to_swap==TRUE, ea      := oa]
  h[ii==TRUE & to_swap==TRUE, oa      := tmp]
  h[ii==TRUE & to_swap==TRUE, beta    := -beta]
  h[ii==TRUE & to_swap==TRUE, eaf     := 1 - eaf]
  h[ii==TRUE & to_swap==TRUE, flipped := TRUE]

  # Any SNPs left with un-matching alleles (e.g. A/C Vs. A/G) need to be removed
  cli::cli_progress_step("checking harmonisation")
  h[, OK   := alleles_ok(ea, oa, ea_ref, oa_ref)]
  h[, keep := keep & OK]
  cli::cli_progress_done()

  # clean up
  cli::cli_progress_step("cleaning up")
  n_all <- nrow(h)
  h <- h[keep == TRUE]
  h[, c("tmp","ii","OK","to_swap","keep") := NULL]
  cli::cli_progress_done()
  if (nrow(h) < n_all) cli::cli_alert_warning("{n_all - nrow(h)} variants failed harmonisation")
  cli::cli_alert_info("{nrow(h)} variants harmonised")

  # order columns
  add_cols <- c("palindromic","flipped","strand_flip")
  data.table::setcolorder(h, c(names(h)[!grepl("_ref$", names(h)) & !names(h) %in% add_cols], names(h)[grepl("_ref$", names(h))], add_cols))

  # finish
  cli::cli_process_done()
  return(h)
}


is.palindromic <- function(A1, A2)
{
  (A1 == "T" & A2 == "A") |
  (A1 == "A" & A2 == "T") |
  (A1 == "G" & A2 == "C") |
  (A1 == "C" & A2 == "G")
}

alleles_ok <- function(A1, A2, B1, B2)
{
  (A1==B1 & A2==B2) & !(A1==A2 | B1==B2)
}

possible_swap <- function(A1, A2, B1, B2)
{
  (A1==B2 & A2==B1) & !(A1==A2 | B1==B2)
}

flip_alleles <- function(x)
{
  x <- toupper(x)
  x <- gsub("^C$", "g", x)
  x <- gsub("^G$", "c", x)
  x <- gsub("^A$", "t", x)
  x <- gsub("^T$", "a", x)
  return(toupper(x))
}

recode_indels <- function(h) {

  # logical mask, whether an indel or not
  indel_flag <- nchar(h$ea)>1 | nchar(h$oa)>1 | h$ea=="D" | h$oa=="I" | nchar(h$ea_ref)>1 | nchar(h$oa_ref)>1 | h$ea_ref=="D" | h$oa_ref=="I"

  # if no indels nothing needed
  if(sum(indel_flag)==0) return(h)

  # indel row indices
  idx <- which(indel_flag)

  # row indices logical mask for each criteria
  EA1_is_I_nchar <- nchar(h[idx, ea])     > nchar(h[idx, oa])     # exposure EA is insertion
  EA1_is_D_nchar <- nchar(h[idx, ea])     < nchar(h[idx, oa])     # exposure EA is deletion
  EA2_is_I_nchar <- nchar(h[idx, ea_ref]) > nchar(h[idx, oa_ref]) # outcome  EA is insertion
  EA2_is_D_nchar <- nchar(h[idx, ea_ref]) < nchar(h[idx, oa_ref]) # outcome  EA is deletion
  EA1_is_I       <- h[idx, ea] == "I"     & h[idx, oa] == "D"     # exposure EA is insertion
  EA1_is_D       <- h[idx, ea] == "D"     & h[idx, oa] == "I"     # exposure EA is deletion
  EA2_is_I       <- h[idx, ea_ref] == "I" & h[idx, oa_ref] == "D" # outcome  EA is insertion
  EA2_is_D       <- h[idx, ea_ref] == "D" & h[idx, oa_ref] == "I" # outcome  EA is deletion

  # recode D/I to allele strings
  h[ idx[EA1_is_I_nchar & EA2_is_I], c("ea_ref","oa_ref") := list(ea,oa)]         # B1[i1] <- A1[i1] ; B2[i1] <- A2[i1]; correct, make outcome alleles the same as exposure (non-D/I)
  h[ idx[EA1_is_D_nchar & EA2_is_I], c("ea_ref","oa_ref") := list(oa,ea)]         # B1[i1] <- A2[i1] ; B2[i1] <- A1[i1]; wrong, flip outcome alleles and make same as exposure (non-D/I)
  h[ idx[EA1_is_I_nchar & EA2_is_D], c("ea_ref","oa_ref") := list(oa,ea)]         # B1[i1] <- A2[i1] ; B2[i1] <- A1[i1]; wrong, flip outcome alleles and make same as exposure (non-D/I)
  h[ idx[EA1_is_D_nchar & EA2_is_D], c("ea_ref","oa_ref") := list(ea,oa)]         # B1[i1] <- A1[i1] ; B2[i1] <- A2[i1]; correct, make outcome alleles the same as exposure (non-D/I)
  h[ idx[EA2_is_I_nchar & EA1_is_I], c("ea","oa")         := list(ea_ref,oa_ref)] # A1[i1] <- B1[i1] ; A2[i1] <- B2[i1]; correct, make exposure alleles the same as outcome (non-D/I)
  h[ idx[EA2_is_D_nchar & EA1_is_I], c("oa","ea")         := list(ea_ref,oa_ref)] # A2[i1] <- B1[i1] ; A1[i1] <- B2[i1]; wrong, flip exposure alleles and make same as outcome (non-D/I)
  h[ idx[EA2_is_I_nchar & EA1_is_D], c("oa","ea")         := list(ea_ref,oa_ref)] # A2[i1] <- B1[i1] ; A1[i1] <- B2[i1]; wrong, flip exposure alleles and make same as outcome (non-D/I)
  h[ idx[EA2_is_D_nchar & EA1_is_D], c("ea","oa")         := list(ea_ref,oa_ref)] # A1[i1] <- B1[i1] | A2[i1] <- B2[i1]; correct, make exposure alleles the same as outcome (non-D/I)

  # flag the ones not to keep i.e. that failed
  EA1_nchar  <- nchar(h[idx, ea])
  EA2_nchar  <- nchar(h[idx, ea])
  EA1_IorD   <- h[idx, ea]     == "D" | h[idx, ea] == "D"
  EA2_IorD   <- h[idx, ea_ref] == "D" | h[idx, ea_ref] == "D"
  EA1_eq_OA1 <- h[idx, ea]     == h[idx, oa]
  EA2_eq_OA2 <- h[idx, ea_ref] == h[idx, oa_ref]

  h[idx[EA1_nchar>1 & EA1_nchar==EA2_nchar & EA2_IorD], keep := FALSE]
  h[idx[EA2_nchar>1 & EA2_nchar==EA2_nchar & EA1_IorD], keep := FALSE]
  h[idx[EA1_eq_OA1], keep := FALSE]
  h[idx[EA2_eq_OA2], keep := FALSE]

  # return
  return(h)
}




validate_ref <- function(ref, map, join) {
  join <- match.arg(join, choices = c("chr:bp", "rsid"))
  if (join == "chr:bp") {
    if (is.null(map)) map <- c(chr='chr',bp='bp',ea='ea',oa='oa')
    stopifnot("rmap must be a named list with at least names `chr`, `bp`, `ea`, `oa`" = all(sapply(c('chr','bp','ea','oa'), function(x) x %in% names(map))))
  } else if (join == "rsid") {
    if (is.null(map)) map <- c(rsid='rsid',ea='ea',oa='oa')
    stopifnot("rmap must be a named list with at least names `rsid`, `ea`, `oa`" = all(sapply(c('rsid','ea','oa'), function(x) x %in% names(map))))
  }
  stopifnot("all `rmap` mapping columns must be present in the reference data - adjust the mapping or data" = all(sapply(map, function(x) x %in% names(ref))))
  data.table::setnames(ref, unname(map), names(map))
  return(ref)
}

validate_gwas <- function(gwas, map, join) {
  join <- match.arg(join, choices = c("chr:bp", "rsid"))
  if (join == "chr:bp") {
    if (is.null(map)) map <- c(chr='chr',bp='bp',ea='ea',oa='oa',eaf='eaf',beta='beta')
    stopifnot("gmap must be a named list with at least names `chr`, `bp`, `ea`, `oa`, `eaf`, `beta`" = all(sapply(c('chr','bp','ea','oa','eaf','beta'), function(x) x %in% names(map))))
  } else if (join == "rsid") {
    if (is.null(map)) map <- c(rsid='rsid',ea='ea',oa='oa',eaf='eaf',beta='beta')
    stopifnot("gmap must be a named list with at least names `rsid`, `ea`, `oa`, `eaf`, `beta`" = all(sapply(c('rsid','ea','oa','eaf','beta'), function(x) x %in% names(map))))
  }
  stopifnot("all `gmap` mapping columns must be present in the gwas data - adjust the mapping or data" = all(sapply(map, function(x) x %in% names(gwas))))
  data.table::setnames(gwas, unname(map), names(map))
  return(gwas)
}
