#' @title Harmonise GWAS effects
#' @description
#' Harmonise effects across two GWAS datasets. The `gwas` objects are expected to be in standard format - see `standardise_gwas()`.
#' Corrects strand for non-palindromic SNPs (flip BETAs). Drop palindromic variants. Drop non-matching variant alleles. \cr
#' This function is based on the same function in the [Slopehunter](https://github.com/Osmahmoud/SlopeHunter/blob/master/R/harmonise_effects.R)
#' and [TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR/blob/master/R/harmonise.R) packages. I have re-written it here to use
#' the standardised column names and use pure data.table syntax for increased processing speed.
#' @references [Slopehunter](https://github.com/Osmahmoud/SlopeHunter/blob/master/R/harmonise_effects.R) - https://github.com/Osmahmoud/SlopeHunter/blob/master/R/harmonise_effects.R
#' @references [TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR/blob/master/R/harmonise.R) - https://github.com/MRCIEU/TwoSampleMR/blob/master/R/harmonise.R
#' @param gwas1 a data.frame like object or file path, GWAS 1
#' @param gwas1_trait a string, suffix to add to gwas1 column names e.g. progression --> columns CHR_progression, BP_progression, ..., etc.
#' @param gwas2_trait a string, suffix to add to gwas2 column names e.g. progression --> columns CHR_progression, BP_progression, ..., etc.
#' @param merge a named character vector c(gwas1_col=gwas2_col) of the columns to join on; e.g. c("CHR1"="CHR2","BP1"="BP2"); can be NULL, in which case the gwases must be data.tables and already have keys set
#' @param gwas2 a data.frame like object or file path, GWAS 2
#' @return a data.table, the harmonised dataset
#' @export
#'
harmonise <- function(gwas1, gwas2, gwas1_trait="incidence", gwas2_trait="progression", merge=NULL) {

  CHR_1 = BP_1 = SNP_2 = SNP_1 = keep = to_swap = EA_1 = OA_1 = EA_2 = OA_2 = EA2tmp = BETA_2 = EAF_2 = palindromic = OK = ii = NULL

  # get the gwas from file or object
  gwas1 <- import_table(gwas1)
  gwas2 <- import_table(gwas2)

  # checks
  if(is.null(merge)) {
    stopifnot("If `merge` is NULL then `gwas1` and `gwas2` must be data.table and have appropriate keys set" = data.table::haskey(gwas1) & data.table::haskey(gwas2))
    gwas1_key <- data.table::key(gwas1)
    gwas2_key <- data.table::key(gwas2)
  } else {
    gwas1_key <- names(merge)
    gwas2_key <- unname(merge)
    stopifnot("`merge` column(s) not found in `gwas1`" = all(gwas1_key %in% colnames(gwas1)))
    stopifnot("`merge` column(s) not found in `gwas2`" = all(gwas2_key %in% names(gwas2)))
  }
  req_cols <- c("SNP","CHR","BP","EA","OA","EAF","BETA","P")
  stopifnot("At least columns SNP, CHR, BP, EA, OA, EAF, BETA, and P must be present in `gwas1`" = all(req_cols %in% colnames(gwas1)))
  stopifnot("At least columns SNP, CHR, BP, EA, OA, EAF, BETA, and P must be present in `gwas2`" = all(req_cols %in% colnames(gwas2)))

  # join; exclude all variants not in both datasets
  h <- data.table::merge.data.table(gwas1, gwas2, by.x=gwas1_key, by.y=gwas2_key, all=FALSE, suffixes=c("_1","_2"))

  # duplicate the join cols - one for each gwas, if the same
  if(all(gwas1_key == gwas2_key)) {
    data.table::setnames(h, (gwas1_key), paste0(gwas1_key, "_1"))
    h[, paste0(gwas2_key, "_2") := .SD, .SDcols = paste0(gwas1_key, "_1")]
  }

  # report loss of variants
  if(nrow(h)==0) {

    warning("No matching variants between gwas inputs")
    return(NULL)

  } else if(nrow(h) != nrow(gwas1) | nrow(h) != nrow(gwas2)) {

    msg <- paste0("Harmonisation: removed ", max(c(nrow(gwas1), nrow(gwas2))) - nrow(h), " variants. ",
                  "gwas1 (n=", nrow(gwas1),")", " gwas2 (n=", nrow(gwas2),")", " harmonised (n=", nrow(h),")\n")
    message(msg)

  }

  # logging information
  jinfo <- list()

  # before recoding add a keep column, the below will alter this
  h[, keep := TRUE]

  # indel recoding
  h <- recode_indels(h)
  jinfo[['indel_kept']] <- sum(h$keep)
  jinfo[['indel_removed']] <- sum(!h$keep)

  # SNPs with alleles in 2 need to swap (e.g. A/C Vs. C/A)
  h[, to_swap := possible_swap(EA_1, OA_1, EA_2, OA_2)]
  jinfo[['switched_alleles']] <- sum(h$to_swap)

  # For B's alleles that need to swap, do swap
  h[to_swap==TRUE, EA2tmp := EA_2]
  h[to_swap==TRUE, EA_2   := OA_2]
  h[to_swap==TRUE, OA_2   := EA2tmp]
  h[to_swap==TRUE, BETA_2 := BETA_2 * -1]
  h[to_swap==TRUE, EAF_2  := 1 - EAF_2]
  h[             , EA2tmp := NULL]

  # Palindromics
  h[, palindromic := is.palindromic(EA_1, OA_1)]
  jinfo[['palindromic']] <- sum(h$palindromic)

  # For 'NON-palindromic and alleles still DON'T match' (e.g. A/C Vs. T/G), do try flipping
  h[        , OK   := alleles_ok(EA_1, OA_1, EA_2, OA_2)]
  h[        , ii   := to_swap & !palindromic & !OK]
  h[ii==TRUE, EA_2 := flip_alleles(EA_2)]
  h[ii==TRUE, OA_2 := flip_alleles(OA_2)]
  jinfo[['flipped_alleles']] <- sum(h$ii)

  # If still don't match (e.g. A/C Vs. G/T, which - after flipping in the previous step - becomes A/C Vs. C/A) then try swapping
  h[                        , OK      := alleles_ok(EA_1, OA_1, EA_2, OA_2)]
  h[                        , ii      := !palindromic & !OK]
  h[                        , to_swap := possible_swap(EA_1, OA_1, EA_2, OA_2)]
  h[ii==TRUE & to_swap==TRUE, EA2tmp  := EA_2]
  h[ii==TRUE & to_swap==TRUE, EA_2    := OA_2]
  h[ii==TRUE & to_swap==TRUE, OA_2    := EA2tmp]
  h[ii==TRUE & to_swap==TRUE, BETA_2  := BETA_2 * -1]
  h[ii==TRUE & to_swap==TRUE, EAF_2   := 1 - EAF_2]
  jinfo[['flipped_then_switched_alleles']] <- sum(h$ii & h$to_swap)

  # Any SNPs left with un-matching alleles (e.g. A/C Vs. A/G) need to be removed
  h[, OK   := alleles_ok(EA_1, OA_1, EA_2, OA_2)]
  h[, keep := keep & OK & !palindromic]

  # recode SNP/chr:pos:allele ID string if this contains allele info i.e. capture things like:
  # 12:432412[b37]A,T
  # 12:432412:A:T
  # 12_432412_A_T
  # TODO: a function that will determine the input pattern and return the same pattern with updated alleles
  snp_id_regex <- "([0-9]+|X|Y)[:_][0-9]+([:_]|\\[(b38|b37)\\])([ACTG]+|[DI])[,_:]([ACTG]+|[DI])"
  if(any(grepl(snp_id_regex, h[, SNP_1]))) {

    # TODO: h[to_swap==TRUE, SNP := pattern_recognition_function(SNP)]
    h[, SNP_1 := paste0(CHR_1,":",BP_1,"[b37]",OA_1,",",EA_1)]
    h[, SNP_2 := SNP_1]

  }

  # clean up
  h[, EA2tmp  := NULL]
  h[, ii      := NULL]
  h[, OK      := NULL]
  h[, to_swap := NULL]

  # set output names and order
  corder <- c(t(outer(c("SNP","CHR","BP","EA","OA","EAF","BETA","P"), c("_1","_2"), paste0)), "palindromic","keep")
  data.table::setcolorder(h, corder)
  if(gwas1_trait!=""){
    gwas1_trait <- paste0("_",gwas1_trait)
  }
  if(gwas2_trait!=""){
    gwas2_trait <- paste0("_",gwas2_trait)
  }
  data.table::setnames(h, names(h), sub("_1$", gwas1_trait, names(h)))
  data.table::setnames(h, names(h), sub("_2$", gwas2_trait, names(h)))

  # set logging attribute
  data.table::setattr(h, "info", jinfo)

  # return
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
  (A1==B1 & A2==B2) &
 !(A1==A2 | B1==B2)
}

possible_swap <- function(A1, A2, B1, B2)
{
  (A1==B2 & A2==B1) &
 !(A1==A2 | B1==B2)
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

recode_indels <- function(h)
{
  EA_1 = OA_1 = EA_2 = OA_2 = keep = NULL

  # logical mask, whether an indel or not
  indel_flag <- nchar(h[,EA_1]) > 1 | nchar(h[,OA_1]) > 1 | h[,EA_1]=="D" | h[,EA_1]=="I"

  # if no indels nothing needed
  if(sum(indel_flag)==0) return(h)

  # indel row indices
  idx <- which(indel_flag)

  # row indices logical mask for each criteria
  EA1_is_I_nchar <- nchar(h[idx, EA_1]) > nchar(h[idx, OA_1]) # exposure EA is insertion
  EA1_is_D_nchar <- nchar(h[idx, EA_1]) < nchar(h[idx, OA_1]) # exposure EA is deletion
  EA2_is_I_nchar <- nchar(h[idx, EA_2]) > nchar(h[idx, OA_2]) # outcome  EA is insertion
  EA2_is_D_nchar <- nchar(h[idx, EA_2]) < nchar(h[idx, OA_2]) # outcome  EA is deletion
  EA1_is_I       <- h[idx, EA_1] == "I" & h[idx, OA_1] == "D" # exposure EA is insertion
  EA1_is_D       <- h[idx, EA_1] == "D" & h[idx, OA_1] == "I" # exposure EA is deletion
  EA2_is_I       <- h[idx, EA_2] == "I" & h[idx, OA_2] == "D" # outcome  EA is insertion
  EA2_is_D       <- h[idx, EA_2] == "D" & h[idx, OA_2] == "I" # outcome  EA is deletion


  # recode D/I to allele strings
  h[ idx[EA1_is_I_nchar & EA2_is_I], c("EA_2","OA_2") := list(EA_1,OA_1)] # B1[i1] <- A1[i1] ; B2[i1] <- A2[i1]; correct, make outcome alleles the same as exposure (non-D/I)
  h[ idx[EA1_is_D_nchar & EA2_is_I], c("EA_2","OA_2") := list(OA_1,EA_1)] # B1[i1] <- A2[i1] ; B2[i1] <- A1[i1]; wrong, flip outcome alleles and make same as exposure (non-D/I)
  h[ idx[EA1_is_I_nchar & EA2_is_D], c("EA_2","OA_2") := list(OA_1,EA_1)] # B1[i1] <- A2[i1] ; B2[i1] <- A1[i1]; wrong, flip outcome alleles and make same as exposure (non-D/I)
  h[ idx[EA1_is_D_nchar & EA2_is_D], c("EA_2","OA_2") := list(EA_1,OA_1)] # B1[i1] <- A1[i1] ; B2[i1] <- A2[i1]; correct, make outcome alleles the same as exposure (non-D/I)
  h[ idx[EA2_is_I_nchar & EA1_is_I], c("EA_1","OA_1") := list(EA_2,OA_2)] # A1[i1] <- B1[i1] ; A2[i1] <- B2[i1]; correct, make exposure alleles the same as outcome (non-D/I)
  h[ idx[EA2_is_D_nchar & EA1_is_I], c("OA_1","EA_1") := list(EA_2,OA_2)] # A2[i1] <- B1[i1] ; A1[i1] <- B2[i1]; wrong, flip exposure alleles and make same as outcome (non-D/I)
  h[ idx[EA2_is_I_nchar & EA1_is_D], c("OA_1","EA_1") := list(EA_2,OA_2)] # A2[i1] <- B1[i1] ; A1[i1] <- B2[i1]; wrong, flip exposure alleles and make same as outcome (non-D/I)
  h[ idx[EA2_is_D_nchar & EA1_is_D], c("EA_1","OA_1") := list(EA_2,OA_2)] # A1[i1] <- B1[i1] | A2[i1] <- B2[i1]; correct, make exposure alleles the same as outcome (non-D/I)

  # flag the ones not to keep i.e. that failed
  EA1_nchar  <- nchar(h[idx, EA_1])
  EA2_nchar  <- nchar(h[idx, EA_2])
  EA1_IorD   <- h[idx, EA_1] == "D" | h[idx, EA_1] == "D"
  EA2_IorD   <- h[idx, EA_2] == "D" | h[idx, EA_2] == "D"
  EA1_eq_OA1 <- h[idx, EA_1] == h[idx, OA_1]
  EA2_eq_OA2 <- h[idx, EA_2] == h[idx, OA_2]

  h[idx[EA1_nchar>1 & EA1_nchar==EA2_nchar & EA2_IorD], keep := FALSE]
  h[idx[EA2_nchar>1 & EA2_nchar==EA2_nchar & EA1_IorD], keep := FALSE]
  h[idx[EA1_eq_OA1], keep := FALSE]
  h[idx[EA2_eq_OA2], keep := FALSE]

  # return
  return(h)
}
