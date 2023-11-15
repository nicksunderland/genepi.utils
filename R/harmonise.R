
#' @title Harmonise GWAS effects
#' @description
#' Harmonise effects across two GWAS datasets. The `gwas` objects are expected to be in standard format - see `standardise_gwas()`.
#' Corrects strand for non-palindromic SNPs (flip BETAs). Drop palindromic variants. Drop non-matching variant alleles. \cr
#' This function is based on the same function in the [Slopehunter](https://github.com/Osmahmoud/SlopeHunter/blob/master/R/harmonise_effects.R)
#' and [TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR/blob/master/R/harmonise.R) packages. I have re-written it here to use
#' the standardised column names and use pure data.table syntax for increased processing speed.
#' @references [Slopehunter](https://github.com/Osmahmoud/SlopeHunter/blob/master/R/harmonise_effects.R) https://github.com/Osmahmoud/SlopeHunter/blob/master/R/harmonise_effects.R
#' @references [TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR/blob/master/R/harmonise.R) https://github.com/MRCIEU/TwoSampleMR/blob/master/R/harmonise.R
#' @param gwas1 a data.frame like object, GWAS 1
#' @param gwas2 a data.frame like object, GWAS 2
#' @return a data.table, the harmonised dataset
#' @export
#'
harmonise <- function(gwas1, gwas2, gwas1_trait="incidence", gwas2_trait="progression", merge="chrpos") {

  # checks
  merge <- match.arg(merge, c("chrpos", "snp"))
  stopifnot("`gwas[1|2]` must be a data.frame like object" = inherits(gwas1, "data.frame") & inherits(gwas2, "data.frame"))
  req_cols <- c("SNP","CHR","BP","EA","OA","EAF","BETA","P")
  stopifnot("At least columns EA, OA, EAF, BETA, and P must be present in `gwas1`" = all(req_cols %in% colnames(gwas1)))
  stopifnot("At least columns EA, OA, EAF, BETA, and P must be present in `gwas2`" = all(req_cols %in% colnames(gwas2)))

  # to data.table format
  gwas1 <- data.table::as.data.table(gwas1)
  gwas2 <- data.table::as.data.table(gwas2)

  # trim columns to only those needed
  gwas1[, names(gwas1)[!names(gwas1) %in% req_cols] := NULL]
  gwas2[, names(gwas2)[!names(gwas2) %in% req_cols] := NULL]

  # join
  if(merge=="chrpos") {

    data.table::setkey(gwas1, CHR, BP)
    data.table::setkey(gwas2, CHR, BP)
    h <- data.table::merge.data.table(gwas1, gwas2, all=FALSE, suffixes=c("_1","_2"))

  } else if(merge=="chrpos") {

    data.table::setkey(gwas1, SNP)
    data.table::setkey(gwas2, SNP)
    h <- data.table::merge.data.table(gwas1, gwas2, all=FALSE, suffixes=c("_1","_2"))

  }



  # report loss of variants
  if(nrow(h)==0) {

    stop("No matching variants between gwas inputs")

  } else if(nrow(h) != nrow(gwas1) | nrow(h) != nrow(gwas2)) {

    msg <- paste0("Harmonisation: removed ", max(c(nrow(gwas1), nrow(gwas2))) - nrow(h), " variants. ",
                  "gwas1 (n=", nrow(gwas1),")", " gwas2 (n=", nrow(gwas2),")", " harmonised (n=", nrow(h),")\n")
    warning(msg)

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
  h[, to_swap := EA_1==OA_2 & OA_1==EA_2]
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
  h[        , OK   := EA_1==EA_2 & OA_1==OA_2]
  h[        , ii   := !palindromic & !OK]
  h[ii==TRUE, EA_2 := flip_alleles(EA_2)]
  h[ii==TRUE, OA_2 := flip_alleles(OA_2)]
  jinfo[['flipped_alleles']] <- sum(h$ii)

  # If still DON'T match (e.g. A/C Vs. G/T, which - after flipping in the previous step - becomes A/C Vs. C/A) then try swapping
  h[                        , OK      := EA_1==EA_2 & OA_1==OA_2]
  h[                        , ii      := !palindromic & !OK]
  h[                        , to_swap := EA_1==OA_2 & OA_1==EA_2]
  h[ii==TRUE & to_swap==TRUE, EA2tmp  := EA_2]
  h[ii==TRUE & to_swap==TRUE, EA_2    := OA_2]
  h[ii==TRUE & to_swap==TRUE, OA_2    := EA2tmp]
  h[ii==TRUE & to_swap==TRUE, BETA_2  := BETA_2 * -1]
  h[ii==TRUE & to_swap==TRUE, EAF_2   := 1 - EAF_2]
  jinfo[['flipped_then_switched_alleles']] <- sum(h$ii & h$to_swap)

  # Any SNPs left with un-matching alleles (e.g. A/C Vs. A/G) need to be removed
  h[, OK   := EA_1==EA_2 & OA_1==OA_2]
  h[, keep := keep & OK & !palindromic]

  # clean up
  h[, EA2tmp  := NULL]
  h[, ii      := NULL]
  h[, OK      := NULL]
  h[, to_swap := NULL]

  # set output names
  data.table::setnames(h, names(h), sub("1$", gwas1_trait, names(h)))
  data.table::setnames(h, names(h), sub("2$", gwas2_trait, names(h)))

  # set logging attribute
  attr(h, "info") <- jinfo

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


flip_alleles <- function(x)
{
  x <- toupper(x)
  x <- gsub("C", "G", x)
  x <- gsub("G", "C", x)
  x <- gsub("A", "T", x)
  x <- gsub("T", "A", x)
  return(x)
}

recode_indels <- function(h)
{
  # logical mask, whether an indel or not
  indel_flag <- nchar(h[,EA_1]) > 1 | nchar(h[,OA_1]) > 1 | h[,EA_1]=="D" | h[,EA_1]=="I"

  # if no indels nothing needed
  if(sum(indel_flag)==0) return(h)

  # indel row indices
  idx <- which(indel_flag)

  # row indices logical mask for each criteria
  EA1_gt_OA1 <- nchar(h[idx, EA_1]) > nchar(h[idx, OA_1])
  EA1_lt_OA1 <- nchar(h[idx, EA_1]) < nchar(h[idx, OA_1])
  EA2_gt_OA2 <- nchar(h[idx, EA_2]) > nchar(h[idx, OA_2])
  EA2_lt_OA2 <- nchar(h[idx, EA_2]) < nchar(h[idx, OA_2])
  EA2_is_I   <- h[idx, EA_2] == "I" & h[idx, OA_2] == "D"
  EA2_is_D   <- h[idx, EA_2] == "D" & h[idx, OA_2] == "I"
  EA1_is_I   <- h[idx, EA_1] == "I" & h[idx, OA_1] == "D"
  EA1_is_D   <- h[idx, EA_1] == "D" & h[idx, OA_1] == "I"

  # recode
  h[ idx[EA1_gt_OA1 & EA2_is_I], c("EA_2","OA_2") := list(EA_1,OA_1)] # B1[i1] <- A1[i1] ; B2[i1] <- A2[i1]
  h[ idx[EA1_lt_OA1 & EA2_is_I], c("EA_2","OA_2") := list(OA_1,EA_1)] # B1[i1] <- A2[i1] ; B2[i1] <- A1[i1]
  h[ idx[EA1_gt_OA1 & EA2_is_D], c("EA_2","OA_2") := list(OA_1,EA_1)] # B1[i1] <- A2[i1] ; B2[i1] <- A1[i1]
  h[ idx[EA1_lt_OA1 & EA2_is_D], c("EA_2","OA_2") := list(EA_1,OA_1)] # B1[i1] <- A1[i1] ; B2[i1] <- A2[i1]
  h[ idx[EA2_gt_OA2 & EA1_is_I], c("EA_1","OA_1") := list(EA_2,OA_2)] # A1[i1] <- B1[i1] ; A2[i1] <- B2[i1]
  h[ idx[EA2_lt_OA2 & EA1_is_I], c("OA_1","EA_1") := list(EA_2,OA_2)] # A2[i1] <- B1[i1] ; A1[i1] <- B2[i1]
  h[ idx[EA2_gt_OA2 & EA1_is_D], c("OA_1","EA_1") := list(EA_2,OA_2)] # A2[i1] <- B1[i1] ; A1[i1] <- B2[i1]
  h[ idx[EA2_lt_OA2 & EA1_is_D], c("EA_1","OA_1") := list(EA_2,OA_2)] # A1[i1] <- B1[i1] | A2[i1] <- B2[i1]

  # flag the ones not to keep
  h[idx[nchar(EA_1)>1 & nchar(EA_1)==nchar(OA_1) & (EA_2=="D" | EA_2=="I")], keep := FALSE]
  h[idx[nchar(EA_2)>1 & nchar(EA_2)==nchar(OA_2) & (EA_1=="D" | EA_1=="I")], keep := FALSE]
  h[idx[EA_1==OA_1], keep := FALSE]
  h[idx[EA_2==OA_2], keep := FALSE]

  # return
  return(h)
}
