# Silence R CMD check
globalVariables(c(),
                package = "genepi.utils")

#' @title MR object
#' @description
#' An MR object is a container for vectors and matrices of 2 or more GWAS data. \cr\cr
#'
#' @param exposure a `GWAS` object or list of `GWAS` objects
#' @param outcome a `GWAS` object
#' @param harmonise_strictness an integer (1,2,3) corresponding to the TwoSampleMR harmonisation options of the same name.
#' @param with_proxies a logical, whether to search for outcome proxy SNPs if the original variant is not present in the data.
#' @param proxy_r2 a numeric (0-1), the minimum r2 threshold used to define a suitable proxy variant.
#' @param proxy_kb an integer, the maximum distance in kb that a proxy variant can be from the index variant.
#' @param proxy_window an integer, the maximum number of variants to search either side of the index variant.
#' @param proxy_eaf a numeric, the minimum minor allele frequency required to define a suitable proxy.
#'
#' @slot snps character, variant ID
#' @slot chr character, chromosome identifier
#' @slot bp integer, base position
#' @slot ea character, effect allele
#' @slot oa character, other allele
#' @slot eafx numeric, exposure effect allele frequency
#' @slot nx integer, exposure total number of samples
#' @slot ncasex integer, exposure number of cases
#' @slot bx numeric, exposure effect size
#' @slot bxse numeric, exposure effect size standard error
#' @slot px numeric, exposure p-value
#' @slot eafy numeric, exposure effect allele frequency
#' @slot ny integer, exposure total number of samples
#' @slot ncasey integer, exposure number of cases
#' @slot by numeric, exposure effect size
#' @slot byse numeric, exposure effect size standard error
#' @slot py numeric, exposure p-value
#' @slot exposure_id character, the GWAS identifier
#' @slot exposure character, the GWAS exposure
#' @slot outcome_id character, the GWAS identifier
#' @slot outcome character, the GWAS outcome
#' @slot group integer, grouping variable used for plotting
#' @slot index_snp logical, whether the variant is an index variant (via clumping)
#' @slot proxy_snp character, the id of the proxy snp
#' @slot ld_info logical, whether there is LD information
#' @slot correlation matrix, a correlation matrix of signed R values between variants
#'
#' @return an S7 class genepi.utils::MR object
#'
#' @import S7
#' @export
MR <- new_class(
  #==============================
  # MR class name
  #==============================
  name    = "MR",
  package = "genepi.utils",

  #==============================
  # MR class properties
  #==============================
  properties = list(
    #----------------------------
    # column names / data vectors
    #----------------------------
    snps        = class_character,
    chr         = class_character,
    bp          = class_integer,
    ea          = class_character,
    oa          = class_character,
    eafx        = new_S3_class('matrix'),
    nx          = new_S3_class('matrix'),
    ncasex      = new_S3_class('matrix'),
    bx          = new_S3_class('matrix'),
    bxse        = new_S3_class('matrix'),
    px          = new_S3_class('matrix'),
    eafy        = class_numeric,
    ny          = class_integer,
    ncasey      = class_integer,
    by          = class_numeric,
    byse        = class_numeric,
    py          = class_numeric,
    exposure_id = class_character,
    exposure    = class_character,
    outcome_id  = class_character,
    outcome     = class_character,
    group       = class_integer,
    index_snp   = class_logical,
    proxy_snp   = class_character,
    ld_info     = class_logical,
    #----------------------------
    # correlation matrix
    #----------------------------
    correlation = new_S3_class('matrix')
  ),

  #==============================
  # MR class constructor func.
  #==============================
  constructor = function(exposure,
                         outcome,
                         harmonise_strictness = 2,
                         with_proxies         = FALSE,
                         proxy_r2             = 0.8,
                         proxy_kb             = 250,
                         proxy_window         = 15,
                         proxy_eaf            = 0.01,
                         correlation          = NULL,
                         verbose              = TRUE) {

    if(verbose) message("Loading MR ...")

#
#     # testing
#     mapping <- list(
#       "rsid"=    list("name"= "rsid",    "alias"= "MARKER",        "type"= "character"),
#       "chr"=     list("name"= "chr",     "alias"= "CHR",           "type"= "character"),
#       "bp"=      list("name"= "bp",      "alias"= "POS",           "type"= "integer"),
#       "se"=      list("name"= "se",      "alias"= "SE",            "type"= "numeric"),
#       "p"=       list("name"= "p",       "alias"= "P",             "type"= "numeric"),
#       "strand"=  list("name"= "strand",  "alias"= "STRAND",        "type"= "character"),
#       "n"=       list("name"= "n",       "alias"= "N_SAMPLE",      "type"= "integer"),
#       "ncase"=   list("name"= "ncase",   "alias"= "N_EVENTS",      "type"= "integer"),
#       "imputed"= list("name"= "imputed", "alias"= "Imputed",       "type"= "logical"),
#       "info"=    list("name"= "info",    "alias"= "QUAL_SCORE",    "type"= "numeric"),
#       "ea"=      list("name"= "ea",      "alias"= c("ea", "ORI_EFFECT_ALLELE", "EFFECT_ALLELE"), "type"= "character"),
#       "oa"=      list("name"= "oa",      "alias"= c("oa", "ORI_OTHER_ALLELE", "OTHER_ALLELE"),   "type"= "character"),
#       "beta"=    list("name"= "beta",    "alias"= c("ORI_BETA", "BETA"),                   "type"= "numeric"),
#       "eaf"=     list("name"= "eaf",     "alias"= c("ORI_EAF", "EAF"),                     "type"= "numeric")
#     )
#     map <- ColumnMap(lapply(mapping, function(x) do.call(genepi.utils::Column, x)))
#     dat = "/Users/xx20081/Documents/local_data/hermes_progression/biostat_disc/raw/BIOSTAT_Discovery.females.allcause.gz"
#     exposure = GWAS(dat,map=map,fill=T,drop=T,fill_rsid = FALSE,missing_rsid = "fill_CHR:BP")
#     outcome = GWAS(dat,map=map,fill=T,drop=T,fill_rsid = FALSE,missing_rsid = "fill_CHR:BP")
#     verbose = TRUE
#     harmonise_strictness=2
#     with_proxies=F
    # correlation = UKBB_GLP1R_LDMAT
    # proxy_r2=0.8
    # proxy_kb=250
    # proxy_window=15
    # proxy_eaf=0.01



    # checks
    if (!is.list(exposure)) exposure <- list(exposure)
    stopifnot("Exposure input list must all be of GWAS class" = all(sapply(exposure, function(x) inherits(x, "genepi.utils::GWAS"))))
    stopifnot("Outcome must be of GWAS class" = inherits(outcome, "genepi.utils::GWAS"))

    # format exposure
    if(verbose) message("[i] processing exposure(s)")
    e <- as.twosample.mr(exposure, "exposure")

    # pre-harmonise exposures if >1 exposure
    if(length(exposure) > 1) {

      # have the reference dataset as the first one in the table
      ref <- e[id.exposure==id.exposure[1], ]

      # TwoSampleMR::harmonise_data flips the outcome data, so have the full data as `outcome` type initially
      names(e) <- gsub("exposure", "outcome", names(e))

      # harmonise exposures against first exposure (ref) dataset
      e <- TwoSampleMR::harmonise_data(ref, e, action = harmonise_strictness)
      e <- subset(e, mr_keep)

      # only keep harmonised SNPs common to all exposures
      tab <- table(e$SNP)
      keepsnps <- names(tab)[tab == length(unique(e$id.outcome))]
      e <- e[e$SNP %in% keepsnps, ]

      # drop dummy `exposure` columns and rename outcome to exposure
      e <- e[, !grepl("exposure|action|mr_keep|SNP_index", names(e))]
      names(e) <- gsub("outcome","exposure", names(e))

      # as a data.table
      e <- data.table::as.data.table(e)
    }

    # now get the real outcome data
    if(verbose) message("[i] processing outcome")
    o <- as.twosample.mr(outcome, "outcome")

    # ?get proxies
    if(with_proxies) {
      warning("`with_proxies` not implemented yet")
      # stopifnot("LD matrix must be supplied for proxy search" = !is.null(correlation))
      # o <- get_proxies(e, o, correlation, proxy_r2, proxy_kb, proxy_window, proxy_eaf)
    }

    # harmonise the (pre-harmonised) exposure data with the outcome data
    if(verbose) message("[i] harmonising GWASs")
    h <- TwoSampleMR::harmonise_data(e, o, action = harmonise_strictness) |> data.table::as.data.table()
    h <- h[mr_keep == TRUE, ]

    # rename SNPs now harmonised
    h[, SNP := sub("(.*)_(?:[AaCcTtGg]+|[DdIi])_(?:[AaCcTtGg]+|[DdIi])$", "\\1", SNP)]
    h[, SNP := paste0(SNP, "_", other_allele.outcome, "_", effect_allele.outcome)]

    # exposure matrices
    exp_beta  <- data.table::dcast(h, SNP + chr.exposure + pos.exposure + effect_allele.exposure + other_allele.exposure ~ id.exposure, value.var = "beta.exposure")
    exp_pval  <- data.table::dcast(h, SNP ~ id.exposure, value.var='pval.exposure')
    exp_se    <- data.table::dcast(h, SNP ~ id.exposure, value.var='se.exposure')
    exp_eaf   <- data.table::dcast(h, SNP ~ id.exposure, value.var='eaf.exposure')
    exp_n     <- data.table::dcast(h, SNP ~ id.exposure, value.var='samplesize.exposure')
    exp_ncase <- data.table::dcast(h, SNP ~ id.exposure, value.var='ncase.exposure')

    # ?find proxies
    if (with_proxies && "proxy_snp.outcome" %in% names(h)) {
      warning("`with_proxies` not implemented yet")
      proxy_snps <- h[id.exposure==id.exposure[[1]], proxy_snp.outcome]
    } else {
      proxy_snps <- character()
    }

    # process correlation matrix
    if (!is.null(correlation)) {
      warning("`correlation` not implemented yet")
      # if(!is.null(correlation)) {
      #   corr_mat <- set_ld_mat(object, correlation)
      # }
    } else {
      corr_mat <- matrix()
      ld_info  <- rep(FALSE, length(exp_beta$SNP))
    }

    # assign to the class object
    object <- new_object(S7::S7_object(),
                         snps        = exp_beta$SNP,
                         chr         = exp_beta$chr.exposure,
                         bp          = as.integer(exp_beta$pos.exposure),
                         ea          = exp_beta$effect_allele.exposure,
                         oa          = exp_beta$other_allele.exposure,
                         eafx        = as.matrix(exp_eaf[,-1]),
                         nx          = as.matrix(exp_n[,-1]),
                         ncasex      = as.matrix(exp_ncase[,-1]),
                         bx          = as.matrix(exp_beta[,-c(1:5)]),
                         bxse        = as.matrix(exp_se[,-1]),
                         px          = as.matrix(exp_pval[,-1]),
                         eafy        = h[id.exposure==id.exposure[[1]], eaf.outcome],
                         ny          = as.integer(h[id.exposure==id.exposure[[1]], samplesize.outcome]),
                         ncasey      = as.integer(h[id.exposure==id.exposure[[1]], ncase.outcome]),
                         by          = h[id.exposure==id.exposure[[1]], beta.outcome],
                         byse        = h[id.exposure==id.exposure[[1]], se.outcome],
                         py          = h[id.exposure==id.exposure[[1]], pval.outcome],
                         exposure_id = h[, unique(id.exposure), by=id.exposure]$V1,
                         exposure    = h[, unique(exposure), by=id.exposure]$V1,
                         outcome_id  = h$id.outcome[[1]],
                         outcome     = h$outcome[[1]],
                         group       = integer(),
                         index_snp   = rep(TRUE, length(exp_beta$SNP)),
                         proxy_snp   = proxy_snps,
                         ld_info     = ld_info,
                         correlation = corr_mat)

    # return the object
    return(object)
  },

  #==============================
  # MR class validator func.
  #==============================
  validator = function(self) {

    # stopifnot("Unequal vector lengths" = sapply(lengths(list(self@snps, self@chr, self@bp, self@ea, self@oa, self@eafx[,1], self@bx[,1], self@bxse[,1], self@px[,1], self@eafy, self@by, self@byse, self@py)), function(x) x == length(self@snps)))
    # stopifnot("Invalid sample size `n`"       = length(self@nx[,1]) == length(self@snps))
    # stopifnot("Invalid sample size `ncase`"   = length(self@ncasex) == length(self@snps))
  }
  # old S4 class - need to adapt
  # correct lengths
  # l <- list(object@snps, object@ea, object@oa, object@chr, object@bp, object@eafx[,1], object@bx[,1], object@bxse[,1], object@px[,1], object@by, object@byse, object@py)
  # stopifnot("Inconsistent lengths" = all(lengths(l)==length(l[[1]])))
  #
  # # if LD matrix ensure symmetric and same size
  # if(!all(is.na(object@correlation))) {
  #   stopifnot("LD matrix names do not match" = all(object@snps==rownames(object@correlation)))
  #   stopifnot("LD matrix doesn't match" = nrow(object@correlation)==length(l[[1]]))
  #   stopifnot("LD info vector inconsistent length" = length(object@ld_info)==length(l[[1]]))
  #   stopifnot("LD matrix must be symmetric" = isSymmetric(object@correlation))
  # }
  #
  # # all positive might be wrong / r2 values
  # if(min(object@correlation, na.rm=T)>0) {
  #   message("LD matrix - all values found to be +ve, please check matrix is `r` not `r2` values")
  # }
)


#' @title Test if MR object is multivariable
#' @noRd
is_multivariable <- new_generic("is_multivariable", "x", function(x) { S7_dispatch() })
method(is_multivariable, MR) <- function(x) { ncol(x@bx) > 1 }


#' @title Get number of exposures loaded in MR object
#' @noRd
num_exposures <- new_generic("num_exposures", "x", function(x) { S7_dispatch() })
method(num_exposures, MR) <- function(x) { ncol(x@bx) }


#' @title Strip variant IDs of any trailing allele information
#' @noRd
snps_no_alleles <- new_generic("snps_no_alleles", "x", function(x) { S7_dispatch() })
method(snps_no_alleles, MR) <- function(x) {
  return( sub("(.*)_(?:[ACTG]+|[DI])_(?:[ACTG]+|[DI])$", "\\1", x@snps, ignore.case = TRUE) )
}

#' @title Reset index SNP
#' @export
reset_index_snp <- new_generic("reset_index_snp", "x", function(x) { S7_dispatch() })
method(reset_index_snp, MR) <- function(x) {
  x@index_snp <- rep(TRUE, length(x@index_snp))
  return(x)
}


#' @title Convert to MendelianRandomization::MRInput object
#' @export
to_MRInput <- new_generic("to_MRInput", "x", function(x, corr = FALSE) { S7_dispatch() })
#' @name to_MRInput
method(to_MRInput, MR) <- function(x, corr = FALSE) {

  # check LD matrix exists
  if(corr && length(mr_obj@ld_info) == 0) warning("Correlation MR requested but no LD info / matrix found")

  mr_input <- MendelianRandomization::mr_input(
    bx          = x@bx[which(x@index_snp == TRUE), 1],
    bxse        = x@bxse[which(x@index_snp == TRUE), 1],
    by          = x@by[which(x@index_snp == TRUE)],
    byse        = x@byse[which(x@index_snp == TRUE)],
    exposure    = x@exposure[1],
    outcome     = x@outcome,
    snps        = x@snps[which(x@index_snp == TRUE)],
    correlation = if(!corr || length(x@ld_info)==0) { matrix() } else { as.matrix(x@correlation[which(x@index_snp == TRUE), which(x@index_snp == TRUE)]) }  # as.matrix in case just one value
  )

  return(mr_input)
}


#' @title Convert to MendelianRandomization::MRMVInput object
#' @export
to_MRMVInput <- new_generic("to_MRMVInput", "x", function(x, corr = FALSE) { S7_dispatch() })
#' @name to_MRMVInput
method(to_MRMVInput, MR) <- function(x, corr = FALSE) {

  # check LD matrix exists
  if(corr && length(x@ld_info) == 0) warning("Correlation MR requested but no LD info / matrix found")

  mr_input <- MendelianRandomization::mr_mvinput(
    bx          = x@bx[which(x@index_snp == TRUE), ],
    bxse        = x@bxse[which(x@index_snp == TRUE), ],
    by          = x@by[which(x@index_snp == TRUE)],
    byse        = x@byse[which(x@index_snp == TRUE)],
    exposure    = x@exposure,
    outcome     = x@outcome,
    snps        = x@snps[which(x@index_snp == TRUE)],
    correlation = if(!corr || length(x@ld_info)==0) { matrix() } else { as.matrix(x@correlation[which(x@index_snp == TRUE),which(x@index_snp == TRUE)]) }
  )

  return(mr_input)
}


#' @title Run MR
#' @export
mr <- new_generic("mr", "x", function(x, corr = FALSE, methods=c('mr_ivw','mr_egger','mr_weighted_median','mr_weighted_mode'), ...) { S7_dispatch() })
#' @name mr
method(mr, MR) <- function(x, corr = FALSE, methods=c('mr_ivw','mr_egger','mr_weighted_median','mr_weighted_mode'), ...) {

  res <- lapply(methods, function(method) do.call(method, args=c(list(x=x, corr=corr), list(...)))) |> `names<-`(methods)
  dat <- mr_results_to_data_table(res)
  return(dat)

}


#' @title Run IVW MR
#' @export
mr_ivw <- new_generic("mr_ivw", "x", function(x, corr = FALSE, ...) { S7_dispatch() })
#' @name mr_ivw
method(mr_ivw, MR) <- function(x, corr = FALSE, ...) {

  tryCatch({
    if(is_multivariable(x)) {
      method <- "mr_mvivw"
      res    <- MendelianRandomization::mr_mvivw(to_MRMVInput(x, corr), nx=apply(x@nx,2,max,na.rm=TRUE), ...)
    } else {
      method <- "mr_ivw"
      res <- MendelianRandomization::mr_ivw(to_MRInput(x, corr), ...)
    }
  },
  error=function(e) {
    warning(paste0("`mr_",ifelse(is_multivariable(x),"mv",""),"ivw(x, corr=",corr,")` failed: ", conditionMessage(e)))
    res <- NULL
  },
  finally={
    return(MRResult(res=res, mr_obj=x, method=method))
  })

}


#' @title Run Egger MR
#' @export
mr_egger <- new_generic("mr_egger", "x", function(x, corr = FALSE, ...) { S7_dispatch() })
#' @name mr_egger
method(mr_egger, MR) <- function(x, corr = FALSE, ...) {

  tryCatch({

    # if(sum(x@index_snp, na.rm=TRUE)) stop("less than 3 data points detected, unable to run.") # leave errors to MR package

    if(is_multivariable(x)) {
      method <- "mr_mvegger"
      res    <- MendelianRandomization::mr_mvegger(to_MRMVInput(x, corr), ...)
    } else {
      method <- "mr_egger"
      res    <- MendelianRandomization::mr_egger(to_MRInput(x, corr), ...)
    }
  },
  error=function(e) {
    warning(paste0("`mr_",ifelse(is_multivariable(x),"mv",""),"egger(x, corr=",corr,")` failed: ", conditionMessage(e)))
    res <- NULL
  },
  finally={
    return(MRResult(res=res, mr_obj=x, method=method))
  })

}


#' @title Run weighted median MR
#' @export
mr_weighted_median <- new_generic("mr_weighted_median", "x", function(x, corr = FALSE, ...) { S7_dispatch() })
#' @name mr_egger
method(mr_weighted_median, MR) <- function(x, corr = FALSE, ...) {

  tryCatch({

    # if(sum(x@index_snp, na.rm=TRUE)) stop("less than 3 data points detected, unable to run.")# leave errors to MR package

    if(is_multivariable(x)) {
      method <- "mr_mvweighted_median"
      res    <- MendelianRandomization::mr_mvmedian(to_MRMVInput(x, corr), ...)
      return(MRResult(res=res, mr_obj=x, ))
    } else {
      method = "mr_weighted_median"
      res <- MendelianRandomization::mr_median(to_MRInput(x, corr), weighting="weighted", ...)
    }
  },
  error=function(e) {
    warning(paste0("`mr_",ifelse(is_multivariable(x),"mv",""),"weighted_median(x, corr=",corr,")` failed: ", conditionMessage(e)))
    res <- NULL
  },
  finally={
    return(MRResult(res=res, mr_obj=x, method=method))
  })

}


#' @title Run weighted mode MR
#' @export
mr_weighted_mode <- new_generic("mr_weighted_mode", "x", function(x, corr = FALSE, ...) { S7_dispatch() })
#' @name mr_weighted_mode
method(mr_weighted_mode, MR) <- function(x, corr = FALSE, ...) {

  tryCatch({
    if(is_multivariable(x)) {
      warning("No multivariable mode based function")
      method <- "mr_mvweighted_mode"
      res    <- NULL
    } else {
      method <- "mr_weighted_mode"
      res    <- MendelianRandomization::mr_mbe(to_MRInput(x, corr), weighting='weighted', ...)
    }
  },
  error=function(e) {
    warning(paste0("`mr_",ifelse(is_multivariable(x),"mv",""),"weighted_mode(x, corr=",corr,")` failed: ", conditionMessage(e)))
    res <- NULL
  },
  finally={
    return(MRResult(res=res, mr_obj=x, method=method))
  })

}


#' @title Run PC-GMM MR
#' @export
mr_pcgmm <- new_generic("mr_pcgmm", "x", function(x, corr = TRUE, ...) { S7_dispatch() })
#' @name mr_pcgmm
method(mr_pcgmm, MR) <- function(x, corr = TRUE, ...) {

  tryCatch({
    if(!corr) warning("corr=FALSE ignored as LD correlation matrix must alway be set with mr_pcgmm")

    if(is_multivariable(x)) {
      method <- "mr_mvpcgmm"
      res    <- MendelianRandomization::mr_mvpcgmm(to_MRMVInput(x, TRUE), nx=apply(x@nx,2,max,na.rm=TRUE), ny=max(x@ny,na.rm=TRUE), ...)
    } else {
      method <- "mr_pcgmm"
      res    <- MendelianRandomization::mr_pcgmm(to_MRInput(x, TRUE), nx=max(x@nx[,1],na.rm=TRUE), ny=max(x@ny,na.rm=TRUE), ...)
    }
  },
  error=function(e) {
    warning(paste0("`mr_",ifelse(is_multivariable(x),"mv",""),"pcgmm(x)` failed: ", conditionMessage(e)))
    res <- NULL
  },
  finally={
    return(MRResult(res=res, mr_obj=x, method=method))
  })

}


#' @title Helper MRResult class
#' @noRd
MRResult <- new_class(
  #==============================
  # MRResult class name
  #==============================
  name    = "MRResult",
  package = "genepi.utils",

  #==============================
  # MRResult class properties
  #==============================
  properties = list(
    #----------------------------
    # column names / data vectors
    #----------------------------
    method         = class_character,
    correlation    = class_logical,
    exposure       = class_character,
    outcome        = class_character,
    n_snp          = class_integer,
    b              = class_numeric,
    b_se           = class_numeric,
    p              = class_numeric,
    intercept      = class_numeric,
    int_se         = class_numeric,
    int_p          = class_numeric,
    qstat          = class_numeric,
    qstat_p        = class_numeric,
    fstat          = class_numeric,
    condfstat      = class_numeric,
    overdispersion = class_numeric,
    n_pc           = class_integer,
    n_hunted       = class_integer,
    slopehunter_pi = class_numeric,
    slopehunter_ent= class_numeric
  ),
  #==============================
  # MRResult class constructor func.
  #==============================
  constructor = function(res = NULL, mr_obj = NULL, method = NULL) {

    object <- new_object(S7::S7_object(),
                         method         = ifelse(is.null(method), NA_character_, method),
                         correlation    = tryCatch({ifelse(all(is.na(res@Correlation)),F,T)}, error=function(e){ F }),
                         exposure       = tryCatch({mr_obj@exposure},   error=function(e){ NA_character_ }),
                         outcome        = tryCatch({mr_obj@outcome},    error=function(e){ NA_character_ }),
                         n_snp          = tryCatch({res@SNPs},          error=function(e){ NA_integer_ }),
                         b              = tryCatch({res@Estimate},      error=function(e){ NA_real_ }),
                         b_se           = tryCatch({res@StdError},      error=function(e){ tryCatch({ res@StdError.Est }, error=function(e){ NA_real_ }) }),
                         p              = tryCatch({res@Pvalue},        error=function(e){ tryCatch({ res@Pvalue.Est   }, error=function(e){ NA_real_ }) }),
                         intercept      = tryCatch({res@Intercept},     error=function(e){ 0 }),
                         int_se         = tryCatch({res@StdError.Int},  error=function(e){ NA_real_ }),
                         int_p          = tryCatch({res@Pvalue.Int},    error=function(e){ NA_real_ }),
                         qstat          = tryCatch({res@Heter.Stat[1]}, error=function(e){ NA_real_ }),
                         qstat_p        = tryCatch({res@Heter.Stat[2]}, error=function(e){ NA_real_ }),
                         fstat          = tryCatch({res@Fstat},         error=function(e){ NA_real_ }),
                         condfstat      = tryCatch({res@CondFstat},     error=function(e){ NA_real_ }),
                         overdispersion = tryCatch({res@Overdispersion},error=function(e){ NA_real_ }),
                         n_pc           = tryCatch({res@PCs },          error=function(e){ NA_integer_ }),
                         n_hunted       = tryCatch({res@n_hunted},      error=function(e){ NA_integer_ }),
                         slopehunter_pi = tryCatch({res@slopehunter_pi},error=function(e){ NA_real_ }),
                         slopehunter_ent= tryCatch({res@slopehunter_ent},error=function(e){ NA_real_ }))

    # return the object
    return(object)
  }
)


#' @title MR results to data.table
#' @param x MRResult object to covert to data.table
#' @param ... argument for data.table generic, ignored in this implementation
#' @export
mr_results_to_data_table <- new_generic("mr_results_to_data_table", "x", function(x) { S7_dispatch() })
#' @name mr_results_to_data_table
method(mr_results_to_data_table, MRResult) <- function(x) {

  d <- data.table::data.table(
    method      = x@method,
    correlation = x@correlation,
    exposure    = x@exposure,
    outcome     = x@outcome,
    n_snp       = x@n_snp,
    b           = x@b,
    b_se        = x@b_se,
    p           = x@p,
    intercept   = x@intercept,
    int_se      = x@int_se,
    int_p       = x@int_p,
    qstat       = x@qstat,
    qstat_p     = x@qstat_p,
    fstat       = x@fstat,
    condfstat   = x@condfstat,
    overdispersion = x@overdispersion,
    n_pc        = x@n_pc,
    n_hunted    = x@n_hunted,
    slopehunter_pi = x@slopehunter_pi,
    slopehunter_ent= x@slopehunter_ent
  )
}

#' @name mr_results_to_data_table
method(mr_results_to_data_table, class_list) <- function(x) {
  stopifnot("Should be a list of MRResult objects" = sapply(x, function(r) inherits(r, "genepi.utils::MRResult")))
  lapply(x, function(r) mr_results_to_data_table(r)) |> data.table::rbindlist()
}


#' @title Set the LD matrix
#' @export
set_ld_mat <- new_generic("set_ld_mat", "x", function(x, correlation) { S7_dispatch() })
method(set_ld_mat, MR) <- function(x, correlation) {

  # LD variants info
  ld_info <- data.table::data.table(
    idx  = 1:length(rownames(correlation)),
    rsid = sub("^(rs[0-9]+|[0-9XY]+:[0-9]+).*","\\1",rownames(correlation)),
    oa   = sub("^(?:rs[0-9]+|[0-9XY]+:[0-9]+)_([ACTG]+|[DI])_.*","\\1",rownames(correlation)),
    ea   = sub("^(?:rs[0-9]+|[0-9XY]+:[0-9]+)_(?:[ACTG]+|[DI])_([ACTG]+|[DI]).*","\\1",rownames(correlation))
  )
  stopifnot("Incorrect correlation matrix naming detected" = all(grepl("^(rs[0-9]+|[0-9XY]+:[0-9]+)$", ld_info$rsid) &
                                                                   grepl("^([ACTG]+|[DI])$", ld_info$oa) &
                                                                   grepl("^([ACTG]+|[DI])$", ld_info$ea)))

  # current object's (harmonised) mr datasets (just take first exposure if multiple)
  current <- as.data.table(x, exposure=1)

  # strip the alleles from the RSID
  current[, rsid := snps_no_alleles(x)]

  # join and flag those that need flipping
  current[ld_info, c('flip','ld_mat_idx') := list(FALSE, i.idx), on=c('rsid','ea','oa')]
  current[ld_info, c('flip','ld_mat_idx') := list(TRUE,  i.idx) , on=c('rsid'='rsid','ea'='oa','oa'='ea')]

  # flag those with LD data
  x@ld_info <- !is.na(current$flip)

  # the indices of those that need flipping
  flip <- which(current$flip==TRUE)

  # do the flipping
  ea_store       <- x@ea
  x@ea[flip]     <- x@oa[flip]
  x@oa[flip]     <- ea_store[flip]
  x@snps         <- paste0(snps_no_alleles(x),"_",x@ea,"_",x@oa)
  x@eafx[flip, ] <- as.matrix(1-x@eafx[flip, ])
  x@bx[flip, ]   <- as.matrix(x@bx[flip, ]*-1)
  x@by[flip]     <- x@by[flip]*-1
  x@eafy[flip]   <- 1-x@eafy[flip]

  # set the LD matrix and ensure order of ld matrix matches data (this will probably have NAs - but deal with at use time)
  x@correlation <- correlation[current$ld_mat_idx, current$ld_mat_idx]
  rownames(x@correlation) <- colnames(x@correlation) <- x@snps

  # return
  return(x)
}


#' @title Clump MR object exposure
#' @export
clump_mr <- new_generic("clump_mr", "x", function(x,
                                            p1 = 1,
                                            p2 = 1,
                                            r2 = 0.001,
                                            kb = 250,
                                            plink2    = genepi.utils::which_plink2(),
                                            plink_ref = genepi.utils::which_1000G_reference(build="GRCh37")) { S7_dispatch() })
#' @name clump_mr
method(clump_mr, MR) <- function(x,
                              p1 = 1,
                              p2 = 1,
                              r2 = 0.001,
                              kb = 250,
                              plink2    = genepi.utils::which_plink2(),
                              plink_ref = genepi.utils::which_1000G_reference(build="GRCh37")) {

  # extract data to work with - might be multiple exposures with multivariable MR
  if(ncol(x@bx)==1) {

    d <- as.data.table(x, exposure=1)

  } else {
    warning('Experimental, multi-exposure clumping')
    warning('Need to make sure finding min P doesnt mess order up')
    d <- lapply(1:ncol(x@bx), function(exp_idx) as.data.table(x, exposure=exp_idx)) |> rbindlist()
    d <- d[ , .SD[which.min(px)], by='rsid']

  }

  # use object's LD data if available
  if (!all(is.na(x@correlation))) {

    # subset with LD data, maintain index, then order by exposure p-value
    d[, i := 1:.N] # i is now the index of `d`, into the LD matrix
    data.table::setorder(d, px)

    # setup data.table
    d[, c('index_snp','group') := list(FALSE, NA_integer_)]

    # loop
    flank = kb*1000
    row = 1
    grp = 1
    while(TRUE) {

      if(d$px[row]<p1 & is.na(d$group[row]) & !d$index_snp[row] & d$ld_info[row]) {

        # Found a significant hit, flag it
        d[row, index_snp := TRUE]

        # Take region and r2/p2 threshold it, ignore previously grouped SNPs in the region
        region <- d[bp >= d[row,bp]-flank & bp <= d[row,bp]+flank & ld_info, ]
        region[, rsq := x@correlation[d[row,i], region$i]^2, drop=TRUE]
        region <- region[!is.na(rsq) & rsq>r2 & px<p2 & is.na(group), ]

        # mutate the group back to d
        d[match(region$i, d$i), group := grp]

      }

      # find the next row
      tbc <- which(d$px<p1 & is.na(d$group) & !d$index_snp & d$ld_info)
      if(length(tbc)>0) {

        row <- tbc[1]
        grp <- grp + 1

      } else {

        # reset the order of d and put the grouping back
        data.table::setorder(d, i)
        x@index_snp <- d$index_snp
        x@group <- d$group
        break

      }
    }

  # if no LD matrix, try to use plink2
  } else {

    # remove allele info from rsids
    d[grepl("rs[0-9]+", rsid), rsid := sub(".*?(rs[0-9]+).*", "\\1", rsid)]

    # rename p for clump function
    data.table::setnames(d, "px", "p")

    # run clumping and rejoin result
    clumps <- clump(d, p1 = p1, p2 = p2, r2 = r2, kb = kb, plink2 = plink2, plink_ref = plink_ref)
    d[clumps, c("index_snp", "group") := list(i.index, i.group), on = "rsid"]

    # add back to object
    x@index_snp <- d$index
    x@group <- d$group
  }

  # return
  return(x)
}


#' @include class_gwas.R
method(as.data.table, MR) <- function(object, exposure = 1) {

  d <- data.table::data.table(
    rsid      = object@snps,
    chr       = object@chr,
    bp        = object@bp,
    ea        = object@ea,
    oa        = object@oa,
    bx        = object@bx[,exposure],
    bxse      = object@bxse[,exposure],
    px        = object@px[,exposure],
    by        = object@by,
    byse      = object@byse,
    py        = object@py,
    proxy_snp = if(length(object@proxy_snp)==0) { rep(NA_character_, length(object@snps)) } else { object@proxy_snp },
    index_snp = if(length(object@index_snp)==0) { rep(TRUE,          length(object@snps)) } else { object@index_snp },
    group     = if(length(object@group)==0)     { rep(NA_integer_,   length(object@snps)) } else { object@group },
    ld_info   = if(length(object@ld_info)==0)   { rep(NA_integer_,   length(object@snps)) } else { object@ld_info },
    exposure  = if(length(object@exposure)==0)  { rep(NA_character_, length(object@snps)) } else if(length(object@exposure)==length(object@snps)) { object@exposure } else { rep(object@exposure[exposure], length(object@snps)) },
    outcome   = if(length(object@outcome)==0)   { rep(NA_character_, length(object@snps)) } else if(length(object@outcome)==1)  { rep(object@outcome,  length(object@snps)) } else { object@outcome }
  )

}



#'
#' #' @title Plot clumps
#' #' @export
#' setGeneric("plot_clumps", function(object, gene_start=NULL, gene_end=NULL, gene_flanks=NULL, with_corr=FALSE) standardGeneric("plot_clumps"))
#' #' @rdname plot_clumps-methods
#' #' @aliases plot_clumps,MR,ANY-method
#' setMethod("plot_clumps", "MR", function(object, gene_start, gene_end, gene_flanks, with_corr) {
#'
#'   # convert object to dt
#'   plot_dat <- as.data.table(object)
#'
#'   # base plot
#'   p <- ggplot2::ggplot(mapping = ggplot2::aes(x=bp, y=-log10(px))) +
#'     ggplot2::geom_hline(yintercept=7.3, linetype="dotted", color="darkgrey")
#'
#'   # add gene region if requested
#'   if(!is.null(gene_start) && !is.null(gene_end) && !is.null(gene_flanks)) {
#'     p <- p +
#'       ggplot2::annotate(geom="rect", xmin=gene_start, xmax=gene_end, ymin=0, ymax=Inf, fill="blue", alpha=0.1) +
#'       ggplot2::scale_x_continuous(labels=scales::label_number(scale_cut=scales::cut_long_scale(),accuracy=1), limits=c(gene_start-gene_flanks, gene_end+gene_flanks))
#'   } else {
#'     p <- p +
#'       ggplot2::scale_x_continuous(labels=scales::label_number(scale_cut=scales::cut_long_scale(),accuracy=1))
#'   }
#'
#'   # add the points
#'   p <- p +
#'     ggplot2::geom_point(data=plot_dat, ggplot2::aes(color=as.factor(group))) +
#'     ggplot2::geom_point(data=plot_dat[index_snp==TRUE,], fill="red", shape=23, size=3) +
#'     ggplot2::labs(x="Chromosome position", y=paste0("-log10P - ", object@exposure), color="Clumps") +
#'     ggrepel::geom_label_repel(data=plot_dat[index_snp==TRUE,], ggplot2::aes(label=rsid))
#'
#'   # ?with outcome correlation plots
#'   if(with_corr) {
#'
#'     # clump names
#'     cn   <- plot_dat[index_snp==TRUE, list(rsid,group)]
#'     labs <- c(cn$rsid,NA_character_) |> `names<-`(c(cn$group, 'NA'))
#'
#'     # the plots
#'     p_corr <- ggplot2::ggplot() +
#'       ggplot2::geom_point(data=plot_dat, mapping=ggplot2::aes(x=bx, y=by, color=as.factor(group))) +
#'       ggplot2::geom_smooth(data=plot_dat, method="lm", formula= y ~ x - 1, mapping=ggplot2::aes(x=bx, y=by, weight=1/(byse^2)), color="red") +
#'       ggplot2::geom_point(data=plot_dat[index_snp==TRUE,], mapping=ggplot2::aes(x=bx, y=by), fill="red", shape=23, size=3) +
#'       ggplot2::theme_light() +
#'       ggplot2::facet_wrap(~group, ncol=4, scales="free", labeller = ggplot2::labeller(group=labs)) +
#'       ggplot2::labs(x=paste0("\u03B2 ", object@exposure), y=paste0("\u03B2 ", object@outcome), color="Clumps")
#'
#'     # combine under
#'     p <- ggpubr::ggarrange(p, p_corr, ncol=1)
#'
#'   }
#'
#'   # return the plot
#'   return(p)
#' })
#'
#'
#' #' @title Plot MR
#' #' @export
#' setGeneric("plot_mr", function(object, result=NULL, labels=TRUE, orientate=TRUE) standardGeneric("plot_mr"))
#' #' @rdname plot_mr-methods
#' #' @aliases plot_mr,MR,ANY-method
#' setMethod("plot_mr", "MR", function(object, result, labels, orientate) {
#'
#'   # checks
#'   stopifnot("results must be an `MRResult` class object, or the data.table output from mr_results_to_data_table(), or NULL" = is.null(result) || inherits(result, "MRResult") || inherits(result, "data.table"))
#'
#'   # convert MR object
#'   p_dat <- as.data.table(object)[index_snp==TRUE, ]
#'
#'   # plot orientated if requested
#'   if(orientate) { p_dat[bx<0, c('bx','by') := list(bx*-1, by*-1)] }
#'
#'   # the plot
#'   p <- ggplot2::ggplot(data    = p_dat,
#'                        mapping = ggplot2::aes(x=bx, y=by)) +
#'     ggplot2::geom_errorbar( mapping = ggplot2::aes(ymin=by-byse, ymax=by+byse), width=0, color="grey") +
#'     ggplot2::geom_errorbarh(mapping = ggplot2::aes(xmin=bx-bxse, xmax=bx+bxse), height=0, color="grey") +
#'     ggplot2::geom_point() +
#'     ggplot2::theme_linedraw() +
#'     ggplot2::labs(colour = "MR method",
#'                   x      = paste("SNP effect on", p_dat[1,exposure]),
#'                   y      = paste("SNP effect on", p_dat[1,outcome]))
#'
#'   # set x min at zero if orientated
#'   if(orientate) {
#'     p <- p + ggplot2::expand_limits(x=0, y=0)
#'   }
#'
#'   # plot labels
#'   if(labels) {
#'     p <- p +
#'       ggrepel::geom_label_repel(mapping = ggplot2::aes(label=rsid))
#'   }
#'
#'   # the MR result data / the line
#'   if(!is.null(result)) {
#'
#'     # if a simple MRResult object, convert to data table; otherwise assume multiple rows in a data.table
#'     if(inherits(result, "MRResult")) {
#'       l_dat <- mr_results_to_data_table(result)
#'     } else {
#'       l_dat <- result
#'     }
#'
#'     # colours for each test
#'     color_map = c("mr_egger"="#1f78b4",
#'                   "mr_mvegger" = "#a6cee3",
#'                   "mr_weighted_mode"="#6a3d9a",
#'                   "mr_mvweighted_mode"="#cab2d6",
#'                   "mr_weighted_median"="#33a02c",
#'                   "mr_mvweighted_median"="#b2df8a",
#'                   "mr_ivw"="#e31a1c",
#'                   "mr_mvivw" = "#fb9a99",
#'                   "mr_pcgmm"="#ff7f00",
#'                   "mr_mvpcgmm"="#fdbf6f")
#'
#'     # add the lines
#'     p <- p +
#'       ggplot2::geom_abline(data    = l_dat,
#'                            mapping = ggplot2::aes(intercept=Intercept, slope=Estimate, color=Method), show.legend=TRUE) +
#'       ggplot2::scale_colour_manual(values=color_map) +
#'       ggplot2::theme(legend.position="top", legend.direction="vertical") +
#'       ggplot2::guides(colour=ggplot2::guide_legend(ncol=2))
#'
#'     # simple IVW line through all the points
#'   } else {
#'
#'     p <- p +
#'       #geom_smooth(method="lm", mapping = aes(weight = 1/(byse^2)), color="red")
#'       geom_smooth(method="lm", formula= y ~ x - 1, mapping=aes(x=bx, y=by, weight=1/(byse^2)), color="red")
#'
#'   }
#'
#'   return(p)
#' })
#'
#'
#'
#' setGeneric("get_proxies", function(exposure, outcome, correlation, proxy_r2, proxy_kb, proxy_window, proxy_eaf) standardGeneric("get_proxies"))
#' setMethod("get_proxies", c("data.table","data.table","matrix"), function(exposure, outcome, correlation, proxy_r2, proxy_kb, proxy_window, proxy_eaf) {
#'
#'   # testing
#'   # exposure = to_TSMR(GWAS(GTEx_v8_GLP1R_ALL[trait %in% tissues[1] & bp>=bp_start-flank & bp<=bp_end+flank, ]), 'exposure')
#'   # outcome  = to_TSMR(GWAS(UKBB_GLP1R_HBA1C[bp>=bp_start-flank & bp<=bp_end+flank, ]), 'outcome')
#'   # correlation = UKBB_GLP1R_LDMAT
#'   # proxy_r2=0.8
#'   # proxy_kb=250
#'   # proxy_window=15
#'   # proxy_eaf=0.01
#'
#'   message("Searching for outcome proxies")
#'
#'   # which exposure variants don't have a partner?
#'   exposure[, rsid := sub("^(rs[0-9]+|[0-9XY]+:[0-9]+).*", "\\1", SNP)]
#'   outcome[ , rsid := sub("^(rs[0-9]+|[0-9XY]+:[0-9]+).*", "\\1", SNP)]
#'   snp_miss <- exposure[!rsid %in% outcome$rsid, ]
#'
#'   # which variants are available in the LD reference matrix?
#'   snp_ld <- data.table::data.table(
#'     id   = colnames(correlation),
#'     rsid = sub("(.*)_[ACTG]+_[ACTG]+$","\\1",colnames(correlation)),
#'     oa   = sub(".*_([ACTG]+)_[ACTG]+$","\\1",colnames(correlation)),
#'     ea   = sub(".*_[ACTG]+_([ACTG]+)$","\\1",colnames(correlation))
#'   )
#'
#'   # join the variants with LD info to the missing outcome SNPs that we want proxies for
#'   snp_miss[snp_ld, ld_id := i.id, on=c('rsid'='rsid','effect_allele.exposure'='ea','other_allele.exposure'='oa')]
#'   snp_miss[snp_ld, ld_id := i.id, on=c('rsid'='rsid','effect_allele.exposure'='oa','other_allele.exposure'='ea')]
#'   snp_miss <- snp_miss[!is.na(ld_id), ]
#'
#'   # join the outcome variants withthe LD info and harmonise to the LD panel
#'   outcome[snp_ld, c('proxy_flip','ld_id') := list(FALSE, i.id), on=c('rsid'='rsid','effect_allele.outcome'='ea','other_allele.outcome'='oa')]
#'   outcome[snp_ld, c('proxy_flip','ld_id') := list(TRUE,i.id), on=c('rsid'='rsid','effect_allele.outcome'='oa','other_allele.outcome'='ea')]
#'   outcome <- outcome[!is.na(ld_id), ]
#'   outcome[, ea_store := effect_allele.outcome]
#'   outcome[proxy_flip==TRUE, c('effect_allele.outcome','other_allele.outcome') := list(other_allele.outcome, ea_store)]
#'   outcome[proxy_flip==TRUE, c('eaf.outcome','beta.outcome') := list(1-eaf.outcome, beta.outcome*-1)]
#'   outcome[, c('ea_store','proxy_flip') := NULL]
#'
#'   # ensure at least 1 row to loop over (as we add the proxy columns)
#'   if(nrow(snp_miss)==0) snp_miss <- snp_miss[, lapply(.SD, function(x) NA)]
#'
#'   # assess each exposure SNP without an outcome SNP pair (cycle the rows; bp0=row BP, ld_id0=row ID)
#'   proxies <- mapply(function(bp0, ld_id0, proxy_kb, proxy_window, proxy_r2, outcome) {
#'
#'     # define candidate outcome proxies (distance, maf, snp-window; exclude self comparison dist==0)
#'     valid <- outcome[pos.outcome >= bp0-proxy_kb*1000 & pos.outcome <= bp0+proxy_kb*1000, ][, dist := abs(pos.outcome-bp0)][dist!=0, ][eaf.outcome>0.01 & eaf.outcome<0.99, ]
#'     valid <- head(valid[order(dist, pval.outcome), ], n=proxy_window)
#'
#'     # get the R2 values for the outcome SNP candidates vs the exposure SNP with a missing partner
#'     r2vals <- correlation[valid$ld_id, valid$ld_id0, drop=TRUE]^2
#'     valid[, c('proxy_snp.outcome','target_snp.outcome','proxy_r2.outcome') := list(ld_id, ld_id0, r2vals)]
#'     valid[, rsid := sub("^(rs[0-9]+|[0-9XY]+:[0-9]+)_.*","\\1", ld_id0)]
#'     valid[, SNP  := paste0(c(rsid,sort(c(other_allele.outcome,effect_allele.outcome))), collapse="_"), by=seq_len(nrow(valid))] # sorted SNP id - see load_GWAS
#'     valid[, effect_allele.outcome := sub("^(?:rs[0-9]+|[0-9XY]+:[0-9]+)_([ACTG]+|[DI])_.*","\\1", ld_id0)]
#'     valid[, other_allele.outcome  := sub("^(?:rs[0-9]+|[0-9XY]+:[0-9]+)_(?:[ACTG]+|[DI])_([ACTG]+|[DI]).*","\\1", ld_id0)]
#'     valid[, c('ld_id','dist','rsid') := NULL]
#'
#'     # pick the proxy with the largest R2 which reaches threshold
#'     valid <- valid[proxy_r2.outcome > proxy_r2, ]
#'     valid <- head(valid[order(-proxy_r2.outcome), ], n=1)
#'
#'     # return row
#'     return(valid)
#'
#'   }, bp0=snp_miss$pos.exposure, ld_id0=snp_miss$ld_id, MoreArgs=list(proxy_kb=proxy_kb,
#'                                                                      proxy_window=proxy_window,
#'                                                                      proxy_r2=proxy_r2,
#'                                                                      outcome=outcome), SIMPLIFY=FALSE)
#'
#'   # clean and combine
#'   exposure[, rsid := NULL]
#'   outcome[ , rsid := NULL]
#'   outcome[, ld_id := NULL]
#'   outcome <- data.table::rbindlist(c(list(outcome), proxies), fill=TRUE)
#'
#'   # return
#'   return(outcome)
#' })
#'
#'
#'
#'

#'
#' #' @title Reset index SNP
#' #' @export
#' setGeneric("reset_index_snp", function(x) standardGeneric("reset_index_snp"))
#' #' @rdname reset_index_snp-methods
#' #' @aliases reset_index_snp,MR,ANY-method
#' setMethod("reset_index_snp", "MR", function(x) {
#'   x@index_snp <- rep(TRUE, length(x@index_snp))
#'   validObject(x)
#'   return(x)
#' })
#'
#' #' @title Prune LD matrix
#' #' @export
#' setGeneric("prune_ld", function(x, r2_thresh=0.95, seed=2024) standardGeneric("prune_ld"))
#' #' @rdname prune_ld-methods
#' #' @aliases prune_ld,MR,ANY-method
#' setMethod("prune_ld", "MR", function(x, r2_thresh, seed) {
#'
#'   # checks
#'   stopifnot("LD matrix is not set" = !all(is.na(x@correlation)))
#'
#'   # just set index SNPs as those with LD info if no pruning
#'   if(r2_thresh==1) {
#'     x@index_snp <- x@ld_info
#'     return(x)
#'   }
#'
#'   # adapted from https://wellcomeopenresearch.org/articles/8-449
#'   set.seed(seed)
#'
#'   # only consider upper triangle of correlations
#'   rho.upper = x@correlation
#'   rho.upper[lower.tri(x@correlation, diag=TRUE)] <- 0
#'
#'   # indices of cells over the threshold
#'   over <- which(abs(rho.upper) > sqrt(r2_thresh), arr.ind=TRUE) |> as.data.frame()
#'
#'   # store the indices to omit
#'   omit <- c()
#'
#'   # cycle down the indices
#'   while (nrow(over)>0) {
#'
#'     # pick a row in over
#'     i <- sample(1:nrow(over), 1)
#'
#'     # randomly pick whether going to exclude matrix row or col
#'     j <- sample(1:2, 1)
#'
#'     # the matrix index to omit
#'     omit <- c(omit, over[i,j])
#'     over_idx <- which(over[,j]==over[i,j])
#'
#'     # filter out these row (or col) indices from the table
#'     over <- over[-over_idx, ]
#'   }
#'
#'   # just set index SNPs as those with LD info if no pruning
#'   x@index_snp <- x@ld_info
#'   x@index_snp[omit] <- FALSE
#'   return(x)
#' })
#'
#' setGeneric("mr_plotting_points", function(x, id_idx=1, orientate=TRUE) standardGeneric("mr_plotting_points"))
#' setMethod("mr_plotting_points", "MR", function(x, id_idx=1, orientate=TRUE) {
#'   points_df <- data.frame(
#'     snp = x@snps,
#'     x   = x@bx[,id_idx],
#'     y   = x@by,
#'     xse = x@bxse[,id_idx],
#'     yse = x@byse
#'   )
#'   if(orientate) {
#'     neg_x <- points_df$x < 0
#'     points_df$x[neg_x] <- points_df$x[neg_x] * -1
#'     points_df$y[neg_x] <- points_df$y[neg_x] * -1
#'   }
#'   return(points_df)
#' })

