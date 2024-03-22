#' @title Run collider bias assessment
#' @param x an object of class MR
#' @param bias_method a character or character vector, one or more of c("dudbridge", "slopehunter", "mr_ivw", "mr_egger", "mr_weighted_median", "mr_weighted_mode")
#' @param r2 a numeric 0-1, r2 used for clumping - set all clumping params to NA to turn off
#' @param p1 a numeric 0-1, p1 used for clumping - set all clumping params to NA to turn off
#' @param kb an integer, kb used for clumping - set all clumping params to NA to turn off
#' @param plink2 a path, the plink2 binary
#' @param plink_ref a path, the reference genome pfile
#' @param ip a numeric 0-1, threshold for removing incidence variants
#' @param pi0 a numeric 0-1, proportion of SNPs in the incidence only cluster (see Slope-Hunter)
#' @param sxy1 a numeric, the covariance between incidence and progression Gip SNPs (see Slope-Hunter)
#' @param bootstraps an integer, number of bootstraps to estimate SE (see Slope-Hunter)
#' @export
collider_bias <- new_generic("collider_bias", "x", function(x,
                                                            bias_method= "dudbridge",
                                                            # clumping
                                                            r2         = 0.001,
                                                            p1         = 5e-8,
                                                            kb         = 250,
                                                            plink2     = genepi.utils::which_plink2(),
                                                            plink_ref  = genepi.utils::which_1000G_reference(build="GRCh37"),
                                                            # slopehunter
                                                            ip         = 0.001,
                                                            pi0        = 0.6,
                                                            sxy1       = 1e-5,
                                                            bootstraps = 100,
                                                            # dudbridge
                                                            weighted   = TRUE,
                                                            method     = "Simex",
                                                            B          = 1000,
                                                            lambda     = seq(0.25, 5, 0.25),
                                                            # cwls
                                                            # ip       = 0.001 - use the same as SH above
                                                            seed       = 2023
                                                            ) { S7_dispatch() })
#' @name collider_bias
method(collider_bias, MR) <- function(x,
                                      bias_method= "dudbridge",
                                      # clumping
                                      r2         = 0.001,
                                      p1         = 5e-8,
                                      kb         = 250,
                                      plink2     = genepi.utils::which_plink2(),
                                      plink_ref  = genepi.utils::which_1000G_reference(build="GRCh37"),
                                      # slopehunter
                                      ip         = 0.001,
                                      pi0        = 0.6,
                                      sxy1       = 1e-5,
                                      bootstraps = 100,
                                      # dudbridge
                                      weighted   = TRUE,
                                      method     = "Simex",
                                      B          = 1000,
                                      lambda     = seq(0.25, 5, 0.25),
                                      # cwls
                                      # ip       = 0.001 - use the same as SH above
                                      seed       = 2023
                                      ) {

  cat("Running collider bias assessment...\n")

  # params: check method; multivariable not supported
  bias_method <- match.arg(bias_method, choices = c("dudbridge", "slopehunter", "mr_ivw", "mr_egger", "mr_weighted_median", "mr_weighted_mode"), several.ok = TRUE)
  method <- match.arg(method, choices = c("Hedges-Olkin", "Simex"))
  if (is_multivariable(x)) {
    stop("collider bias assessment only supported between one incidence and one progression GWAS, MR object is multivariable (multiple exposures)")
  }

  # expand the parameters
  params <- expand.grid(bias_method= bias_method,
                        Exposure   = x@exposure,
                        Outcome    = x@outcome,
                        r2         = r2,
                        p1         = p1,
                        kb         = kb,
                        ip         = ip,
                        pi0        = pi0,
                        sxy1       = sxy1,
                        bootstraps = bootstraps,
                        weighted   = weighted,
                        method     = method,
                        B          = B,
                        stringsAsFactors = FALSE)

  # results from each combination of parameters
  results <- list()
  for(i in 1:nrow(params)) {
    cat("Active parameters:\n")
    print(params[i, ])

    # clump if parameters provided
    if (all(!is.na(c(r2, p1, kb)))) {
      x <- clump_mr(x, r2 = r2, p1 = p1, kb = kb, plink2 = plink2, plink_ref = plink_ref)
    } else if(any(!x@index_snp)) {
      warning("No clumping requested but the MR object already has altered index_snp field. If you have processed this prior then this is fine.")
    }

    # run the method
    res <- do.call(what = as.character(params$bias_method[i]),
                   args = c(list(x = x), as.list(params[i, ])))

    # bind
    results <- c(results, list(res))
  }

  # combine parameters and results
  results_dt <- mr_results_to_data_table(results)
  results_dt <- cbind(results_dt, params[, names(params)[!names(params) %in% c("Exposure", "Outcome")]])
  results_dt[, c("Method", "Correl") := NULL]

  # return
  return(results_dt)
}



#' @title Dudbridge collider bias method
#' @param ... parameter sink, additional ignored parameters
#' @return an object of class MRResult
#' @export
#' @importFrom indexevent indexevent
dudbridge <- new_generic("dudbridge", "x", function(x, weighted = TRUE, prune = NULL, method = "Simex", B = 1000, lambda = seq(0.25, 5, 0.25), seed = 2018, ...) { S7_dispatch() })
#' @name dudbridge
method(dudbridge, MR) <- function(x,
                                  weighted = TRUE,
                                  prune    = NULL,
                                  method   = "Simex",
                                  B        = 1000,
                                  lambda   = seq(0.25, 5, 0.25),
                                  seed     = 2018,
                                  ...) {

  # Outcome and Exposure params in here
  extra_args <- list(...)

  # to_MRInput filters by index_snp
  dat <- to_MRInput(x, corr = FALSE)

  # Dudbridge indexevent
  res <- indexevent::indexevent(xbeta    = dat@betaX,
                                xse      = dat@betaXse,
                                ybeta    = dat@betaY,
                                yse      = dat@betaYse,
                                weighted = weighted,
                                prune    = prune,
                                method   = method,
                                B        = B,
                                lambda   = lambda,
                                seed     = seed)

  # populate an MRResult object with the results
  res_obj <- MRResult()
  res_obj@Method    <- "dudbridge"
  res_obj@Estimate  <- res$b
  res_obj@StdError  <- res$b.se
  res_obj@Intercept <- 0
  res_obj@SNPs      <- length(dat@snps)
  res_obj@Method    <- "dudbridge"
  res_obj@Outcome   <- extra_args[["Outcome"]]
  res_obj@Exposure  <- extra_args[["Exposure"]]

  # return
  return(res_obj)
}


#' @title Slope-Hunter collider bias method
#' @param ... parameter sink, additional ignored parameters
#' @return an object of class MRResult
#' @export
#' @importFrom SlopeHunter hunt
slopehunter <- new_generic("slopehunter", "x", function(x, ip = 0.001, pi0 = 0.6, sxy1 = 1e-5, bootstraps = 100, seed = 777, ...) { S7_dispatch() })
#' @name slopehunter
method(slopehunter, MR) <- function(x,
                                    ip         = 0.001,
                                    pi0        = 0.6,
                                    sxy1       = 1e-5,
                                    bootstraps = 100,
                                    seed       = 777,
                                    ...) {

  # Outcome and Exposure params in here
  extra_args <- list(...)

  # to_MRInput filters by index_snp
  dat <- to_MRInput(x, corr = FALSE)
  dat_df <- data.frame(snps = dat@snps,
                       bx   = dat@betaX,
                       bxse = dat@betaXse,
                       by   = dat@betaY,
                       byse = dat@betaYse)

  # SlopeHunter hunt
  res <- SlopeHunter::hunt(dat = dat_df,
                           xbeta_col     = "bx",
                           xse_col       = "bxse",
                           ybeta_col     = "by",
                           yse_col       = "byse",
                           xp_thresh     = ip,
                           init_pi       = pi0,
                           init_sigmaIP  = sxy1,
                           Bootstrapping = TRUE,
                           M             = bootstraps,
                           seed          = seed,
                           Plot          = FALSE)

  # populate an MRResult object with the results
  res_obj <- MRResult()
  res_obj@Method          <- "slopehunter"
  res_obj@Estimate        <- res$b
  res_obj@StdError        <- res$bse
  res_obj@Intercept       <- 0
  res_obj@SNPs            <- nrow(res$Fit)
  res_obj@Hunted          <- sum(res$Fit$clusters == "Hunted", na.rm = TRUE)
  res_obj@SlopeHunter.Pi  <- res$pi
  res_obj@SlopeHunter.Ent <- res$entropy
  res_obj@Outcome         <- extra_args[["Outcome"]]
  res_obj@Exposure        <- extra_args[["Exposure"]]

  # return
  return(res_obj)
}

