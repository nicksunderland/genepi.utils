#' @title Collider bias result object
#' @description
#' An S3 class which is the standard output from various collider bias analyses.
#' @param method a string, the method used
#' @param ip a numeric, the initialising p-value
#' @param pi0 a numeric, the initialising Pi value for the Slope-hunter method
#' @param sxy1 a numeric, the initialising SNPip Cov value for the Slope-hunter method
#' @param b a numeric, the estimated correction factor
#' @param bse a numeric, the standard error of the correction factor
#' @param pi a numeric, the final Pi value, proportion of SNPs in the Gi SNP cluster
#' @param fit a data.table, the raw points data for plotting table columns `c("SNP_incidence","BETA_incidence","BETA_progression","CLUSTER")`
#' @param entropy a numeric, the entropy - output from the Slope-hunter method
#' @return an S3 object of class ColliderBiasResult
#' @export
#'
ColliderBiasResult <- function(method = NA_character_,
                               ip     = NA_real_,
                               pi0    = NA_real_,
                               pi     = NA_real_,
                               sxy1   = NA_real_,
                               b      = NA_real_,
                               intercept=NA_real_,
                               bse    = NA_real_,
                               entropy= NA_real_,
                               fit    = data.table::data.table()) {
  s3_struct <- structure(
    .Data = list(
      method = method,
      ip     = ip,
      pi0    = pi0,
      sxy1   = sxy1,
      b      = b,
      intercept=intercept,
      bse    = bse,
      entropy= entropy,
      pi     = pi,
      fit    = fit
    ),
    class = "ColliderBiasResult"
  )
}

#' @title ColliderBiasResult to data.table
#' @param results a ColliderBiasResult object, or list of ColliderBiasResult objects
#' @export
collider_results_to_data_table <- function(results) UseMethod("collider_results_to_data_table")
#' @rdname collider_results_to_data_table
#' @export
collider_results_to_data_table.ColliderBiasResult <- function(results) {
  # convert to list if singular
  collider_results_to_data_table(list(results))
}
#' @rdname collider_results_to_data_table
#' @export
collider_results_to_data_table.list <- function(results) {

  # if not all class ColliderBiasResult, error
  if (!is.list(results) || !all(sapply(results, inherits, "ColliderBiasResult"))) {

    stop("Input must be a ColliderBiasResult object or a list of ColliderBiasResult objects.")

  }

  # convert
  dt <- lapply(results, function(res) {

    # remove the large data.table as don't want this
    res$fit <- NULL

    # remove class, make just a list
    class(res) <- "list"

    # return as a data.table the rest of the list elements
    return(data.table::as.data.table(res))

  }) |>
    # rbind
    data.table::rbindlist()

  # return the full data.table
  return(dt)
}


#' @title Slope-hunter collider bias method
#' @description
#' The Slope-hunter method has it's own R package and does not need to be re-implemented.
#' I have re-written it here primarily for my own learning and to implement new plotting functions,
#' such as the iteration plotting output. \cr
#' Run the SlopeHunter expectation-maximization method to estimate the bias adjustment
#' factor. The algorithm uses model-based clustering and proposes that the distributions
#' of the incidence(GI) and progression (GP) BETAs can be written like so:  \cr
#' \cr
#' \eqn{
#' \binom{\beta_{GI}}{\beta'_{GI}}
#' \sim
#' \pi_{1}N
#' \begin{pmatrix}
#' \underline{0},
#' \begin{bmatrix}
#' \sigma_{I}^{2} & b_{1}\sigma_{I}^{2} \\
#' b_1\sigma_{I}^{2} & b_{1}^{2}\sigma_{I}^{2} \\
#' \end{bmatrix}
#' \end{pmatrix}
#' +
#' \pi_{2}N
#' \begin{pmatrix}
#' \underline{0},
#' \begin{bmatrix}
#' \sigma_{I}^{2} & b_{1}\sigma_{I}^{2} + \sigma_{IP} \\
#' b_1\sigma_{I}^{2} + \sigma_{IP} & b_{1}^{2}\sigma_{I}^{2} + \sigma_{P}^{2} + 2b_{1}\sigma_{IP} \\
#' \end{bmatrix}
#' \end{pmatrix}
#' +
#' \pi_{3}
#' \begin{pmatrix}
#' \eta_{0} \\
#' N(0, \sigma_{P}^{2})
#' \end{pmatrix}
#' +
#' \pi_{4}
#' \begin{pmatrix}
#' \eta_{0} \\
#' \eta_{0}
#' \end{pmatrix}
#' } \cr
#' \itemize{
#'   \item Cluster 1: SNPs that cause incidence but not progression
#'   \item Cluster 2: SNPs that cause incidence and progression
#'   \item Cluster 3: SNPs that cause progression but not incidence
#'   \item Cluster 4: SNPs that cause neither incidence or progression
#'   \item The values \eqn{\pi_{1}, \pi_{2}, \pi_{3}, \pi_{4}} are the probabilities that a SNP
#'   belongs to the respective clusters.
#' }
#' The first thing to note is that we can filter out SNPs from cluster
#' 3 and 4 by only including SNPs with a significant association (P-value)
#' with incidence - this is the `ip` parameter. The problem is then reduced
#' to finding two clusters that best fit distributions 1 and 2 of the equation
#' above. The SlopeHunter EM algorithm iteratively determines which SNPs
#' belong to each distribution (probabilistically) and once complete the
#' adjustment factor (slope gradient) can be determined from group 1, i.e. only
#' those SNPs thought to solely cause incidence. \cr
#' This function's code is adapted from the SlopeHunter R package and if using
#' this method the SlopeHunter package should be referenced as the original
#' source (see references).
#' @references https://github.com/Osmahmoud/SlopeHunter/blob/b0686fe0c7a4e6d056d902765ba9cd0a0b35ad5c/R/SH-utils.R#L53
#' @references https://doi.org/10.1038/s41467-022-28119-9
#' @param gwas_i a data.frame like object, clumped incidence GWAS summary statistics in `standardise_gwas()` format
#' @param gwas_p a data.frame like object, progression GWAS summary statistics in `standardise_gwas()` format
#' @inheritParams harmonise
#' @param ip a numeric, range 0-1, the initial p-value by which to filter the incidence GWAS
#' @param pi0 a numeric, range 0-1, the initial weight / percentage of SNPs that affect incidence only
#' @param sxy1 a numeric, the initial (guess) of the covariance between incidence and progression BETAs
#' @param bootstraps a whole-number numeric, the number of bootstrap samples to estimate the SE of the adjustment factor
#' (can be zero, in which case no SE is returned)
#' @param ... parameter sink, additional ignored parameters
#' @param seed an integer, seed for reproducibility
#' @return an S3 ColliderBiasResult object
#' @export
#'
slopehunter <- function(gwas_i,
                        gwas_p,
                        merge      = c("CHR"="CHR","BP"="BP"),
                        ip         = 0.001,
                        pi0        = 0.6,
                        sxy1       = 1e-5,
                        bootstraps = 100,
                        seed       = 2023,
                        ...) {

  # silence RMD checks
  P = keep = NULL

  # as data.tables
  gwas_i <- data.table::as.data.table(gwas_i)
  gwas_p <- data.table::as.data.table(gwas_p)

  # filter incidence gwas by threshold
  gwas_i <- gwas_i[P <= ip, ]
  if(nrow(gwas_i) == 0) {

    # exit if there are no significant incidence variants at this threshold
    warning(paste0("No variants remaining after thresholding incidence SNPs at `ip`=", ip, "\n"))
    return(ColliderBiasResult(method="slopehunter", ip=ip, pi0=pi0, sxy1=sxy1))

  }

  # harmonise and remove invalid alleles and palindromic SNPs
  h <- harmonise(gwas_i, gwas_p, gwas1_trait="incidence", gwas2_trait="progression", merge=merge)
  h <- h[keep==TRUE, ]

  # if there is data after thresholding and harmonising and run SH method
  if(nrow(h) > 0) {

    # run the EM algorithm
    result <- shclust(h, pi0=pi0, sxy1=sxy1, collect_iters=FALSE)

  } else {

    # exit if there are no significant incidence variants at this threshold
    warning(paste0("No variants remaining after merging with progression SNPs. Check
                    that incidence SNPs with P<threshold exist with in the progression GWAS\n"))
    return(ColliderBiasResult(method="slopehunter", ip=ip, pi0=pi0, sxy1=sxy1))

  }

  # run lots more times bootstrapping to estimate SE
  if (bootstraps > 0){
    b.boots = vector("numeric", bootstraps)
    for(i in 1:bootstraps){
      if(i%%10==0) cat("Bootstrapping... ", i, "of", bootstraps, "\n")
      set.seed(seed+i)
      boot_idxs   <- sample(1:nrow(h), size = nrow(h), replace = TRUE)
      h_boot      <- h[boot_idxs,]
      result_boot <- shclust(h_boot, pi0, sxy1, collect_iters=FALSE)
      b.boots[i]  <- result_boot$b
    }

    # If any of the model fits for the bootstrap samples generated NA
    if(any(is.na(b.boots))){
      b.boots <- b.boots[!is.na(b.boots)]
      warning(paste("Only", length(b.boots), "bootstrap samples - out of", bootstraps,
                    "- produced converged models used for estimating the standard error." ))
    }

    # calculate bse
    result$bse  = sd(abs(b.boots))
    result$bse_LB = result$b - 1.96*result$bse
    result$bse_UB = result$b + 1.96*result$bse
  }

  # overall result to return
  res <- ColliderBiasResult(method  = "slopehunter",
                            ip      = ip,
                            pi0     = pi0,
                            sxy1    = sxy1,
                            b       = result$b,
                            intercept=0,
                            bse     = result$bse,
                            entropy = result$entropy,
                            pi      = result$pi,
                            fit     = result$fit)
  return(res)
}


#' @title The SlopeHunter clustering algorithm
#' @inheritParams slopehunter
#' @param collect_iters a logical, whether to collect data from each EM algorithm iteration (for plotting)
#' @return a list with the results of the SH analysis run
#' @importFrom stats sd cov
#' @importFrom mclust dmvnorm
#' @noRd
shclust <- function(d, pi0, sxy1, collect_iters=FALSE) {

  f0 = f1 = pt = CLUSTER = keep = NULL

  d <- data.table::copy(d)

  # set the initial sdev and cov
  sx0 = sx1 = stats::sd(d$BETA_incidence)
  sy0 = sy1 = stats::sd(d$BETA_progression)
  dir0 = sign(stats::cov(d$BETA_incidence, d$BETA_progression))
  if (dir0==0) stop("All associations with at least either x or y are constant")

  # data for contour plotting
  iter_data <- list()
  iter_bs <- list()

  # convergence criterion
  loglkl_ck = 0

  # EM algorithm
  for(iter in 1:50000){
    # The E step:
    sxy0 = sx0*sy0*dir0*0.95 # enforce perfect correlation for distribution 1 (SNP -> incidence only)
    sigma0 = matrix(c(sx0^2,sxy0,sxy0,sy0^2), 2, 2) # covariance matrix for distribution 1 (SNP -> incidence only)
    sigma1 = matrix(c(sx1^2,sxy1,sxy1,sy1^2), 2, 2) # covariance matrix for distribution 2 (SNP -> incidence & progression)

    # get the probabilities that the SNP belongs to distribution 1 (f0) or 2 (f1)
    # 1st component / distribution
    d[, f0 := mclust::dmvnorm(d[, c("BETA_incidence", "BETA_progression")], c(0,0), sigma0)]
    d[, f0 := ifelse(f0<9.99999999999999e-301, 9.99999999999999e-301, f0)]

    # 2nd component
    d[, f1 := mclust::dmvnorm(d[, c("BETA_incidence", "BETA_progression")], c(0,0), sigma1)]
    d[, f1 := ifelse(f1<9.99999999999999e-301, 9.99999999999999e-301, f1)]

    # proportional contribution of density of f0 (for every point) to the total mixture
    d[, pt := pi0*f0/(pi0*f0+(1-pi0)*f1)]
    d[, pt := ifelse(pt>0.9999999, 0.9999999, pt)]

    # define the clusters
    d[, CLUSTER := factor(ifelse(pt >= 0.5, "Hunted", "Pleiotropic"), levels = c("NA", "Pleiotropic","Hunted"))] # hunted last level so plotted on top

    if(collect_iters & (iter %in% 1:20 | iter%%10==0)) {
      tmp <- d[, c("BETA_incidence", "BETA_progression", "pt", "CLUSTER")]
      iter_data <- c(iter_data, list(tmp))
      iter_bs <- c(iter_bs, list(dir0*sy0/sx0))
    }

    # loglik of the mixture model: pi0 * f0 + (1-p0) * f1
    loglkl = sum( log( pi0 * d$f0 + (1-pi0) * d$f1 ) )

    # The M step:
    pi0 = mean(d$pt) # update pi0; what is now the proportion of SNPs belonging to distribution 1 (i.e. with pt>=0.5)
    if (pi0<0.0001) pi0 = 0.0001
    if (pi0>0.9999) pi0 = 0.9999

    # update the standard deviations of distribution 1 now with it's new points (sx0 & sy0)
    sx0  = sqrt( sum(d$pt * (d$BETA_incidence^2))   / sum(d$pt) )
    sy0  = sqrt( sum(d$pt * (d$BETA_progression^2)) / sum(d$pt) )
    dir0 = sign( sum(d$pt * d$BETA_progression * d$BETA_incidence) / sum(d$pt) )
    if (dir0==0) {
      warning("Enforcing dir0=1")
      dir0=sample(c(1,-1), 1)   # avoid slope = 0 (horizontal line)
    }

    # update the standard deviations and cov of distribution 2 (sx1, sy1 & sxy1)
    sx1  = sqrt( sum((1-d$pt) * (d$BETA_incidence^2))   / (length(d$BETA_incidence)   - sum(d$pt)) )
    sy1  = sqrt( sum((1-d$pt) * (d$BETA_progression^2)) / (length(d$BETA_progression) - sum(d$pt)) )
    sxy1 = sum((1-d$pt) * d$BETA_incidence * d$BETA_progression) / (length(d$BETA_incidence) - sum(d$pt))
    if (abs(sxy1) > 0.75*sx1*sy1) {
      warning("Enforcing sxy1=sign(sxy1)*0.75*sx1*sy1")
      sxy1 = sign(sxy1)*0.75*sx1*sy1
    }

    # Check convergence
    if (iter%%10==0){
      if ((loglkl - loglkl_ck)/loglkl < 1e-10){
        break
      } else {
        loglkl_ck = loglkl
      }
    }
  }

  # Diagnosis
  if (iter == 50000) warning("The algorithm may not have converged.\n")

  # work out the results
  b       = dir0*sy0/sx0
  bse     = 0
  entropy = mean(d[CLUSTER=="Hunted", pt], na.rm=TRUE)

  # store the results
  out_cols <- c("SNP","BETA_incidence","BETA_progression","CLUSTER")
  data.table::setnames(d, c("SNP_incidence","BETA_incidence","BETA_progression","CLUSTER"), out_cols)
  d[, names(d)[!names(d)%in%out_cols] := NULL]
  result <- list(fit     = d,
                 b       = b,
                 bse     = bse,
                 bse_LB  = b - 1.96*bse,
                 bse_UB  = b + 1.96*bse,
                 entropy = entropy,
                 pi      = pi0,
                 iters   = iter_data,
                 iters_b = iter_bs)

  return(result)
}


#' @title Apply correction factor to GWAS data
#' @description
#' Adjust progression GWAS data for collider bias using a calculated correction
#' factor and correction factor standard error. The correction factor should be
#' calculated with one of the Slope-hunter, IVW-MR, or Dudbridge methods. \cr
#' This function first harmonises the datasets and removes invalid and palindromic
#' alleles as default
#' @param gwas_i a data.table, incidence GWAS
#' @param gwas_p a data.table, progression GWAS
#' @param b_correction_factor a numeric, the correction factor to apply
#' @param b_std_err a numeric, the correction factor standard error to apply
#' @param keep_palindromic a logical, whether to allow palindromic alleles to be adjusted and returned
#' @inheritParams harmonise
#' @return a data.table, progression GWAS with additional columns c("adjusted_beta","adjusted_se","adjusted_p")
#' @export
#' @importFrom stats pchisq
#'
apply_collider_correction <- function(gwas_i,
                                      gwas_p,
                                      b_correction_factor,
                                      b_std_err,
                                      merge = c("CHR"="CHR","BP"="BP"),
                                      keep_palindromic = FALSE) {

  # silence R CMD checks
  keep = palindromic = BETA_progression = BETA_incidence = SE_progression = SE_incidence = adjusted_beta = adjusted_se = adjusted_p = NULL

  # harmonise the GWASs;
  h <- harmonise(gwas_i, gwas_p, gwas1_trait="incidence", gwas2_trait="progression", merge=merge)

  # remove invalid alleles +/- palindromics
  if(keep_palindromic) {
    h <- h[keep==TRUE | palindromic==TRUE, ]
  } else {
    h <- h[keep==TRUE, ]
  }

  # adjust the beta, se and p
  h[, adjusted_beta := BETA_progression - (b_correction_factor * BETA_incidence)]
  h[, adjusted_se   := sqrt(
                              (SE_progression^2) +
                                ((SE_incidence^2) * (b_correction_factor^2)) +
                                ((SE_incidence^2) * (b_std_err^2)) +
                                ((SE_incidence^2) * (b_std_err^2))
                            )
    ]
  h[, adjusted_p := stats::pchisq( (adjusted_beta / adjusted_se)^2, 1, lower.tail=FALSE)]

  # return the harmonised adjusted data
  return(h)
}


#' @title Dudbridge CWLS collider bias method
#' @description
#' The Dudbridge corrected weighted least-squares method for assessing collider bias.
#' @inheritParams slopehunter
#' @inheritParams harmonise
#' @param ... parameter sink, additional ignored parameters
#' @return an S3 ColliderBiasResult object
#' @export
#' @importFrom MendelianRandomization mr_ivw
#'
dudbridge <- function(gwas_i,
                      gwas_p,
                      merge = c("CHR"="CHR","BP"="BP"),
                      ...) {

  # silence R CMD checks
  keep = NULL

  # as data.tables
  gwas_i <- data.table::as.data.table(gwas_i)
  gwas_p <- data.table::as.data.table(gwas_p)

  # harmonise and remove invalid alleles and palindromic SNPs
  h <- harmonise(gwas_i, gwas_p, gwas1_trait="incidence", gwas2_trait="progression", merge=merge)
  h <- h[keep==TRUE, ]

  # if there is data after thresholding and harmonising and run SH method
  if(nrow(h) == 0) {

    # exit if there are no significant incidence variants at this threshold
    warning(paste0("No variants remaining after merging with progression SNPs. Check
                    that incidence SNPs exist with in the progression GWAS\n"))
    return(ColliderBiasResult(method="dudbridge"))

  }

  # calculate the correction
  cwls_correction <- MendelianRandomization::mr_ivw(
    MendelianRandomization::mr_input(
      bx   = h$BETA_incidence,
      bxse = h$SE_incidence,
      by   = h$BETA_progression,
      byse = h$SE_progression
    )
  )

  # weighting
  weights <- 1 / h$SE_progression^2
  weighting <-  (sum(weights * h$BETA_incidence^2)) /
               ((sum(weights * h$BETA_incidence^2)) - (sum(weights * h$SE_incidence^2)))

  # slope
  cwls_estimated_slope <- cwls_correction$Estimate * weighting
  cwls_estimated_standard_error <- cwls_correction$StdError * weighting

  # overall result to return
  h[, "CLUSTER" := factor("NA")]
  out_cols <- c("SNP","BETA_incidence","BETA_progression", "CLUSTER")
  data.table::setnames(h, c("SNP_incidence","BETA_incidence","BETA_progression", "CLUSTER"), out_cols)
  h[, names(h)[!names(h)%in%out_cols] := NULL]
  res <- ColliderBiasResult(method  = "dudbridge",
                            b       = cwls_estimated_slope,
                            intercept=0,
                            bse     = cwls_estimated_standard_error,
                            fit     = h)
  return(res)
}

#' @title IVW MR collider bias method
#' @description
#' The Inverse variance weighted Mendelian Randomisation method for assessing collider bias.
#' @inheritParams slopehunter
#' @inheritParams harmonise
#' @param tsmr_method a string; a valid `TwoSampleMR::mr()` `method_list` parameter
#' @param ... parameter sink, additional ignored parameters
#' @return an S3 ColliderBiasResult object
#' @export
#' @importFrom TwoSampleMR format_data harmonise_data mr
#'
mr_collider_bias <- function(gwas_i,
                             gwas_p,
                             ip = 0.001,
                             merge = c("CHR"="CHR","BP"="BP"),
                             tsmr_method = "mr_ivw",
                             ...) {

  # silence RMD checks
  P = NULL

  # as data.tables
  gwas_i <- data.table::as.data.table(gwas_i)
  gwas_p <- data.table::as.data.table(gwas_p)

  # filter incidence gwas by threshold
  gwas_i <- gwas_i[P <= ip, ]
  if(nrow(gwas_i) == 0) {

    # exit if there are no significant incidence variants at this threshold
    warning(paste0("No variants remaining after thresholding incidence SNPs at `ip`=", ip, "\n"))
    return(ColliderBiasResult(method=tsmr_method, ip=ip))

  }

  # format for 2SMR
  tsmr_inc <- TwoSampleMR::format_data(gwas_i |> as.data.frame(),
                                       type = "exposure",
                                       snp_col = "SNP",
                                       beta_col = "BETA",
                                       se_col = "SE",
                                       eaf_col = "EAF",
                                       effect_allele_col = "EA",
                                       other_allele_col = "OA",
                                       pval_col = "P",
                                       chr_col = "CHR",
                                       pos_col = "POS")
  tsmr_pro <- TwoSampleMR::format_data(gwas_p |> as.data.frame(),
                                       type = "outcome",
                                       snp_col = "SNP",
                                       beta_col = "BETA",
                                       se_col = "SE",
                                       eaf_col = "EAF",
                                       effect_allele_col = "EA",
                                       other_allele_col = "OA",
                                       pval_col = "P",
                                       chr_col = "CHR",
                                       pos_col = "POS")
  tsmr_h <- TwoSampleMR::harmonise_data(tsmr_inc, tsmr_pro)

  # if there is data after harmonising and run method
  if(nrow(tsmr_h) == 0) {

    # exit if there are no incidence variants left after harmonising
    warning(paste0("No variants remaining after merging with progression SNPs. Check
                    that incidence SNPs exist with in the progression GWAS\n"))
    return(ColliderBiasResult(method=tsmr_method))

  } else {

    # run the MR
    mr_results <- TwoSampleMR::mr(tsmr_h, method_list = tsmr_method)

    # overall result to return
    tsmr_h <- data.table::as.data.table(tsmr_h)
    tsmr_h[, "CLUSTER" := factor("NA")]
    out_cols <- c("SNP","BETA_incidence","BETA_progression","CLUSTER")
    data.table::setnames(tsmr_h, c("SNP","beta.exposure","beta.outcome", "CLUSTER"), out_cols)
    tsmr_h[, names(tsmr_h)[!names(tsmr_h)%in%out_cols] := NULL]
    res <- ColliderBiasResult(method = tsmr_method,
                              b      = mr_results$b,
                              intercept = ifelse(is.null(mr_results$intercept), 0, mr_results$intercept),
                              bse    = mr_results$se,
                              ip     = ip,
                              fit    = tsmr_h)
    return(res)

  }
}


#' @title Run collider bias analysis
#' @description
#' Run a set of collider bias analyses with all combinations of input parameters.
#' Be careful with large numbers of parameters as this will quickly lead to very
#' large numbers of combinations. e.g. \cr
#' \itemize{
#'   \item methods = c("slopehunter")
#'   \item ip = c(0.05,0.001,0.00001)
#'   \item pi0 = c(0.6, 0.65, 0.7)
#'   \item sxy1 = c(1e-5, 1e-4, 1e-3)
#' }
#' ...will lead to 27 separate analyses. \cr
#' @inheritParams slopehunter
#' @param methods a character vector of collider bias method function names to run
#' @return a list of S3 ColliderBiasResult objects
#' @export
#'
analyse_collider_bias <- function(gwas_i,
                                  gwas_p,
                                  merge      = c("CHR"="CHR","BP"="BP"),
                                  methods    = c("slopehunter", "mr_collider_bias", "dudbridge"),
                                  tsmr_method= c("mr_ivw", "mr_egger_regression", "mr_simple_median", "mr_simple_mode", "mr_raps"),
                                  ip         = c(0.9),
                                  pi0        = c(0.6),
                                  sxy1       = c(1e-5),
                                  bootstraps = c(100),
                                  seed       = 2023) {

  # expand the parameters
  params <- lapply(methods, function(m) {
    if(m=="slopehunter") return(expand.grid(method="slopehunter", tsmr_method=NA_character_, ip=ip, pi0=pi0, sxy1=sxy1, bootstraps=bootstraps))
    if(m=="dudbridge") return(expand.grid(method="dudbridge", tsmr_method=NA_character_, ip=NA_real_, pi0=NA_real_, sxy1=NA_real_, bootstraps=NA_real_))
    if(m=="mr_collider_bias") return(expand.grid(method="mr_collider_bias", tsmr_method=tsmr_method, ip=ip, pi0=NA_real_, sxy1=NA_real_, bootstraps=NA_real_))
  }) |> data.table::rbindlist()

  # results from each combination of parameters
  results <- list()
  for(i in 1:nrow(params)) {

    res <- do.call(what = as.character(params$method[i]),
                   args = list(gwas_i = gwas_i,
                               gwas_p = gwas_p,
                               merge  = merge,
                               tsmr_method = as.character(params$tsmr_method[i]),
                               ip          = params$ip[i],
                               pi0         = params$pi0[i],
                               sxy1        = params$sxy1[i],
                               bootstraps  = params$bootstraps[i]))

    results <- c(results, list(res))
  }

  return(results)
}


#' @title Plot the Slope-hunter scatter plot
#' @param results a ColliderBiasResult object, or list of ColliderBiasResult objects, from `analyse_collider_bias()`,
#' `slopehunter()`, `mr_collider_bias()`, or `dudbridge()`
#' @return a ggplot2
#' @export
#' @import ggplot2
#'
plot_slope <- function(results) {

  # silence RCMD checks
  BETA_incidence = BETA_progression = CLUSTER = hjustvar = vjustvar = xpos = ypos = label = slope = intercept = ip = NULL

  # make list of one if a single ColliderBiasResult object
  if(inherits(results, "ColliderBiasResult")) {
    results <- list(results)
  }

  # pvalue levels
  p_levels <- unique(sort(sapply(results, function(res) res$ip)))

  # results for each pvalue
  dat <- lapply(results, function(res) res$fit[, c("ip", "method") := list(res$ip, res$method)]) |> data.table::rbindlist()
  dat[, ip := factor(ip, levels=p_levels)]

  # slope data
  b_dat <- lapply(results, function(res) data.table::data.table(slope=res$b,
                                                                intercept=res$intercept,
                                                                label = paste0("b = ", as.character(round(res$b, digits=2))),
                                                                hunted_snps_label = ifelse(sum(res$fit$CLUSTER=="Hunted", na.rm=T)>0, paste0("Hunted = ", sum(res$fit$CLUSTER=="Hunted", na.rm=T)), NA_character_),
                                                                xpos = Inf,
                                                                ypos = Inf,
                                                                hjustvar = 1.2,
                                                                vjustvar = 1.5,
                                                                ip = res$ip,
                                                                method = res$method)) |> data.table::rbindlist()
  b_dat[, ip := factor(ip, levels=p_levels)]

  # facet plot of different P value thresholds
  g <- ggplot2::ggplot(data=dat, ggplot2::aes(BETA_incidence, BETA_progression, color=CLUSTER, alpha=CLUSTER, size=CLUSTER)) +
    ggplot2::geom_point() +
    ggplot2::geom_abline(data=b_dat, ggplot2::aes(slope=slope, intercept=intercept), color = 'black') +
    ggplot2::geom_text(data=b_dat, ggplot2::aes(x=xpos,y=ypos,hjust=hjustvar+0.2,vjust=vjustvar,label=label), color="black", show.legend = FALSE) +
    ggplot2::geom_text(data=b_dat, ggplot2::aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar+2,label=hunted_snps_label), color="black", show.legend = FALSE) +
    ggplot2::scale_size_manual(values=c("Hunted"=2, "Pleiotropic"=1, "NA"=1)) +
    ggplot2::scale_alpha_manual(values=c("Hunted"=1, "Pleiotropic"=0.3, "NA"=1)) +
    ggplot2::scale_color_manual(values=c("Hunted"="#FC4E07", "Pleiotropic"="#00AFBB", "NA"="#00AFBB")) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "SlopeHunter analysis by p-value threshold (\u03BB)",
                  x = expression("\u03B2"[indicence]),
                  y = expression("\u03B2"[progression]),
                  color = "Cluster",
                  size = "Cluster",
                  alpha = "Cluster") +
    ggplot2::facet_grid(method ~ ip,
                        labeller = labeller(ip = ~ paste0("\u03BB = ", .x)),
                        scales = "free")

  return(g)
}


#' @title Create an animation of Slope-hunter iterations
#' @inheritParams slopehunter
#' @param x_lims x axis limits
#' @param y_lims y axis limits
#' @return a plot
#' @export
#' @importFrom gganimate transition_states
#'
plot_slopehunter_iters <- function(gwas_i,
                                   gwas_p,
                                   merge      = c("CHR"="CHR","BP"="BP"),
                                   ip         = 0.001,
                                   pi0        = 0.6,
                                   sxy1       = 1e-5,
                                   seed       = 2023,
                                   x_lims     = NULL,
                                   y_lims     = NULL) {

  # silence RMD checks
  P = CLUSTER = BETA_incidence = BETA_progression = pt = iter = keep = xpos = ypos = hjustvar = vjustvar = label = slope = intercept = NULL

  # as data.tables
  gwas_i <- data.table::as.data.table(gwas_i)
  gwas_p <- data.table::as.data.table(gwas_p)

  # filter incidence gwas by threshold
  gwas_i <- gwas_i[P <= ip, ]
  if(nrow(gwas_i) == 0) {

    # exit if there are no significant incidence variants at this threshold
    warning(paste0("No variants remaining after thresholding incidence SNPs at `ip`=", ip, "\n"))
    return(ColliderBiasResult(method="slopehunter", ip=ip, pi0=pi0, sxy1=sxy1))

  }

  # harmonise and remove invalid alleles and palindromic SNPs
  h <- harmonise(gwas_i, gwas_p, gwas1_trait="incidence", gwas2_trait="progression", merge=merge)
  h <- h[keep==TRUE, ]

  # if there is data after thresholding and harmonising and run SH method
  if(nrow(h) == 0) {

    # exit if there are no significant incidence variants at this threshold
    warning(paste0("No variants remaining after merging with progression SNPs. Check
                    that incidence SNPs with P<threshold exist with in the progression GWAS\n"))
    return(ColliderBiasResult(method="slopehunter", ip=ip, pi0=pi0, sxy1=sxy1))

  }

  # run the EM algorithm
  result <- shclust(h, pi0=pi0, sxy1=sxy1, collect_iters=TRUE)

  # make the plot data
  result_dat <- data.table::rbindlist(result$iters, idcol="iter")

  # make the slope data
  result_b <- data.table::data.table(slope = unlist(result$iters_b),
                                     intercept = 0,
                                     iter=1:length(result$iters_b),
                                     label = paste0("b = ", as.character(round(unlist(result$iters_b), digits=2))),
                                     xpos = Inf, ypos = Inf, hjustvar = 1.2, vjustvar = 1.4)

  # Make a ggplot, but add frame=year: one image per year
  p <- ggplot2::ggplot() +
    ggplot2::geom_text(data=result_b, mapping=ggplot2::aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=label)) +
    ggplot2::geom_abline(data=result_b, mapping=ggplot2::aes(slope=slope,intercept=intercept), color="darkred") +
    ggplot2::geom_point(data=result_dat, mapping=ggplot2::aes(x=BETA_incidence, y=BETA_progression, fill=pt, color=CLUSTER, size=CLUSTER), shape=21) +
    viridis::scale_fill_viridis(limits=c(0,1)) +
    ggplot2::scale_size_manual(values = c("Hunted"=1.5, "Pleiotropic"=0.5)) +
    ggplot2::scale_alpha_manual(values = c("Hunted"=1.0, "Pleiotropic"=0.1)) +
    ggplot2::scale_shape_manual(values = c("Hunted"=21, "Pleiotropic"=1)) +
    ggplot2::scale_color_manual(values = c("Hunted"="red", "Pleiotropic"="#44015410")) +
    ggplot2::theme_bw() +
    # gganimate specific bits:
    ggplot2::labs(title = 'Slope-hunter iteration: {closest_state}',
                  x = expression("\u03B2"[indicence]),
                  y = expression("\u03B2"[progression]),
                  color = "Cluster",
                  fill = "Likelihood Gi SNP") +
    ggplot2::guides(size = "none")

  # limits
  if(!is.null(x_lims)) {
    p <- p + ggplot2::xlim(x_lims)
  }
  if(!is.null(y_lims)) {
    p <- p + ggplot2::ylim(y_lims)
  }

  p <- p +
    gganimate::transition_states(iter, transition_length=0.1, state_length=1)

  return(p)
}


#' @title Plot the estimated correction factor variability
#' @param results a list of ColliderBiasResult objects, from `analyse_collider_bias()`
#' @return a ggplot2 plot
#' @export
#' @import ggplot2
#'
plot_correction_stability <- function(results) {

  # silence RMD checks
  ip = b = bse = method = NULL

  # the collider results
  plot_dat <- collider_results_to_data_table(results)
  plot_dat[, ip := ifelse(is.na(ip), 0, ip)]

  # labels
  log10P_label <- expression(paste( "-log"[10], plain(P), " `ip` threshold"))

  # b by p-values
  p <- ggplot2::ggplot(data=plot_dat,
                       mapping=ggplot2::aes(x=-log10(ip), y=b)) +
    ggplot2::geom_ribbon(mapping=ggplot2::aes(ymin=b-(1.96*bse), ymax=b+(1.96*bse), fill=method), alpha=0.3) +
    ggplot2::geom_point(mapping=ggplot2::aes(color=method)) +
    ggplot2::geom_line(mapping=ggplot2::aes(color=method)) +
    ggplot2::labs(title = 'Bias correction factor by p-value threshold',
                  x = log10P_label,
                  y = 'b value (correction slope)',
                  color = "Method",
                  fill = "95%CI")

  return(p)
}
