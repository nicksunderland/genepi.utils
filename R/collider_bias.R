#' @title ColliderBias class
#' @description
#' An S3 class for collider bias analysis. This object can be passed around the
#' collider bias analysis tools, used for plotting, and for applying the chosen
#' adjustment to your GWAS summary statistics.
#' @param gwas_i a data.frame like object, incidence GWAS summary statistics in `standardise_gwas()` format
#' @param gwas_p a data.table like object, progression GWAS summary statistics in `standardise_gwas()` format
#' @param methods a character vector, one or more valid collider bias analysis methods `slopehunter`, `dudbridge`, `ivw_mr`
#' @param ip a numeric, range 0-1, the initial p-value by which to filter the incidence GWAS
#' @param pi0 a numeric, range 0-1, the initial weight / percentage of SNPs that affect incidence only
#' @param sxy1 a numeric, the initial (guess) of the covariance between incidence and progression BETAs
#' @param bootstraps a whole-number numeric, the number of bootstrap samples to estimate the SE of the adjustment factor (can be zero)
#' @return an S3 object of class ColliderBias
#' @importFrom data.table as.data.table
#' @importFrom TwoSampleMR format_data harmonise_data
#' @export
#'
ColliderBias <- function(gwas_i,
                         gwas_p,
                         methods    = c("slopehunter","dudbridge","ivw_mr"),
                         ip         = 0.001,
                         pi0        = 0.6,
                         sxy1       = 1e-5,
                         bootstraps = 100) {

  # checks
  stopifnot("Invalid `ColliderBias` parameters" = valid_params(methods, ip, pi0, sxy1, bootstraps))

  # format data - add cols used by TwoSampleMR::format_data
  if(!"phenotype" %in% names(gwas_i)) gwas_i[, "phenotype" := "incidence"]
  if(!"phenotype" %in% names(gwas_p)) gwas_p[, "phenotype" := "progression"]
  if(!"id" %in% names(gwas_i))        gwas_i[, "id" := "incidence"]
  if(!"id" %in% names(gwas_p))        gwas_p[, "id" := "progression"]

  # define the columns for TwoSampleMR::format_data
  cols <- list(snp_col="SNP",chr_col="CHR",pos_col="BP",effect_allele_col="EA",other_allele_col="OA",beta_col="BETA",se_col="SE",pval="P",eaf_col="EAF",phenotype_col="phenotype",id_col="id")
  i_cols = p_cols = cols
  if("N" %in% names(gwas_p)) i_cols <- c(i_cols, list(samplesize_col="N"))
  if("N" %in% names(gwas_p)) p_cols <- c(p_cols, list(samplesize_col="N"))

  # run the formatting
  gwas_i <- do.call(TwoSampleMR::format_data, c(list(dat=gwas_i, type="exposure"), i_cols))
  gwas_p <- do.call(TwoSampleMR::format_data, c(list(dat=gwas_p, type="outcome"), p_cols))

  # harmonise effects
  gwas_harmonised <- TwoSampleMR::harmonise_data(gwas_i, gwas_p) |> data.table::as.data.table()

  # package the data into a structure of class ColliderBias
  s3_struct <- structure(
    .Data = list(
      gwas_i     = gwas_i,
      gwas_p     = gwas_p,
      harmonised = gwas_harmonised,
      methods    = c("slopehunter","dudbridge","ivw_mr"),
      ip         = 0.001,
      pi0        = 0.6,
      sigma0     = 1e-5,
      bootstraps = 100,
      results    = list(),
      plots      = list()
    ),
    class = "ColliderBias"
  )

  return(s3_struct)
}


#' @title Run collider bias analysis
#' @description
#' Run all the analyses defined in the ColliderBias object `x` with all combinations
#' of the listed parameters. Be careful with large numbers of parameters as this will
#' quickly lead to very large numbers of combinations. e.g. \cr
#' \itemize{
#'   \item methods = c("slopehunter")
#'   \item ip = c(0.05,0.001,0.00001)
#'   \item pi0 = c(0.6, 0.65, 0.7)
#'   \item sxy1 = c(1e-5, 1e-4, 1e-3)
#' }
#' ...will lead to 27 separate analyses. \cr
#' @param x a ColliderBias object
#' @inheritParams ColliderBias
#' @return a ColliderBias object with updated parameter and results fields
#' @export
analyse <- function(x,methods,ip,pi0,sxy1) UseMethod("analyse")
#' @rdname analyse
#' @export
analyse.ColliderBias <- function(x,methods,ip,pi0,sxy1) {
  return(x)
}


#' @title SlopeHunter method
#' @description
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
#' @param x a ColliderBias object
#' @inheritParams ColliderBias
#' @param seed an integer, seed for reproducibility
#' @return a ColliderBias object with updated parameter and results fields
#' @export
slopehunter <- function(x,ip,pi0,sxy1,bootstraps,seed) UseMethod("slopehunter")
#' @rdname slopehunter
#' @export
slopehunter.ColliderBias <- function(x,
                                     ip         = 0.9,
                                     pi0        = 0.6,
                                     sxy1       = 1e-5,
                                     bootstraps = 100,
                                     seed       = 2023) {

  pval.exposure = NULL

  # filter by the incidence threshold
  d <- x$harmonised[pval.exposure < ip, ]

  # exit if there are no significant incidence variants at this threshold
  if(nrow(d)==0) {
    warning(paste0("No significant incidence SNPs at `ip` ", ip, "\n"))
    return(x)
  }

  # run the EM algorithm
  result <- shclust(d, pi0, sxy1, collect_iters=FALSE)

  # run lots more times bootstrapping to estimate SE
  if (bootstraps > 0){
    b.boots = vector("numeric", bootstraps)
    for(i in 1:bootstraps){
      cat("Bootstrap sample", i, "of", bootstraps, "samples ...\n")
      set.seed(seed+i)
      boot_idxs = sample(1:nrow(d), size = nrow(d), replace = TRUE)
      d_boot = d[boot_idxs,]
      result_boot = shclust(d_boot, pi0, sxy1, collect_iters=FALSE)
      b.boots[i] = result_boot$b
    }

    # If any of the model fits for the bootstrap samples generated NA
    if(any(is.na(b.boots))){
      b.boots <- b.boots[!is.na(b.boots)]
      warning(paste("Only", length(b.boots), "bootstrap samples - out of", bootstraps, "- produced converged models used for estimating the standard error." ))
    }

    # calculate bse
    result$bse  = sd(abs(b.boots))
    result$b_CI = c(result$b - 1.96*result$bse, result$b + 1.96*result$bse)
  }

  # print the results
  cat("SlopeHunter:  b =", result$b, " ( SE", result$bse, ")\n")

  # add the result to the collider object
  x$results[["slopehunter"]] <- result

  # return the updated object
  return(x)
}


#' @title The SlopeHunter clustering algorithm
#' @inheritParams slopehunter
#' @param collect_iters a logical, whether to collect data from each EM algorithm iteration (for plotting)
#' @return a list with the results of the SH analysis run
#' @importFrom stats sd cov
#' @importFrom mclust dmvnorm
#' @noRd
shclust <- function(d, pi0, sxy1, collect_iters=FALSE) {

  f0 = f1 = pt = clusters = NULL

  d <- data.table::copy(d)

  # set the initial sdev and cov
  sx0 = sx1 = stats::sd( d[["beta.exposure"]] )
  sy0 = sy1 = stats::sd( d[["beta.outcome"]] )
  dir0 = sign(  stats::cov(d[["beta.exposure"]], d[["beta.outcome"]]) )
  if (dir0==0) stop("All associations with at least either x or y are constant")

  # data for contour plotting
  iter_data <- list()
  b_iter_data <- list()

  # convergence criterion
  loglkl_ck = 0

  # EM algorithm
  for(iter in 1:50000){
    # The E step:
    sxy0 = sx0*sy0*dir0*0.95 # enforce perfect correlation for distribution 1 (SNP -> incidence only)
    sigma0 = matrix(c(sx0^2,sxy0,sxy0,sy0^2), 2, 2) # covariance matrix for distribution 1 (SNP -> incidence only)
    sigma1 = matrix(c(sx1^2,sxy1,sxy1,sy1^2), 2, 2) # covariance matrix for distribution 2 (SNP -> incidence & progression)

    # get the probabilities that the SNP belongs to distribution 1 (f0) or 2 (f1)
    cols <- c("beta.exposure", "beta.outcome")

    # 1st component / distribution
    d[, f0 := mclust::dmvnorm(d[, cols, with=FALSE], c(0,0), sigma0)]
    d[, f0 := ifelse(f0<1e-300, 1e-300, f0)]

    # 2nd component
    d[, f1 := mclust::dmvnorm(d[, cols, with=FALSE], c(0,0), sigma1)]
    d[, f1 := ifelse(f1<1e-300, 1e-300, f1)]

    # proportional contribution of density of f0 (for every point) to the total mixture
    d[, pt := pi0*f0/(pi0*f0+(1-pi0)*f1)]
    d[, pt := ifelse(pt>0.9999999, 0.9999999, pt)]

    # define the clusters
    d[, clusters := factor(ifelse(pt >= 0.5, "Hunted", "Pleiotropic"))]

    if(collect_iters & (iter %in% 1:20 | iter%%10==0)) {
      tmp <- d[, c("beta.exposure", "beta.outcome", "pt", "clusters")]
      iter_data <- c(iter_data, list(tmp))
    }

    # loglik of the mixture model: pi0 * f0 + (1-p0) * f1
    loglkl = sum(log( pi0*d[["f0"]]+ (1-pi0)*d[["f1"]] ))

    # The M step:
    pi0 = mean(d[['pt']]) # update pi0; what is now the proportion of SNPs belonging to distribution 1 (i.e. with pt>=0.5)
    if (pi0<0.0001) pi0 = 0.0001
    if (pi0>0.9999) pi0 = 0.9999

    # update the standard deviations of distribution 1 now with it's new points (sx0 & sy0)
    sx0  = sqrt(sum(d[['pt']]*(d[["beta.exposure"]]^2))/sum(d[['pt']]))
    sy0  = sqrt(sum(d[['pt']]*(d[["beta.outcome"]]^2))/sum(d[['pt']]))
    dir0 = sign(sum(d[['pt']]*d[["beta.outcome"]]*d[["beta.exposure"]])/sum(d[['pt']]))
    if (dir0==0) dir0=sample(c(1,-1), 1)   # avoid slope = 0 (horizontal line)

    # update the standard deviations and cov of distribution 2 (sx1, sy1 & sxy1)
    sx1  = sqrt(sum((1-d[['pt']])*(d[["beta.exposure"]]^2))/(length(d[["beta.exposure"]])-sum(d[['pt']])))
    sy1  = sqrt(sum((1-d[['pt']])*(d[["beta.outcome"]]^2))/(length(d[["beta.outcome"]])-sum(d[['pt']])))
    sxy1 = sum((1-d[['pt']])*d[["beta.exposure"]]*d[["beta.outcome"]])/(length(d[["beta.exposure"]])-sum(d[['pt']]))
    if (abs(sxy1) > 0.75*sx1*sy1)  sxy1 = sign(sxy1)*0.75*sx1*sy1

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
  entropy = mean(d[clusters=="Hunted", pt], na.rm=TRUE)

  # store the results
  out_cols <- c("SNP","BETAi","BETAp","CLUSTER")
  data.table::setnames(d, c("SNP","beta.exposure","beta.outcome","clusters"), out_cols)
  d[, names(d)[!names(d)%in%out_cols] := NULL]
  result <- list(fit     = d,
                 b       = b,
                 bse     = bse,
                 b_CI    = c(b - 1.96*bse, b + 1.96*bse),
                 entropy = entropy,
                 iters   = iter_data)

  return(result)
}



#' @title Create an animation of Slope-hunter iterations
#' @param x a ColliderBias object
#' @return a plot
#' @export
#' @importFrom gganimate transition_states
#'
plot_slopehunter_iters <- function(x) UseMethod("plot_slopehunter_iters")
#' @rdname plot_slopehunter_iters
#' @export
plot_slopehunter_iters.ColliderBias <- function(x) {

  pval.exposure = beta.exposure = beta.outcome = pt = clusters = iter = NULL

  # filter by the incidence threshold
  d <- x$harmonised[pval.exposure < x$ip, ]

  # run the EM algorithm
  result <- shclust(d, x$pi0, x$sigma0, collect_iters=TRUE)

  # make the plot
  result <- data.table::rbindlist(result$iters, idcol="iter")

  # Make a ggplot, but add frame=year: one image per year
  p <- ggplot2::ggplot() +
    geom_point(data=result, mapping=ggplot2::aes(x=beta.exposure, y=beta.outcome, fill=pt, color=clusters, size=clusters), shape=21) +
    viridis::scale_fill_viridis(limits=c(0,1)) +
    ggplot2::scale_size_manual(values = c("Hunted"=1.5, "Pleiotropic"=0.5)) +
    ggplot2::scale_alpha_manual(values = c("Hunted"=1.0, "Pleiotropic"=0.1)) +
    ggplot2::scale_shape_manual(values = c("Hunted"=21, "Pleiotropic"=1)) +
    ggplot2::scale_color_manual(values = c("Hunted"="red", "Pleiotropic"="#44015410")) +
    ggplot2::theme_bw() +
    # gganimate specific bits:
    ggplot2::labs(title = 'Slope-hunter iteration: {closest_state}',
                  x = 'BETA incidence',
                  y = 'BETA progression',
                  color = "Cluster",
                  fill = "Likelihood Gi SNP") +
    ggplot2::guides(size = "none") +
    gganimate::transition_states(iter, transition_length=0.1, state_length=1)

  return(p)
}



valid_params <- function(methods, ip, pi0, sigma0, bootstraps) {

  # check methods
  methods <- match.arg(methods, c("slopehunter","dudbridge","ivw_mr"), several.ok=TRUE)

  # check ip
  if(any(c("slopehunter","ivw_mr") %in% methods)) {
    stopifnot("Slopehunter amd IVW-MR require a numeric `ip` p-value, between 0 and 1, by which to filter the incidence variants" = is.numeric(ip) & ip<=1 & ip>=0)
  }

  # check pi0, sigma0, bootstraps
  if("slopehunter" %in% methods) {
    stopifnot("Slopehunter requires a numeric `pi0` value, between 0 and 1, for the weight of the mixture component that represents the cluster of SNPs affecting incidence only" = is.numeric(pi0) & pi0<=1 & pi0>=0)
    stopifnot("Slopehunter requires an initialising numeric `sigma0` value for the for the covariance between incidence and progression BETAs" = is.numeric(sigma0))
    stopifnot("Slopehunter requires a numeric whole-number `bootstraps` value for the number of bootstrap samples drawn to estimate the standard error of the adjustment factor" = is.numeric(bootstraps) & bootstraps%%1==0)
  }

  # if get here, checks ok
  return(TRUE)
}




#' @title Dudbridge CWLS method
#' @description
#' The Dudbridge corrected weighted least-squares method for assessing collider bias.
#'
#' @param x a ColliderBias object
#' @return a ColliderBias object with updated parameter and results fields
#' @export
dudbridge <- function(x) UseMethod("dudbridge")
#' @rdname dudbridge
#' @export
dudbridge.ColliderBias <- function(x) {
#
#   # copy the harmonised data in
#   d <- x$harmonised
#
#
#
#
#
#   cwls_correction <- MendelianRandomization::mr_ivw(MendelianRandomization::mr_input(
#     bx = pruned_harmonised_effects$BETA.incidence,
#     bxse = pruned_harmonised_effects$SE.incidence,
#     by = pruned_harmonised_effects$BETA.prognosis,
#     byse = pruned_harmonised_effects$SE.prognosis
#   ))
#
#   pruned_harmonised_effects$weights <- 1 / pruned_harmonised_effects$SE.prognosis^2
#
#   weighting <- (sum(pruned_harmonised_effects$weights * pruned_harmonised_effects$BETA.incidence^2)) /
#     (
#       (sum(pruned_harmonised_effects$weights * pruned_harmonised_effects$BETA.incidence^2)) -
#         (sum(pruned_harmonised_effects$weights * pruned_harmonised_effects$SE.incidence^2))
#     )
#
#   cwls_estimated_slope <- cwls_correction$Estimate * weighting
#   cwls_estimated_standard_error <- cwls_correction$StdError * weighting
#
#   collider_bias_results <- rbind(collider_bias_results,
#                                  data.table::data.table(
#                                    METHOD = "cwls",
#                                    P_VALUE_THRESHOLD = NA,
#                                    SNPS_USED = nrow(clumped_snps),
#                                    BETA = cwls_estimated_slope,
#                                    SE = cwls_estimated_standard_error,
#                                    ENTROPY = NA,
#                                    PLEIOTROPIC = NA
#                                  )
#   )
#


}



#' @title IVW MR
#' @description
#' The Inverse variance weighted Mendelian Randomisation method for assessing collider bias.
#' @param x a ColliderBias object
#' @inheritParams ColliderBias
#' @return a ColliderBias object with updated parameter and results fields
#' @export
ivw_mr <- function(x,ip) UseMethod("ivw_mr")
#' @rdname ivw_mr
#' @export
ivw_mr.ColliderBias <- function(x, ip = 0.001) {
#
#   # filter by the incidence threshold
#   d <- x$harmonised["pval.exposure" <= ip, ]
#
#   # exit if there are no significant incidence variants at this threshold
#   if(nrow(d)==0) {
#     warning(paste0("No significant incidence SNPs at `ip` ", ip, "\n"))
#     return(x)
#   }
#

    #
    # # determine which to use for MR
    # x <- subset(x1, mr_keep)
    #
    # # apply the method
    # res <- lapply(method_list, function(meth)
    # {
    #   get(meth)(x$beta.exposure, x$beta.outcome, x$se.exposure, x$se.outcome, parameters)
    # })

    # the method
    # Inverse variance weighted regression
    #
    # The default multiplicative random effects IVW estimate.
    # The standard error is corrected for under dispersion
    # Use the  mr_ivw_mre  function for estimates that don't correct for under dispersion.
    # #    }
    # mr_ivw <- function(b_exp, b_out, se_exp, se_out, parameters=default_parameters())
    # {
      # if(sum(!is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) & !is.na(se_out)) < 2)
      #   return(list(b=NA, se=NA, pval=NA, nsnp=NA))
      #
      # ivw.res <- summary(stats::lm(b_out ~ -1 + b_exp, weights = 1/se_out^2))
      # b <- ivw.res$coef
      # se <- ivw.res$coef["b_exp","Std. Error"]/min(1,ivw.res$sigma) #sigma is the residual standard error
      # pval <- 2 * stats::pnorm(abs(b/se), lower.tail=FALSE)
      # Q_df <- length(b_exp) - 1
      # Q <- ivw.res$sigma^2 * Q_df
      # Q_pval <- stats::pchisq(Q, Q_df, lower.tail=FALSE)
      # # from formula phi =  Q/DF rearranged to to Q = phi*DF, where phi is sigma^2
      # # Q.ivw<-sum((1/(se_out/b_exp)^2)*(b_out/b_exp-ivw.reg.beta)^2)
      # return(list(b = b, se = se, pval = pval, nsnp=length(b_exp), Q = Q, Q_df = Q_df, Q_pval = Q_pval))
    # }




    # default_parameters <- function()
    # {
    #   list(
    #     test_dist = "z",
    #     nboot = 1000,
    #     Cov = 0,
    #     penk = 20,
    #     phi = 1,
    #     alpha = 0.05,
    #     Qthresh = 0.05,
    #     over.dispersion = TRUE,
    #     loss.function = "huber",
    #     shrinkage = FALSE
    #   )
    # }

  # }


}




