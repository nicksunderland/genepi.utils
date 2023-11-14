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

  message("Creating ColliderBias object...")

  # checks
  stopifnot("Invalid `ColliderBias` parameters" = valid_params(methods, ip, pi0, sxy1, bootstraps))

  # format
  formatted_gwases <- format_gwas_for_collider_analysis(gwas_i, gwas_p)

  # harmonise
  harmonised <- TwoSampleMR::harmonise_data(formatted_gwases[["gwas_i"]], formatted_gwases[["gwas_p"]]) |> data.table::as.data.table()

  # package the data into a structure of class ColliderBias
  s3_struct <- structure(
    .Data = list(
      gwas_i     = formatted_gwases[["gwas_i"]],
      gwas_p     = formatted_gwases[["gwas_p"]],
      harmonised = harmonised,
      methods    = c("slopehunter","dudbridge","ivw_mr"),
      ip         = ip,
      pi0        = pi0,
      sxy1       = sxy1,
      bootstraps = bootstraps,
      results    = data.table::data.table()
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
analyse <- function(x,methods,ip,pi0,sxy1,bootstraps) UseMethod("analyse")
#' @rdname analyse
#' @export
analyse.ColliderBias <- function(x,methods=NULL,ip=NULL,pi0=NULL,sxy1=NULL,bootstraps=NULL) {

  if(is.null(methods)) methods <- x$methods
  if(is.null(ip)) ip <- x$ip
  if(is.null(pi0)) pi0 <- x$pi0
  if(is.null(sxy1)) sxy1 <- x$sxy1
  if(is.null(bootstraps)) bootstraps <- x$bootstraps

  params <- lapply(methods, function(m) {
    if(m=="slopehunter") return(expand.grid(method="slopehunter", ip=ip, pi0=pi0, sxy1=sxy1, bootstraps=bootstraps))
    if(m=="dudbridge") return(expand.grid(method="dudbridge", ip=NA_real_, pi0=NA_real_, sxy1=NA_real_, bootstraps=NA_real_))
    if(m=="ivw_mr") return(expand.grid(method="ivw_mr", ip=ip, pi0=NA_real_, sxy1=NA_real_, bootstraps=NA_real_))
  }) |> data.table::rbindlist()

  # results from each combination of parameters
  results <- data.table::data.table()
  for(i in 1:nrow(params)) {

    res <- do.call(what = as.character(params$method[[i]]),
                   args = list(x    = x,
                               ip   = params$ip[[i]],
                               pi0  = params$pi0[[i]],
                               sxy1 = params$sxy1[[i]],
                               bootstraps = params$bootstraps[[i]]))

    results <- rbind(results, as.data.frame(res))
  }

  # add the results dt to the object
  x[["results"]] <- results

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
#' @return a list with the method parameters and results
#' @export
slopehunter <- function(x,ip,pi0,sxy1,bootstraps,seed) UseMethod("slopehunter")
#' @rdname slopehunter
#' @export
slopehunter.ColliderBias <- function(x          = NULL,
                                     ip         = 0.001,
                                     pi0        = 0.6,
                                     sxy1       = 1e-5,
                                     bootstraps = 100,
                                     seed       = 2023) {


  # inc <- SlopeHunter::format_data(x$gwas_i,
  #                          type = "incidence",
  #                          snp_col = "SNP",
  #                          beta_col = "beta.exposure",
  #                          se_col = "se.exposure",
  #                          pval_col = "pval.exposure",
  #                          eaf_col = "eaf.exposure",
  #                          effect_allele_col = "effect_allele.exposure",
  #                          other_allele_col = "other_allele.exposure",
  #                          chr_col = "chr.exposure",
  #                          pos_col = "pos.exposure")
  #
  # pro <- SlopeHunter::format_data(x$gwas_p,
  #                               type = "prognosis",
  #                               snp_col = "SNP",
  #                               beta_col = "beta.outcome",
  #                               se_col = "se.outcome",
  #                               pval_col = "pval.outcome",
  #                               eaf_col = "eaf.outcome",
  #                               effect_allele_col = "effect_allele.outcome",
  #                               other_allele_col = "other_allele.outcome",
  #                               chr_col = "chr.outcome",
  #                               pos_col = "pos.outcome")
  #
  # h <- SlopeHunter::harmonise_effects(inc,pro) |> as.data.frame()
  #
  # sh <- SlopeHunter:: hunt(h, xp_thresh=ip, init_pi=pi0, init_sigmaIP=sxy1)
  #
  # res <- analysis_result(method ="slopehunter",
  #                        ip     = ip,
  #                        pi0    = pi0,
  #                        sxy1   = sxy1,
  #                        b      = sh$b,
  #                        bse    = sh$bse,
  #                        entropy= sh$entropy)







  pval.exposure = NULL

  # filter by the incidence threshold
  d <- x$harmonised[pval.exposure < ip, ]

  # exit if there are no significant incidence variants at this threshold
  if(nrow(d)==0) {
    warning(paste0("No significant incidence SNPs at `ip` ", ip, "\n"))
    return(analysis_result(method ="slopehunter",
                           ip     = ip,
                           pi0    = pi0,
                           sxy1   = sxy1))
  }

  # run the EM algorithm
  result <- shclust(d, pi0=pi0, sxy1=sxy1, collect_iters=FALSE)

  # run lots more times bootstrapping to estimate SE
  if (bootstraps > 0){
    b.boots = vector("numeric", bootstraps)
    for(i in 1:bootstraps){
      if(i%%10==0) cat("Bootstrapping... ", i, "of", bootstraps, "\n")
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
    result$bse_LB = result$b - 1.96*result$bse
    result$bse_UB = result$b + 1.96*result$bse
  }

  # overall result to return
  res <- analysis_result(method ="slopehunter",
                         ip     = ip,
                         pi0    = pi0,
                         sxy1   = sxy1,
                         b      = result$b,
                         bse    = result$bse,
                         entropy= result$entropy)
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

  f0 = f1 = pt = clusters = NULL

  d <- data.table::copy(d)

  # set the initial sdev and cov
  sx0 = sx1 = stats::sd( d[["beta.exposure"]] )
  sy0 = sy1 = stats::sd( d[["beta.outcome"]] )
  dir0 = sign(  stats::cov(d[["beta.exposure"]], d[["beta.outcome"]]) )
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
      iter_bs <- c(iter_bs, list(dir0*sy0/sx0))
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
                 bse_LB  = b - 1.96*bse,
                 bse_UB  = b + 1.96*bse,
                 entropy = entropy,
                 iters   = iter_data,
                 iters_b = iter_bs)

  return(result)
}


#' @title Dudbridge CWLS method
#' @description
#' The Dudbridge corrected weighted least-squares method for assessing collider bias.
#' @param ... sink for unneeded params
#' @param x a ColliderBias object
#' @return a ColliderBias object with updated parameter and results fields
#' @export
#' @importFrom MendelianRandomization mr_ivw
dudbridge <- function(x, ...) UseMethod("dudbridge")
#' @rdname dudbridge
#' @export
dudbridge.ColliderBias <- function(x, ...) {

  # dots just so we can call all the methods with the same parameters
  ignore <- list(...)

  # harmonise the filtered incidence variants against the progression
  d <- TwoSampleMR::harmonise_data(x$gwas_i, x$gwas_p) |> data.table::as.data.table()

  # calculate the correction
  cwls_correction <- MendelianRandomization::mr_ivw(
    MendelianRandomization::mr_input(
      bx   = d[["beta.exposure"]],
      bxse = d[["se.exposure"]],
      by   = d[["beta.outcome"]],
      byse = d[["se.outcome"]]
    )
  )

  weights <- 1 / d[["se.outcome"]]^2

  weighting <-  (sum(weights*d[["beta.exposure"]]^2)) /
               ((sum(weights*d[["beta.exposure"]]^2)) - (sum(weights*d[["se.exposure"]]^2)))

  cwls_estimated_slope <- cwls_correction$Estimate * weighting
  cwls_estimated_standard_error <- cwls_correction$StdError * weighting

  # overall result to return
  res <- analysis_result(method ="dudbridge",
                         b      = cwls_estimated_slope,
                         bse    = cwls_estimated_standard_error)

  return(res)

}



#' @title IVW MR
#' @description
#' The Inverse variance weighted Mendelian Randomisation method for assessing collider bias.
#' @param x a ColliderBias object
#' @param ... sink for unneeded params
#' @inheritParams ColliderBias
#' @return a ColliderBias object with updated parameter and results fields
#' @export
#' @importFrom TwoSampleMR harmonise_data mr
ivw_mr <- function(x,ip, ...) UseMethod("ivw_mr")
#' @rdname ivw_mr
#' @export
ivw_mr.ColliderBias <- function(x,ip = 0.001, ...) {

  ignore <- list(...)

  pval.exposure = NULL

  # filter by the incidence threshold
  d <- x$harmonised[pval.exposure < ip, ]

  # exit if there are no significant incidence variants at this threshold
  if(nrow(d)==0) {
    warning(paste0("No significant incidence SNPs at `ip` ", ip, "\n"))
    return(analysis_result(method ="ivw_mr", ip = ip))
  }

  # run the MR
  mr_results <- TwoSampleMR::mr(d, method_list = c("mr_ivw"))

  # overall result to return
  res <- analysis_result(method ="ivw_mr",
                         b   = mr_results$b,
                         bse = mr_results$se,
                         ip  = ip)

  return(res)
}


#' @title Apply collider bias correction factor
#' @description
#' Apply collider bias adjustment.
#' @param x a ColliderBias object
#' @param b the calculated adjustment factor `b`
#' @param se the calculated b standard error
#' @return a data.table
#' @export
#' @importFrom stats pchisq
#'
apply_correction <- function(x,b,se) UseMethod("apply_correction")
#' @rdname apply_correction
#' @export
apply_correction.ColliderBias <- function(x,b,se) {

  adjusted_beta = adjusted_se = adjusted_p = NULL

  x$harmonised[, adjusted_beta := beta.outcome - (b * beta.exposure)]
  x$harmonised[, adjusted_se   := sqrt(
    (se.outcome^2) +
    ((b^2) * (se.exposure^2)) +
    ((se.exposure^2) * (se^2)) +
    ((se.exposure^2) * (se^2))
    )
  ]

  x$harmonised[, adjusted_p := stats::pchisq( (adjusted_beta / adjusted_se)^2, 1, lower.tail=FALSE)]

  data.table::setkey(x$harmonised, "SNP")
  data.table::setkey(x$gwas_p, "SNP")

  gwas <- x$gwas_p[x$harmonised, c("adjusted_beta", "adjusted_se", "adjusted_p") := .(i.adjusted_beta,
                                                                                      i.adjusted_se,
                                                                                      i.adjusted_p), on="SNP"]

  return(gwas)
}

#' @title Create an animation of Slope-hunter iterations
#' @inheritParams slopehunter
#' @param x_lims x axis limits
#' @param y_lims y axis limits
#' @return a plot
#' @export
#' @importFrom gganimate transition_states
#'
plot_slopehunter_iters <- function(x,ip,pi0,sxy1,bootstraps,seed,x_lims,y_lims) UseMethod("plot_slopehunter_iters")
#' @rdname plot_slopehunter_iters
#' @export
plot_slopehunter_iters.ColliderBias <- function(x,
                                                ip         = 0.001,
                                                pi0        = 0.6,
                                                sxy1       = 1e-5,
                                                bootstraps = 0,
                                                seed       = 2023,
                                                x_lims     = NULL,
                                                y_lims     = NULL) {

  pval.exposure = beta.exposure = beta.outcome = pt = clusters = iter = NULL

  # filter by the incidence threshold
  d <- x$harmonised[pval.exposure < x$ip, ]

  # run the EM algorithm
  result <- shclust(d, x$pi0, x$sxy1, collect_iters=TRUE)

  # make the plot
  result_dat <- data.table::rbindlist(result$iters, idcol="iter")

  # slope data
  result_b <- data.table::data.table(slope = unlist(result$iters_b),
                                     intercept = 0,
                                     iter=1:length(result$iters_b),
                                     label = paste0("b = ", as.character(round(unlist(result$iters_b), digits=2))),
                                     xpos = Inf, ypos = Inf, hjustvar = 1.2, vjustvar = 1.4)

  # Make a ggplot, but add frame=year: one image per year
  p <- ggplot2::ggplot() +
    ggplot2::geom_text(data=result_b, mapping=ggplot2::aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=label)) +
    ggplot2::geom_abline(data=result_b, mapping=ggplot2::aes(slope=slope,intercept=intercept), color="darkred") +
    ggplot2::geom_point(data=result_dat, mapping=ggplot2::aes(x=beta.exposure, y=beta.outcome, fill=pt, color=clusters, size=clusters), shape=21) +
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
#' @param x a ColliderBias object with populated results field
#' @return a plot
#' @export
#' @import ggplot2
#'
plot_p_threshold <- function(x) UseMethod("plot_p_threshold")
#' @rdname plot_p_threshold
#' @export
plot_p_threshold.ColliderBias <- function(x) {

  # b by p-values
  p <- ggplot2::ggplot(data=x$results,
                  mapping=ggplot2::aes(x=ip, y=b)) +
    ggplot2::geom_ribbon(mapping=ggplot2::aes(ymin=b-(1.96*bse), ymax=b+(1.96*bse), fill=method), alpha=0.3) +
    ggplot2::geom_point(mapping=ggplot2::aes(color=method)) +
    ggplot2::geom_line(mapping=ggplot2::aes(color=method)) +
    # viridis::scale_fill_viridis(discrete = TRUE) +
    # viridis::scale_color_viridis(discrete = TRUE) +
    ggplot2::scale_x_continuous(trans='log10') +
    ggplot2::labs(title = 'Bias correction factor by p-value threshold',
                  x = 'P-value `ip` threshold',
                  y = 'b value (correction slope)',
                  color = "Method",
                  fill = "95%CI")

  return(p)
}


#' @title Plot the Slope-hunter scatter plot
#' @param x a ColliderBias object with populated results field
#' @return a plot
#' @export
#' @import ggplot2
#'
plot_slopehunter <- function(x) UseMethod("plot_slopehunter")
#' @rdname plot_slopehunter
#' @export
plot_slopehunter.ColliderBias <- function(x) {

  BETAi = BETAp = CLUSTER = hjustvar = vjustvar = xpos = ypos = label = slope = intercept = NULL

  # results for each pvalue
  dat <- lapply(x$ip, function(thresh) {
    d <- x$harmonised[pval.exposure < thresh, ]
    result <- shclust(d, pi0=x$pi0, sxy1=x$sxy1, collect_iters=FALSE)
    return(result)
  })

  # slope data
  b_dat <- lapply(dat, function(x) data.table::data.table(slope=x$b,
                                                          intercept=0,
                                                          label = paste0("b = ", as.character(round(x$b, digits=2))),
                                                          xpos = Inf,
                                                          ypos = Inf,
                                                          hjustvar = 1.2,
                                                          vjustvar = 1.4)) |> setNames(x$ip) |> data.table::rbindlist(idcol="ip")
  b_dat[, ip := factor(ip, levels=x$ip)]

  # point data
  dat <- lapply(dat, function(x) x$fit) |> setNames(x$ip) |> data.table::rbindlist(idcol="ip")
  dat[, ip := factor(ip, levels=x$ip)]

  # facet plot of different P value thresholds
  g <- ggplot2::ggplot(data=dat, ggplot2::aes(BETAi, BETAp, color=CLUSTER)) +
    ggplot2::geom_point() +
    ggplot2::geom_abline(data=b_dat, ggplot2::aes(slope=slope, intercept=intercept), color = 'black') +
    ggplot2::geom_text(data=b_dat, ggplot2::aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=label), color="black") +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "SlopeHunter analysis - by p-value threshold (\u03BB)",
                  x = "beta (indicence)",
                  y = "beta (progression)") +
    ggplot2::facet_wrap(~ ip)

  return(g)
}



format_gwas_for_collider_analysis <- function(gwas_i, gwas_p) {

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
  message("Formatting GWAS inputs")
  gwas_i <- do.call(TwoSampleMR::format_data, c(list(dat=gwas_i, type="exposure"), i_cols)) |> data.table::as.data.table()
  gwas_p <- do.call(TwoSampleMR::format_data, c(list(dat=gwas_p, type="outcome"), p_cols)) |> data.table::as.data.table()

  # return
  return(list("gwas_i"=gwas_i, "gwas_p"=gwas_p))

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


analysis_result <- function(method=NA_character_,
                            ip     = NA_real_,
                            pi0    = NA_real_,
                            sxy1   = NA_real_,
                            b      = NA_real_,
                            bse    = NA_real_,
                            entropy= NA_real_) {
  s3_struct <- structure(
    .Data = list(
      method = method,
      ip     = ip,
      pi0    = pi0,
      sxy1   = sxy1,
      b      = b,
      bse    = bse,
      entropy= entropy
    ),
    class = "list"
  )
}








