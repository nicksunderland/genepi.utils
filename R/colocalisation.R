#' @title Coloc probability plot
#' @description
#' A plotting wrapper for the `coloc` package. Produces a ggplot for either
#' the prior or posterior probability sensitivity analyses. See the
#' [coloc](https://chr1swallace.github.io/coloc/articles/a04_sensitivity.html)
#' package vignettes for details.
#' @param coloc coloc object, output from `coloc::coloc.abf()`
#' @param rule a string, a valid rule indicating success e.g. "H4 > 0.5"
#' @param type a string, either `prior` or `posterior`
#' @param row an integer, row in a `coloc.susie` or `coloc.signals` object
#' @return a ggplot
#' @export
#' @references [coloc](https://chr1swallace.github.io/coloc/articles/a04_sensitivity.html)
plot_coloc_probabilities <- function(coloc, rule="H4 > 0.5", type="prior", row=1) {

  # RCMD check warnings
  h <- x <- NULL

  # checks
  type <- match.arg(type, choices=c('prior','posterior'))
  stopifnot("list" %in% class(coloc))
  stopifnot("priors" %in% names(coloc))
  stopifnot("summary" %in% names(coloc))
  rule.init <- rule
  rule <- gsub("(H.)","PP.\\1.abf",rule,perl=TRUE)

  ## extract the results
  # multiple signals?
  if(data.table::is.data.table(coloc$summary)) {
    if(!(row %in% 1:nrow(coloc$summary)))
      stop("row must be between 1 and ",nrow(coloc$summary))
    pp <- unlist(c(coloc$summary[row,grep("PP|nsnp",names(coloc$summary)),with=FALSE]))
    if(paste0("SNP.PP.H4.row",row) %in% names(coloc$results)) {
      coloc$results[["SNP.PP.H4"]] <- coloc$results[[paste0("SNP.PP.H4.row",row)]]
    }
    pp <- unlist(c(coloc$summary[row,grep("PP|nsnp",names(coloc$summary)),with=FALSE]))
  } else {
    pp <- coloc$summary
  }
  results <- coloc$results
  p12     <- coloc$priors["p12"]
  p1      <- coloc$priors["p1"]
  p2      <- coloc$priors["p2"]

  # do the tests
  check   <- function(pp) { with(as.list(pp),eval(parse(text=rule))) }
  testp12 <- 10^seq(log10(p1*p2),log10(min(p1,p1)),length.out=100)
  testH   <- prior.snp2hyp(pp["nsnps"],p12=testp12,p1=p1,p2=p2)
  testpp  <- data.table::as.data.table(prior.adjust(summ=pp,newp12=testp12,p1=p1,p2=p2,p12=p12))
  names(testpp) <- gsub("(H.)","PP.\\1.abf",colnames(testpp),perl=TRUE)
  pass    <- check(testpp)
  w       <- which(pass)

  # base plot to add to
  p <- ggplot2::ggplot()

  # which plot to create
  if(type=="prior") {

    # create the data.frame for prior probabilities
    prior_dat <- data.table::data.table(testH)
    prior_dat[, x := testp12]
    prior_dat <- data.table::melt(prior_dat,
                                  id.vars = c("x"),
                                  measure.vars = paste0("H",0:4),
                                  variable.name = "h",
                                  value.name = "p")
    prior_dat[, h := factor(h, levels=paste0("H",0:4))]

    p_max <- max(prior_dat$p[prior_dat$h != "H0"], na.rm=T)

    if(any(pass)) {
      p <- p + ggplot2::geom_rect(ggplot2::aes(xmin=testp12[min(w)], xmax=testp12[max(w)], ymin=0, ymax=p_max), fill="green", alpha=0.1, color="green")
    }

    p <- p +
      ggplot2::geom_point(data = prior_dat, mapping = ggplot2::aes(x=x, y=p, fill=h), color="#333333", shape=21, size=3, stroke=0.25) +
      ggplot2::lims(y=c(0,p_max)) +
      ggplot2::labs(title = "Prior probabilities", subtitle = paste("shaded region:",rule.init))

  } else if(type=="posterior") {

    # create the data.frame for posterior probabilities
    posterior_dat <- data.table::data.table(as.matrix(testpp))
    posterior_dat[, x := testp12]
    posterior_dat <- data.table::melt(posterior_dat,
                                      id.vars = c("x"),
                                      measure.vars = paste0("PP.H",0:4,".abf"),
                                      variable.name = "h",
                                      value.name = "p")
    posterior_dat[, h := factor(h, levels=paste0("PP.H",0:4,".abf"))]

    p_max <- max(posterior_dat$p, na.rm=T)

    if(any(pass)) {
      p <- p + ggplot2::geom_rect(ggplot2::aes(xmin=testp12[min(w)], xmax=testp12[max(w)], ymin=0, ymax=p_max), fill="green", alpha=0.1, color="green")
    }

    p <- p +
      ggplot2::geom_point(data = posterior_dat, mapping = ggplot2::aes(x=x, y=p, fill=h), color="#333333", shape=21, size=3, stroke=0.25) +
      ggplot2::lims(y=c(0,p_max)) +
      ggplot2::labs(title = "Posterior probabilities", subtitle = paste("shaded region:",rule.init))

  } else {

    stop("type parameter error")

  }

  # add the common elements
  p <- p +
    ggplot2::scale_fill_manual(values = c("#ffffffff",viridis::viridis(5,alpha=1)[-1])) +
    ggplot2::geom_vline(xintercept = p12, linetype="dashed", color="gray") +
    ggplot2::scale_x_continuous(trans='log10') +
    ggplot2::labs(x="p12", y="Prob") +
    ggplot2::annotate("text", x=p12, y=0.5, label="results", angle=90, color="gray40") +
    ggplot2::theme(legend.position = "top",
                   legend.title    = ggplot2::element_blank())

  # return
  return(p)
}


# copied from coloc package to make the above work
prior.adjust <- function(summ,newp12,p1=1e-4,p2=1e-4,p12=1e-6) {
  if(is.list(summ) && "summary" %in% names(summ))
    summ <- summ$summary
  if(!identical(names(summ), c("nsnps", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")))
    stop("not a coloc summary vector")
  ## back calculate likelihoods
  f <- function(p12)
    prior.snp2hyp(summ["nsnps"],p12=p12,p1=p1,p2=p2)
  pr1 <- f(newp12)
  pr0 <- matrix(f(p12),nrow=nrow(pr1),ncol=ncol(pr1),byrow=TRUE)
  newpp <- matrix(summ[-1],nrow=nrow(pr1),ncol=ncol(pr1),byrow=TRUE) * pr1/pr0 # prop to, not equal to
  newpp/rowSums(newpp)
}


# copied from coloc package to make the above work
prior.snp2hyp <- function(nsnp,p12=1e-6,p1=1e-4,p2=1e-4) {
  if(any(p12<p1*p2) || any(p12 > p1) || any(p12 > p2)) {
    stop("Sensitivity plot input probability issue - ensure NOT `any(p12<p1*p2) || any(p12 > p1) || any(p12 > p2)`")
  }
  tmp <- cbind(nsnp * p1,
               nsnp * p2,
               nsnp * (nsnp-1) * p1 * p2,
               nsnp * p12)
  tmp <- cbind(1-rowSums(tmp),tmp)
  colnames(tmp) <- paste0("H",0:4)
  tmp
}
