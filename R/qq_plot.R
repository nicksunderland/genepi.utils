#' @title QQ plot
#' @param gwas a data.frame like object or valid file path
#' @param pval_col the P value column
#' @param title a string, the title for the plot
#' @param subtitle a string, the subtitle for the plot
#' @param facet_grid_row_col a string, the column name in `gwas` by which to facet the plot (rows)
#' @param facet_grid_col_col a string, the column name in `gwas` by which to facet the plot (cols)
#' @param colours a 2 element list of colour codes (1-the uncorrected points, 2-the GC corrected points)
#' @param plot_corrected a logical, whether to apply and plot the lambda correction
#' @param facet_nrow an integer, passed to facet_wrap, the number of rows to facet by (if only facet_grid_row_col is provided)
#' @param facet_ncol an integer, passed to facet_wrap, the number of cols to facet by (if only facet_grid_col_col is provided)
#'
#' @return a ggplot
#' @export
#' @importFrom stats qchisq median as.formula
#'
qq_plot <- function(gwas,
                    pval_col = "p",
                    colours = list(raw="#2166AC"), #, corrected="#B2182B"),
                    title = NULL,
                    subtitle = NULL,
                    plot_corrected = FALSE,
                    facet_grid_row_col = NULL,
                    facet_grid_col_col = NULL,
                    facet_nrow = NULL,
                    facet_ncol = NULL) {

  # conversion from GWAS object
  if (inherits(gwas, "genepi.utils::GWAS")) {
    gwas <- genepi.utils::as.data.table(gwas)
  }

  # silence R CMD check warnings
  P <- adj_P <- adj_chisq <- adj_observed <- chisq <- clower <- corrected <- cupper <- expected <- label <- lambda <- observed <- NULL

  # get the gwas
  gwas <- import_table(gwas)

  # checks
  stopifnot("Colours must be length 1, or length 2 if `plot_corrected` is TRUE" = length(colours)==1 | length(colours)==2 & plot_corrected)
  stopifnot("There must be a `P` column" = pval_col %in% names(gwas))
  data.table::setnames(gwas, pval_col, "P")
  if(any(gwas[, is.na(P) | is.nan(P) | is.null(P) | is.infinite(P) | P<0.0 | P>1.0])) {
    invalid_idx <- which(gwas[, is.na(P) | is.nan(P) | is.null(P) | is.infinite(P) | P<0.0 | P>1.0])
    print(gwas[invalid_idx, ])
    warning("Invalid P values, please remove")
  }
  if(!is.null(facet_grid_row_col)) {
    stopifnot("The row faceting column, if given, must be in the data source" = facet_grid_row_col %in% names(gwas))
  }
  if(!is.null(facet_grid_col_col)) {
    stopifnot("The col faceting column, if given, must be in the data source" = facet_grid_col_col %in% names(gwas))
  }

  # gather the data for the QQ plot
  gwas[, chisq := stats::qchisq(1-P, 1)]

  facets <- c(facet_grid_row_col, facet_grid_col_col)

  # calculate lambda
  gwas[, lambda := stats::median(chisq) / qchisq(0.5, 1), by=facets]

  # new / adjusted chisq statistic and P
  gwas[, adj_chisq := chisq/lambda]
  gwas[, adj_P     := stats::pchisq(adj_chisq, 1, lower.tail=FALSE)]

  # calculate observed and lambda
  plot_data <- gwas[, list(lambda   = median(chisq) / qchisq(0.5, 1),
                           observed = -log10(sort(P, decreasing=FALSE)),
                           adj_observed = -log10(sort(adj_P, decreasing=FALSE)),
                           expected = -log10(stats::ppoints(.N)),
                           clower   = -log10(stats::qbeta(p = (1 - 0.95) / 2, shape1 = 1:.N, shape2 = .N:1)),
                           cupper   = -log10(stats::qbeta(p = (1 + 0.95) / 2, shape1 = 1:.N, shape2 = .N:1))),
                    by=facets]

  # generate QQ-plot

  # axis labels
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))

  # lambda labels
  lambda_labels <- plot_data[, list(lambda = lambda[1]), by=facets]
  lambda_labels[, expected := 3.0]
  lambda_labels[, observed := 0.25]
  lambda_labels[, label    := paste0("\u03BB = ", format(round(lambda, 3), nsmall=3))]

  # downsample for faster plotting
  data.table::setorder(plot_data, expected)
  plot_data <- plot_data[sample(x       = 1:.N,
                                replace = TRUE, # replace samples, for speed
                                size    = ifelse(5000>.N,.N,5000)+ifelse(.N-5000<0,0,floor((.N-5000)*0.05)), # size to return is 5000 + 5% of original data
                                prob    = ifelse(rep(5000>.N,.N),rep(1,.N),c(seq(0,1,length.out=.N-floor(.N*0.05)), rep(1,floor(.N*0.05))))), ] # probability of choosing row is lowest for 0-95th percentile, and certain for the top 5th centile

  # plot
  plot <- plot_data |>
    ggplot2::ggplot(ggplot2::aes(x = expected, y = observed)) +
    ggplot2::geom_point(size = 0.5, color=colours[[1]]) +
    ggplot2::geom_ribbon(mapping = ggplot2::aes(x = expected, ymin = clower, ymax = cupper), alpha = 0.1, color="transparent") +
    ggplot2::geom_abline(slope=1, intercept=0, color="darkgrey", linetype = "dotted") +
    ggplot2::geom_text(data = lambda_labels, ggplot2::aes(x=expected, y=observed, label=label), color="black", show.legend = FALSE) +
    ggplot2::labs(title = title,
                  subtitle = subtitle,
                  x = log10Pe,
                  y = log10Po,
                  color = NULL) +
    ggplot2::theme_light() +
    ggplot2::theme(aspect.ratio = 1)

  if(plot_corrected) {
    plot_data[, corrected := "Corrected"]
    plot <- plot +
      ggplot2::geom_point(data = plot_data,
                          mapping = ggplot2::aes(x = expected, y = adj_observed, color=corrected), size = 0.5, inherit.aes=FALSE) +
      ggplot2::scale_color_manual(values=c(colours[[2]])) +
      ggplot2::theme(legend.position="top") +
      ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=3)))
  }

  if(!is.null(facet_grid_row_col) & !is.null(facet_grid_col_col)) {
    plot <- plot + ggplot2::facet_grid(facet_grid_row_col ~ facet_grid_col_col)
  } else if(!is.null(facet_grid_row_col)) {
    plot <- plot + ggplot2::facet_wrap(stats::as.formula(paste("~", facet_grid_row_col)), nrow=facet_nrow, ncol=facet_ncol)
  } else if(!is.null(facet_grid_col_col)) {
    plot <- plot + ggplot2::facet_wrap(stats::as.formula(paste("~", facet_grid_col_col)), nrow=facet_nrow, ncol=facet_ncol)
  }

  # return
  return(plot)

}
