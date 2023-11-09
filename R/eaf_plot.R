FRQ_DIFF_FCT = NULL

#' @title Effect allele frequency plot
#' @description
#' Plotting reported effect allele frequencies (EAF) against a reference set to identify
#' study variants which significantly deviate from the expected population frequencies.
#' @param gwas a data.table
#' @param eaf_col a string, the column containing the study EAF data
#' @param ref_eaf_col a string, the column containing the reference EAF data
#' @param tolerance a numeric, frequency difference that determines outliers
#' @param colours a 3 element list of colour codes, e.g. list(missing="#5B1A18", outlier="#FD6467", within="#7294D4")
#' @param title a string, the plot title
#' @param facet_grid_row_col (optional), a column by which to facet the plot by rows
#' @param facet_grid_col_col (optional), a column by which to facet the plot by columns
#' @return a ggplot
#' @export
#' @importFrom stats setNames
#' @importFrom data.table fcase
#' @import ggplot2
#'
eaf_plot <- function(gwas,
                     eaf_col = "EAF",
                     ref_eaf_col = "EUR_EAF",
                     tolerance = 0.2,
                     colours = list(missing="#5B1A18", outlier="#FD6467", within="#7294D4"),
                     title = NULL,
                     facet_grid_row_col = NULL,
                     facet_grid_col_col = NULL) {

  # checks
  stopifnot("gwas input must be a data.frame like object" = inherits(gwas, "data.frame"))
  stopifnot("eaf_col and ref_eaf_col must column names within the gwas data.frame" = all(c(eaf_col, ref_eaf_col) %in% colnames(gwas)))
  stopifnot("colours must be a list, length 3, of colour codes" = is.list(colours) & length(colours)==3 & all(sapply(colours, is.character)))

  # factor levels and colour labels
  diff_levels <- c("No reference data", "EAF outlier", "EAF within tolerance")
  colour_labels <- stats::setNames(c(colours[[1]], colours[[2]], colours[[3]]), diff_levels)

  # mutate a factor for the reference vs. study EAF difference
  gwas[, FRQ_DIFF_FCT := data.table::fcase(is.na(get(ref_eaf_col)), factor("No reference data", levels=diff_levels),
                                           abs(get(eaf_col)-get(ref_eaf_col)) > tolerance, factor("EAF outlier", levels=diff_levels),
                                           abs(get(eaf_col)-get(ref_eaf_col)) <= tolerance, factor("EAF within tolerance", levels=diff_levels),
                                           TRUE, factor(NA_character_)) ]

  # set no reference data to zero
  gwas[, (ref_eaf_col) := ifelse(is.na(get(ref_eaf_col)), 0, get(ref_eaf_col))]

  # plot
  plot <- ggplot2::ggplot(gwas, ggplot2::aes(x=get(ref_eaf_col), y=get(eaf_col), color=FRQ_DIFF_FCT)) +
    ggplot2::geom_point(size = 0.1) +
    ggplot2::geom_abline(slope=1, intercept= tolerance, color="darkgrey", linetype = "dotted") +
    ggplot2::geom_abline(slope=1, intercept=-tolerance, color="darkgrey", linetype = "dotted") +
    ggplot2::labs(title = title,
                  x = "Reference EAF",
                  y = "EAF",
                  color = NULL) +
    ggplot2::scale_color_manual(values=colour_labels) +
    ggplot2::theme_classic() +
    ggplot2::theme(aspect.ratio = 1,
                   legend.position = "top",
                   legend.box = "horizontal") +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 5))) +
    ggplot2::scale_y_continuous(limits=c(0,1), expand=c(0, 0)) +
    ggplot2::scale_x_continuous(limits=c(0,1), expand=c(0, 0))

  # facet plot if requested
  if(!is.null(facet_grid_row_col) & is.null(facet_grid_row_col)) {

    plot <- plot + ggplot2::facet_grid(rows = vars(get(facet_grid_row_col)))

  } else if(is.null(facet_grid_row_col) & !is.null(facet_grid_row_col)) {

    plot <- plot + ggplot2::facet_grid(cols = vars(get(facet_grid_col_col)))

  } else if(!all(is.null(facet_grid_row_col), is.null(facet_grid_row_col))) {

    plot <- plot + ggplot2::facet_grid(cols = vars(get(facet_grid_col_col)),
                                       rows = vars(get(facet_grid_row_col)))

  }

  # return the plot
  return(plot)
}

