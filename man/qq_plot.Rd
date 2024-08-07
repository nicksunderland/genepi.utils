% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qq_plot.R
\name{qq_plot}
\alias{qq_plot}
\title{QQ plot}
\usage{
qq_plot(
  gwas,
  pval_col = "p",
  colours = list(raw = "#2166AC"),
  title = NULL,
  subtitle = NULL,
  plot_corrected = FALSE,
  facet_grid_row_col = NULL,
  facet_grid_col_col = NULL,
  facet_nrow = NULL,
  facet_ncol = NULL
)
}
\arguments{
\item{gwas}{a data.frame like object or valid file path}

\item{pval_col}{the P value column}

\item{colours}{a 2 element list of colour codes (1-the uncorrected points, 2-the GC corrected points)}

\item{title}{a string, the title for the plot}

\item{subtitle}{a string, the subtitle for the plot}

\item{plot_corrected}{a logical, whether to apply and plot the lambda correction}

\item{facet_grid_row_col}{a string, the column name in \code{gwas} by which to facet the plot (rows)}

\item{facet_grid_col_col}{a string, the column name in \code{gwas} by which to facet the plot (cols)}

\item{facet_nrow}{an integer, passed to facet_wrap, the number of rows to facet by (if only facet_grid_row_col is provided)}

\item{facet_ncol}{an integer, passed to facet_wrap, the number of cols to facet by (if only facet_grid_col_col is provided)}
}
\value{
a ggplot
}
\description{
QQ plot
}
