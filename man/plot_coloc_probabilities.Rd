% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/colocalisation.R
\name{plot_coloc_probabilities}
\alias{plot_coloc_probabilities}
\title{Coloc probability plot}
\usage{
plot_coloc_probabilities(coloc, rule = "H4 > 0.5", type = "prior", row = 1)
}
\arguments{
\item{coloc}{coloc object, output from \code{coloc::coloc.abf()}}

\item{rule}{a string, a valid rule indicating success e.g. "H4 > 0.5"}

\item{type}{a string, either \code{prior} or \code{posterior}}

\item{row}{an integer, row in a \code{coloc.susie} or \code{coloc.signals} object}
}
\value{
a ggplot
}
\description{
A plotting wrapper for the \code{coloc} package. Produces a ggplot for either
the prior or posterior probability sensitivity analyses. See the
\href{https://chr1swallace.github.io/coloc/articles/a04_sensitivity.html}{coloc}
package vignettes for details.
}
\references{
\href{https://chr1swallace.github.io/coloc/articles/a04_sensitivity.html}{coloc}
}
