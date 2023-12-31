% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/collider_bias.R
\name{apply_collider_correction}
\alias{apply_collider_correction}
\title{Apply correction factor to GWAS data}
\usage{
apply_collider_correction(
  gwas_i,
  gwas_p,
  b_correction_factor,
  b_std_err,
  merge = c(CHR = "CHR", BP = "BP"),
  keep_palindromic = FALSE
)
}
\arguments{
\item{gwas_i}{a data.table, incidence GWAS}

\item{gwas_p}{a data.table, progression GWAS}

\item{b_correction_factor}{a numeric, the correction factor to apply}

\item{b_std_err}{a numeric, the correction factor standard error to apply}

\item{merge}{a named character vector c(gwas1_col=gwas2_col) of the columns to join on; e.g. c("CHR1"="CHR2","BP1"="BP2"); can be NULL, in which case the gwases must be data.tables and already have keys set}

\item{keep_palindromic}{a logical, whether to allow palindromic alleles to be adjusted and returned}
}
\value{
a data.table, progression GWAS with additional columns c("adjusted_beta","adjusted_se","adjusted_p")
}
\description{
Adjust progression GWAS data for collider bias using a calculated correction
factor and correction factor standard error. The correction factor should be
calculated with one of the Slope-hunter, IVW-MR, or Dudbridge methods. \cr
This function first harmonises the datasets and removes invalid and palindromic
alleles as default
}
