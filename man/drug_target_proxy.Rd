% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drug_target_proxy.R
\name{drug_target_proxy}
\alias{drug_target_proxy}
\title{Create a drug target proxy instrument}
\usage{
drug_target_proxy(
  gwas_gene,
  gene_chr,
  gene_start,
  gene_end,
  gene_flanks = 5e+05,
  build = "GRCh37",
  clump = TRUE,
  clump_ref = which_1000G_reference("GRCh37"),
  p1 = 5e-08,
  p2 = 1,
  r2 = 0.2,
  kb = 250,
  join_key = "RSID",
  QTL_list = list(),
  concordance = data.table::data.table(data_name_1 = character(), data_name_2 =
    character(), concordant = logical())
)
}
\arguments{
\item{gwas_gene}{.}

\item{gene_chr}{.}

\item{gene_start}{.}

\item{gene_end}{.}

\item{gene_flanks}{.}

\item{build}{.}

\item{clump}{.}

\item{clump_ref}{.}

\item{p1}{.}

\item{p2}{.}

\item{r2}{.}

\item{kb}{.}

\item{join_key}{.}

\item{QTL_list}{.}

\item{concordance}{.}
}
\value{
a data.table
}
\description{
Create a drug target proxy instrument base on several association statistics
}
