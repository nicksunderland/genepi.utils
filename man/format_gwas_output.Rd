% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/standardise_gwas.R
\name{format_gwas_output}
\alias{format_gwas_output}
\title{Format a GWAS file}
\usage{
format_gwas_output(file_gwas, output_file, output_format = "default")
}
\arguments{
\item{file_gwas}{a string file path, or data.frame like object}

\item{output_file}{a string, a file path}

\item{output_format}{a string, a valid format - these can be obtained with \code{input_formats(ColumnMapping())}}
}
\value{
a file written to output_file
}
\description{
This formats a GWAS file column names and writes out to \code{output_file}
}
