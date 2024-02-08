genepi.utils
================

<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->
<!-- badges: end -->

The `genepi.utils` package is a collection of utility functions for
working with genetic epidemiology data. For common use cases please see
the [vignettes](https://nicksunderland.github.io/genepi.utils/).

## Installation

You can install the development version of genepi.utils like so:

``` r
# install.packages("devtools")
devtools::install_github("nicksunderland/genepi.utils")
```

The package relies on the parallel processing capabilities of the
[`data.table`](https://rdatatable.gitlab.io/data.table/),
[`furrr`](https://furrr.futureverse.org) and
[`fst`](https://www.fstpackage.org) packages. These can be installed
like so:

``` r
install.packages("data.table")
install.packages("furrr")
install.packages("fst")
```

**Important**: on mac it can be more challenging to enable OpenMP
parallel processing as the *clang* compiler does not include an OpenMP
runtime as standard. I recommend following the instructions on the
`data.table` GitHub
[link](https://github.com/Rdatatable/data.table/wiki/Installation). Key
for successful installation on my Macbook M2 Max was creating a Makevars
file in the root directory `~.R/Makevars`, which is a simple text file,
containing the compilation flags below, and then re-installing from
source.

``` bash
LDFLAGS += -L/opt/homebrew/opt/libomp/lib -lomp
CPPFLAGS += -I/opt/homebrew/opt/libomp/include -Xclang -fopenmp
```

Re-install packages from source.

``` r
remove.packages("data.table")
remove.packages("furrr")
remove.packages("fstcore")
remove.packages("fst")

install.packages("https://cran.r-project.org/src/contrib/data.table_1.14.10.tar.gz", type="source", repos=NULL)
install.packages("https://cran.r-project.org/src/contrib/furrr_0.3.1.tar.gz", type="source", repos=NULL)
install.packages("https://cran.r-project.org/src/contrib/fstcore_0.9.18.tar.gz", type="source", repos=NULL)
install.packages("https://cran.r-project.org/src/contrib/fst_0.9.8.tar.gz", type="source", repos=NULL)
```

## Download the dbSNP data repository

For RSID mapping you will need a copy of the dbSNP `.fst` file
directory. Until this is hosted somewhere you will need to ask
<nicholas.sunderland@bristol.ac.uk> for the data files. If you are
working at the Bristol IEU the package should work on BC4. If not, ask
for the location of the dbSNP `.fst` file directory on the HPC. If you
use a custom location for the dbSNP directory you will need to set this
in the package - this only needs to be done once per package install.

``` r
set_dbsnp_directory("/path_to_directory/dbsnp")
```
