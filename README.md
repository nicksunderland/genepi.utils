genepi.utils
================

<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->
<!-- badges: end -->

The `genepi.utils` package is a collection of utility functions for
working with genetic epidemiology data.

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
for successful installation on my Macbook M2 Max was creating a Makecars
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
remove.packages("fst")

install.packages("https://cran.r-project.org/src/contrib/data.table_1.14.8.tar.gz", type="source", repos=NULL)
install.packages("https://cran.r-project.org/src/contrib/furrr_0.3.1.tar.gz", type="source", repos=NULL)
install.packages("https://cran.r-project.org/src/contrib/fst_0.9.8.tar.gz", type="source", repos=NULL)
```

## Download the dbSNP data repository

Until this is hosted you will need to ask
<nicholas.sunderland@bristol.ac.uk> for the data files.

## A basic example

Here is an example of mapping CHR:POS:REF:ALT columns to dbSNP build 156
RSIDs.

``` r
library(genepi.utils)

gwas <- data.table::fread(system.file("extdata", "example_gwas_sumstats.tsv", package="genepi.utils"))

head(gwas)
#>    chromosome base_pair_location effect_allele other_allele
#>         <int>              <int>        <char>       <char>
#> 1:         10          100000625             G            A
#> 2:         10          100000645             C            A
#> 3:         10          100003242             G            T
#> 4:         10          100003785             C            T
#> 5:         10          100004360             A            G
```

``` r
gwas_with_rsids <- genepi.utils::chrpos_to_rsid(gwas,
                                                chr_col = "chromosome",
                                                pos_col = "base_pair_location",
                                                ea_col  = "effect_allele",
                                                nea_col = "other_allele",
                                                build   = "b37_dbsnp156")

head(gwas_with_rsids)

#>    chromosome base_pair_location effect_allele other_allele       RSID
#>         <int>              <int>        <char>       <char>     <char>
#> 1:         10          100000625             G            A  rs7899632
#> 2:         10          100000645             C            A rs61875309
#> 3:         10          100003242             G            T rs12258651
#> 4:         10          100003785             C            T  rs1359508
#> 5:         10          100004360             A            G  rs1048754
```

## Monitoring progress

The mapping can take some time, depending on the number of cores on your
machine. To monitor progress we have to use the form below. First we set
up the parallel processing plan with the number of workers - this should
be \<= the number of cores on your machine, which can be queried with
`parallel::detectCores()`. Then to monitor progress in the console we
need to wrap the function in `progressr::with_progress({ fn() }`.

``` r
future::plan(future::multisession, workers = 6)

progressr::with_progress({
  
  gwas_with_rsids <- genepi.utils::chrpos_to_rsid(gwas,
                                                  chr_col = "chromosome",
                                                  pos_col = "base_pair_location",
                                                  ea_col  = "effect_allele",
                                                  nea_col = "other_allele",
                                                  build   = "b37_dbsnp156")

})

# Chromosomes to process:  10
# |=================                                                                    |  12%
```

## Function options

There are multiple options that can lead to significant changes in the
processing time.

### Alt RSID output

In the dbSNP database there are multiple RSID assigned to the same
CHR:POS:REF:ALT combination. This function will return the first RSID
encountered and filter out the rest. If you want to get the alternative
RSIDs then use the `alt_rsids=TRUE` flag. This will change the structure
of the output to a named `list()` of two `data.table` elements: 1)
“data”=gwas_data and 2) “alt_rsids”=alt_rsid_data.

``` r
gwas_with_rsids <- genepi.utils::chrpos_to_rsid(gwas,
                                                chr_col = "chromosome",
                                                pos_col = "base_pair_location",
                                                ea_col  = "effect_allele",
                                                nea_col = "other_allele",
                                                build   = "b37_dbsnp156", 
                                                alt_rsids=TRUE)
#> Chromosomes to process:  10

str(gwas_with_rsids)
#> List of 2
#>  $ data     :Classes 'data.table' and 'data.frame':  5 obs. of  5 variables:
#>   ..$ chromosome        : int [1:5] 10 10 10 10 10
#>   ..$ base_pair_location: int [1:5] 100000625 100000645 100003242 100003785 100004360
#>   ..$ effect_allele     : chr [1:5] "G" "C" "G" "C" ...
#>   ..$ other_allele      : chr [1:5] "A" "A" "T" "T" ...
#>   ..$ RSID              : chr [1:5] "rs7899632" "rs61875309" "rs12258651" "rs1359508" ...
#>   ..- attr(*, ".internal.selfref")=<externalptr> 
#>  $ alt_rsids:Classes 'data.table' and 'data.frame':  0 obs. of  6 variables:
#>   ..$ RSID    : chr(0) 
#>   ..$ CHR     : chr(0) 
#>   ..$ BP      : int(0) 
#>   ..$ REF     : chr(0) 
#>   ..$ ALT     : chr(0) 
#>   ..$ baseRSID: chr(0) 
#>   ..- attr(*, ".internal.selfref")=<externalptr>
```

### Allele specification

Although matching on alleles is desirable, especially for indels, we can
still match on just chromosome and position by either leaving `ea_col`
and `nea_col` arguments out, or setting them to `NULL`.

``` r
gwas_with_rsids <- genepi.utils::chrpos_to_rsid(gwas,
                                                chr_col = "chromosome",
                                                pos_col = "base_pair_location",
                                                build   = "b37_dbsnp156")
#> Chromosomes to process:  10

str(gwas_with_rsids)
#> Classes 'data.table' and 'data.frame':   5 obs. of  5 variables:
#>  $ chromosome        : chr  "10" "10" "10" "10" ...
#>  $ base_pair_location: int  100000625 100000645 100003242 100003785 100004360
#>  $ effect_allele     : chr  "G" "C" "G" "C" ...
#>  $ other_allele      : chr  "A" "A" "T" "T" ...
#>  $ RSID              : chr  "rs7899632" "rs61875309" "rs12258651" "rs1359508" ...
#>  - attr(*, ".internal.selfref")=<externalptr>
```

## Evaluation speed

The choice of parameters will impact computation speed. Since the dbSNP
data is stored as `.fst` binary files and only the desired rows /
columns are ever read into memory, the more data that is requested the
slower the computation. That said, it is still much faster / feasible
than trying to read in the entire dbSNP database (over 1 billion rsIDs).

<div class="figure" style="text-align: center">

<img src="vignettes/figures/microbenchmark_chrpos_to_rsid.png" alt="Computation speed: Apple M2 Max 96GB 10 cores; 4.3 million SNPs queried against ~1 billion dbSNP156 RSIDs" width="70%" />
<p class="caption">
Computation speed: Apple M2 Max 96GB 10 cores; 4.3 million SNPs queried
against ~1 billion dbSNP156 RSIDs
</p>

</div>

``` r
library(microbenchmark)
library(ggplot2)

# some GWAS data
dt <- data.table::fread( gwas_sumstats_4.3million_rows )

# benchmarking
mbm <- microbenchmark("chrpos_to_rsid: no alleles, no alt rsids" = {
                              future::plan(future::multisession, workers = 10)
                              progressr::with_progress({
                                genepi.utils::chrpos_to_rsid(dt,
                                                             chr_col="CHR",
                                                             pos_col="POS",
                                                             ea_col=NULL,
                                                             nea_col=NULL,
                                                             build="b37_dbsnp156",
                                                             alt_rsids=FALSE,
                                                             flip="report")
                              })
                        },
                      "chrpos_to_rsid: no alleles, with alt rsids" = {
                        future::plan(future::multisession, workers = 10)
                        progressr::with_progress({
                          genepi.utils::chrpos_to_rsid(dt,
                                                       chr_col="CHR",
                                                       pos_col="POS",
                                                       ea_col=NULL,
                                                       nea_col=NULL,
                                                       build="b37_dbsnp156",
                                                       alt_rsids=TRUE,
                                                       flip="report")
                        })
                      },
                      "chrpos_to_rsid: with alleles, no alt rsids" = {
                        future::plan(future::multisession, workers = 10)
                        progressr::with_progress({
                          genepi.utils::chrpos_to_rsid(dt,
                                                       chr_col="CHR",
                                                       pos_col="POS",
                                                       ea_col="EFFECT_ALLELE",
                                                       nea_col="OTHER_ALLELE",
                                                       build="b37_dbsnp156",
                                                       alt_rsids=FALSE,
                                                       flip="report")
                        })
                      },
                      "chrpos_to_rsid: with alleles, with alt rsids" = {
                        future::plan(future::multisession, workers = 10)
                        progressr::with_progress({
                          genepi.utils::chrpos_to_rsid(dt,
                                                       chr_col="CHR",
                                                       pos_col="POS",
                                                       ea_col="EFFECT_ALLELE",
                                                       nea_col="OTHER_ALLELE",
                                                       build="b37_dbsnp156",
                                                       alt_rsids=TRUE,
                                                       flip="report")
                        })
                      },
                      times = 10L)

# plot
autoplot(mbm)
```
