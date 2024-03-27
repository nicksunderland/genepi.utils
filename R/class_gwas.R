# Silence R CMD check
globalVariables(c("p", "double.xmin", "n", "ncase", "trait", "id", "rsid", "chr", "bp", "ea", "oa", "eaf", "se", "correlation", "strand", "imputed", "info"),
                package = "genepi.utils")

#' @title GWAS object
#' @description
#' A GWAS object is a container for vectors of GWAS data, a correlation matrix, and
#' meta-data regarding quality control procedures applied at the point of object creation / data import. \cr\cr
#' \strong{Quality control steps:}
#' \tabular{rr}{
#'   \strong{Step} \tab \strong{Description} \cr
#'   Column mapping \tab Column name mapping using the `ColumnMap` and `Column` classes which define standard column names and a list
#'   of aliases. The set of core data fields are: `rsid`, `chr`, `bp`, `ea`, `oa`, `eaf`, `beta`, `se`, `p`, and `n`. \cr
#'   Column typing \tab Column types are enforced using the `ColumnMap` and `Column` classes which define standard column types (see
#'   return value table below). \cr
#'   Standardise columns \tab Core data fields are checked, if absent, they are calculated from the available data if possible. The
#'   current conversions include: \emph{chr} and \emph{bp} extration from \emph{rsid} column containing \emph{chr:pos +/- allele} id coding;
#'   conversion of odds ratio \emph{or} and odds ratio standard error / upper and lower bounds \emph{or_se/or_ub/or_lb} to \emph{beta} and
#'   \emph{se}; and conversion of \emph{z} scores to \emph{p} values. \cr
#'   Standardise allele coding \tab Alleles are standardised to upper case \emph{ACTGDI} or combination of these letters. `NA` values
#'   are allowed in one of the \emph{ea} or \emph{oa} fields, but not both, where it is taken to represent a deletion. \cr
#'   Reference harmonisation \tab TODO: harmonising against a reference panel is an important step, however has not yet been implemented
#'   for this class constructor (under active development). \cr
#'   Filters \tab Custom filters can be added to the `filters` argument. This should be a list of named string that can be evaluated as R
#'   expressions and when applied to the standard column names result in a logical vector mask of the rows to \emph{keep}. The defaults are:
#'   \itemize{
#'      \item{beta_invalid} = {"!is.infinite(beta) & abs(beta) < 20"}
#'      \item{eaf_invalid} = {"eaf > 0 & eaf < 1"}
#'      \item{p_invalid} = {"!is.infinite(p)"}
#'      \item{se_invalid} = {"!is.infinite(se)"}
#'      \item{chr_missing} = {"!is.na(chr)"}
#'      \item{bp_missing} = {"!is.na(bp)"}
#'      \item{beta_missing} = {"!is.na(beta)"}
#'      \item{se_missing} = {"!is.na(se)"}
#'      \item{p_missing} = {"!is.na(p)"}
#'      \item{eaf_missing} = {"!is.na(eaf)"}
#'   }
#' }
#' @param dat a valid string file path to be read by `data.table::fread` or a `data.table::data.table` object; the GWAS data source
#' @param map a valid input to the `ColumnMap` class constructor (a predefined map string id, a named list or character vector, or a ColumMap object)
#' @param drop a logical, whether to drop data source columns not in the column `map`
#' @param fill a logical, whether to add (NAs) missing columns present in the column `map` but not present in the data source
#' @param fill_rsid either FALSE or a valid argument for the chrpos_to_rsid `build` argument, e.g. "b37_dbsnp156"
#' @param missing_rsid a string, how to handle missing rsids: one of "fill_CHR:BP", "fill_CHR:BP_OA_EA", "overwrite_CHR:BP",
#' "overwrite_CHR:BP:OA:EA", "none", or "leave"
#' @param parallel_cores an integer, number of cores to used for RSID mapping, default is maximum machine cores
#' @param filters a list of named strings, each to be evaluated as an expression to filter the data during the quality control steps (above)
#' @param reference a valid string file path to be read by `data.table::fread` or a `data.table::data.table` object; the reference data
#' @param ref_map a valid input to the `ColumnMap` class constructor (a predefined map id (a string), a named list or character vector, or a
#' ColumMap object) defining at least columns `rsid` (or `chr`, `bp`), `ea`, `oa` and `eaf`.
#' @param verbose a logical, whether to print details
#' @param ... variable capture to be passed to the constructor, e.g. individual vectors for the slots, rather that `dat`
#'
#' @slot rsid character, variant ID - usually in rs12345 format, however this can be changed with the `missing_rsid` argument
#' @slot chr character, chromosome identifier
#' @slot bp integer, base position
#' @slot ea character, effect allele
#' @slot oa character, other allele
#' @slot eaf numeric, effect allele frequency
#' @slot beta numeric, effect size
#' @slot se numeric, effect size standard error
#' @slot p numeric, p-value
#' @slot n integer, total number of samples
#' @slot ncase integer, number of cases
#' @slot trait character, the GWAS trait
#' @slot id character, the GWAS identifier
#' @slot source character, data source; either the file path, or "data.table" if loaded directly
#' @slot correlation matrix, a correlation matrix of signed R values between variants
#' @slot map `ColumnMap`, a mapping of class ColumnMap
#' @slot qc list, a named list of filters; name is the filter expression and value is an integer vector of rows that fail the filter
#'
#' @return an S7 class genepi.utils::GWAS object
#'
#' @import S7
#' @include class_colmap.R
#' @export
GWAS <- new_class(
  #==============================
  # GWAS class name
  #==============================
  name    = "GWAS",
  package = "genepi.utils",

  #==============================
  # GWAS class properties
  #==============================
  properties = list(
    #----------------------------
    # column name mapping
    #----------------------------
    map     = ColumnMap,
    #----------------------------
    # qc procedures
    #----------------------------
    qc      = class_list,
    #----------------------------
    # column names / data vectors
    #----------------------------
    rsid    = class_character,
    chr     = class_character,
    bp      = class_integer,
    ea      = class_character,
    oa      = class_character,
    eaf     = class_numeric,
    beta    = class_numeric,
    se      = class_numeric,
    p       = class_numeric,
    n       = class_integer,
    ncase   = class_integer,
    strand  = class_character,
    imputed = class_logical,
    info    = class_numeric,
    #----------------------------
    # correlation matrix
    #----------------------------
    correlation = new_S3_class('matrix'),
    #----------------------------
    # meta data
    #----------------------------
    trait   = class_character,
    id      = class_character,
    source  = class_character
  ),

  #==============================
  # GWAS class constructor func.
  #==============================
  constructor = function(dat,
                         map            = "default",
                         drop           = FALSE,
                         fill           = FALSE,
                         fill_rsid      = FALSE,
                         missing_rsid   = "fill_CHR:BP",
                         parallel_cores = parallel::detectCores(),
                         filters = list(
                           beta_invalid    = "!is.infinite(beta) & abs(beta) < 20",
                           eaf_invalid     = "eaf > 0 & eaf < 1",
                           p_invalid       = "!is.infinite(p)",
                           se_invalid      = "!is.infinite(se)",
                           alleles_invalid = "!is.na(ea) & !is.na(oa)",
                           chr_missing     = "!is.na(chr)",
                           bp_missing      = "!is.na(bp)",
                           beta_missing    = "!is.na(beta)",
                           se_missing      = "!is.na(se)",
                           p_missing       = "!is.na(p)",
                           eaf_missing     = "!is.na(eaf)"),
                         reference      = NULL,
                         ref_map        = NULL,
                         verbose        = TRUE,
                         ...) {

    if(verbose) message("Loading GWAS ...")

    # parse the map
    map <- ColumnMap(map)

    # as data table
    if(inherits(dat, 'data.frame') && !inherits(dat, 'data.table')) { dat <- data.table::as.data.table(dat) }

    # load (either from file or data.frame)
    gwas <- load_gwas(dat, map=map, drop=drop, fill=fill, verbose=verbose)

    # ensure correct type
    gwas <- type_columns(gwas, map=map, verbose=verbose)

    # standardise columns
    gwas <- standardise_columns(gwas, verbose=verbose)

    # standardise alleles
    gwas <- standardise_alleles(gwas, verbose=verbose)

    # apply filters
    g    <- apply_filters(gwas, filters=filters, verbose=verbose)
    gwas <- g$gwas
    qc   <- g$qc

    # add or reformat id/rsid
    gwas <- populate_rsid(gwas, fill_rsid=fill_rsid, missing_rsid=missing_rsid, parallel_cores=parallel_cores, verbose=verbose)

    # harmonise to reference
    g <- harmonise_reference(gwas, reference=reference, ref_map=ref_map, qc=qc, verbose=verbose)
    gwas <- g$gwas
    qc   <- g$qc

    # other changes
    gwas[p==0, p := .Machine$double.xmin]

    # group together as a big list; add the other properties that may have been passed
    props <- c(gwas, list(map=map, ...))

    # don't need to save all of repetitive columns
    if(!"n"      %in% names(props))                                   { props$n       <- NA_integer_         }
    if(!"ncase"  %in% names(props))                                   { props$ncase   <- NA_integer_         }
    if(!"strand" %in% names(props))                                   { props$strand  <- NA_character_       }
    if(!"imputed"%in% names(props))                                   { props$imputed <- NA                  }
    if(!"info"   %in% names(props))                                   { props$info    <- NA_real_            }
    if("n"       %in% names(props) && length(unique(props$n))==1    ) { props$n       <- unique(props$n)     }
    if("ncase"   %in% names(props) && length(unique(props$ncase))==1) { props$ncase   <- unique(props$ncase) }
    if("strand"  %in% names(props) && length(unique(props$strand))==1){ props$strand  <- unique(props$strand)}
    if(!"trait"  %in% names(props))                                   { props$trait   <- "trait"}
    if(!"id"     %in% names(props))                                   { props$id      <- "id" }
    if(length(unique(props$trait))==1)                                { props$trait   <- unique(props$trait) }
    if(length(unique(props$id))==1)                                   { props$id      <- unique(props$id)    }

    # set filters and qc
    props$qc <- qc

    # set correlation matrix if not provided
    if(!"correlation" %in% names(props)) { props$correlation <- matrix() }

    # set data source
    if(is.character(dat) && file.exists(dat)) {
      props$source <- dat
    } else {
      props$source <- 'data.table'
    }

    # assign to the class object
    object <- new_object(S7::S7_object(),
               map         = props$map,
               qc          = props$qc,
               rsid        = props$rsid,
               chr         = props$chr,
               bp          = props$bp,
               ea          = props$ea,
               oa          = props$oa,
               eaf         = props$eaf,
               beta        = props$beta,
               se          = props$se,
               p           = props$p,
               n           = props$n,
               ncase       = props$ncase,
               strand      = props$strand,
               imputed     = props$imputed,
               info        = props$info,
               trait       = props$trait,
               id          = props$id,
               source      = props$source,
               correlation = props$correlation)

    # return the object
    return(object)
  },

  #==============================
  # GWAS class validator func.
  #==============================
  validator = function(self) {
    stopifnot("Unequal vector lengths" = sapply(lengths(list(self@rsid, self@chr, self@bp, self@ea, self@oa, self@eaf, self@beta, self@se, self@p)), function(x) x == length(self@rsid)))
    stopifnot("Invalid sample size `n`"       = length(self@n)<=1 || length(self@n)==length(self@rsid))
    stopifnot("Invalid sample size `ncase`"   = length(self@ncase)<=1 || length(self@ncase)==length(self@rsid))
    stopifnot("Invalid sample size `strand`"  = length(self@strand)<=1 || length(self@strand)==length(self@rsid))
    stopifnot("Invalid sample size `imputed`" = length(self@imputed)<=1 || length(self@imputed)==length(self@rsid))
    stopifnot("Invalid sample size `info`"    = length(self@info)<=1 || length(self@info)==length(self@rsid))
    stopifnot("Invalid `trait` field length"  = length(self@trait)<=1)
    stopifnot("Invalid `id` field length"     = length(self@id)<=1)
  }
)



load_gwas <- new_generic('load_gwas', c('dat','map'), function(dat, map, drop=FALSE, fill=FALSE, verbose=TRUE) { S7_dispatch() })
method(load_gwas, list(class_character, ColumnMap)) <- function(dat, map, drop=FALSE, fill=FALSE, verbose=TRUE) {
  stopifnot("If `dat` is a character() it must be a valid file path" = file.exists(dat))
  if(verbose) message("[i] loading data from: ", basename(dat))

  h <- data.table::fread(dat, nrows=1, nThread=parallel::detectCores())
  n <- raw_names(map, names(h))
  n <- n[!is.na(n)]
  d <- data.table::fread(dat, select=unname(n), nThread=parallel::detectCores())
  d <- load_gwas(dat=d, map=map, drop=drop, fill=fill, verbose=verbose)

  return(d)
}
method(load_gwas, list(new_S3_class('data.table'), ColumnMap)) <- function(dat, map, drop=FALSE, fill=FALSE, verbose=TRUE) {

  if(verbose) message("[i] applying column mapping")

  n    <- raw_names(map, names(dat))
  here <- map[!is.na(n)]
  miss <- map[ is.na(n)]

  # rename to standard names
  data.table::setnames(dat, unname(n[!is.na(n)]), names(n[!is.na(n)]))

  # add requested columns if not present
  stopifnot("Mapping column name(s) not found in data.table" = fill || length(miss@map)==0)
  if(fill && length(miss@map) > 0) {
    if(verbose) {
      message(paste0("\t[+] adding missing column(s) [", paste0(sapply(miss@map, function(x) x@name), collapse=", "), "] specified in mapping but not found in gwas data  (see `fill` option)"))
    }
    dat[, sapply(miss@map, function(x) x@name) := lapply(miss@map, function(x) vector(mode=x@type, length=nrow(dat)))] # force col type
    dat[, sapply(miss@map, function(x) x@name) := NA] # set to NA
  }

  # drop columns if requested
  if(drop && length(names(dat)[!names(dat) %in% names(n)]) > 0) {
    if(verbose) {
      message(paste0("\t[-] dropping unwanted column(s) [", paste0(names(dat)[!names(dat) %in% names(n)], collapse=", "), "] found in gwas data but not specified in mapping (see `drop` option)"))
    }
    dat[, names(dat)[!names(dat) %in% names(n)] := NULL]
  }

  return(dat)
}


standardise_columns <- new_generic('standardise_columns', 'gwas', function(gwas, verbose=TRUE) { S7_dispatch() })
method(standardise_columns, new_S3_class('data.table')) <- function(gwas, verbose=TRUE) {

  if(verbose) message("[i] standardising columns")

  # if chromosome and position column not present, try to parse from rsid column
  if (!all(c("chr", "bp") %in% names(gwas))) {

    if(verbose) message("\t[i] attempting to parse `chr` and `bp` columns from `rsid` column")

    valid_pattern <- grepl("(?i)(\\d+|X|Y)[:_]\\d+", gwas$rsid)
    if(all(valid_pattern)) {
      gwas[, chr := sub("(?i)(\\d+|X|Y)[:_]\\d+(?:.*)", "\\1", rsid)]
      gwas[, bp  := sub("(?i)(?:\\d+|X|Y)[:_](\\d+)(?:.*)", "\\1", rsid)]
    } else {
      example_failure <- gwas[[which(!valid_pattern)[[1]], "rsid"]]
      stop(paste0("problem standardising `chr` and `bp` columns. Cannot parse rsid` column with regex \\d[:_]\\d for enteries such as: ", example_failure))
    }
  }

  # check chromosome coding
  invalid_chr <- which(gwas[, !chr %in% as.character(1:25)])
  if (length(invalid_chr) > 0) {

    message(paste0("\t[i] ", length( invalid_chr ), " invalid chromosome codings e.g. '", gwas[invalid_chr[1], chr], "', attempting recoding..."))

    # recode chromosomes character integers 1:25
    gwas[invalid_chr, chr := data.table::fcase(grepl("X|PAR1", chr, ignore.case = TRUE), "23",
                                               grepl("Y|PAR2", chr, ignore.case = TRUE), "24",
                                               grepl("M|MT",   chr, ignore.case = TRUE), "25",
                                               default = NA_character_)]

    # recheck
    invalid_chr <- which(gwas[, !chr %in% as.character(1:25)])
    if (length(invalid_chr) > 0) {

      message(paste0("\t[i] ", length( invalid_chr ), " chromosome codings still invalid after recoding to 1:25 e.g. '", gwas[invalid_chr[1], chr], "'"))

    }

  }

  # convert odds ratio to betas and standard error, if beta, se not provided
  if(!all(c("beta","se") %in% names(gwas))) {

    if(verbose) message("\t[i] attempting to covert `or` column to `beta` and `se` columns")

    if ("or" %in% names(gwas) && (all(c("or_ub", "or_lb") %in% names(gwas)) | "or_se" %in% names(gwas))) {
      gwas[, beta := log(or)]
      if (all(c("or_ub", "or_lb") %in% names(gwas)) && !"or_se" %in% names(gwas)) {
        gwas[, or_ub := or + or_se * 1.96]
        gwas[, or_lb := or - or_se * 1.96]
      }
      gwas[, se := (or_ub - or_lb) / (2 * 1.95996)]
    } else {
      stop("problem standardising/creating `beta` and `se` columns. No `beta` `se` columns provided and no `or` `or_lb` `or_ub` or `or` `or_se` columns provided to compute them.")
    }
  }

  # convert Z-score to Pvalue if P not provided but Z score is
  if(!"p" %in% names(gwas)) {
    if("z" %in% names(gwas)) {

      if(verbose) message("\t[i] coverting `z` column to `p` column")

      gwas[, p := 2 * stats::pnorm(-abs(z))]
    } else {
      stop("problem standardising/creating `p` value column. No `p` value column found and no `z` column provided to compute it.")
    }
  }

  return(gwas)
}


type_columns <- new_generic('type_columns', c('gwas','map'), function(gwas, map, verbose=TRUE) { S7_dispatch() })
method(type_columns, list(new_S3_class('data.table'), ColumnMap)) <- function(gwas, map, verbose=TRUE) {

  if(verbose) message("[i] enforcing column types")

  here <- map[sapply(map@map, function(x) x@name) %in% names(gwas)]

  for(j in 1:ncol(gwas)) {
    col_name <- names(gwas)[j]
    if(typeof(gwas[, get(col_name)])==here[[col_name]]@type) next
    gwas[, (col_name) := do.call(paste0('as.',here[[col_name]]@type), list(x=get(col_name)))]
  }

  return(gwas)
}


apply_filters <- new_generic('apply_filters', c('gwas','filters'), function(gwas, filters, remove=FALSE, verbose=TRUE) { S7_dispatch() })
method(apply_filters, list(new_S3_class('data.table'), class_list)) <- function(gwas, filters, remove=FALSE, verbose=TRUE) {

  if(verbose) message("[i] applying filters")

  if(length(filters) == 0) return(list(gwas=gwas, qc=filters))

  if(is.null(names(filters))) { names(filters) <- paste0('filter_', 1:length(filters)) }

  # data.table to store QC results
  qc <- list()

  # assess each filter
  for(i in seq_along(filters)) {

    # parse the string / filter expression
    filter_expr <- parse(text = filters[[i]])

    # apply filter / create row index logical mask
    qc[[ filters[[i]] ]] <- which(!gwas[, eval(filter_expr)])

    if(verbose && length( qc[[ filters[[i]] ]] )) {
      message(paste0("\t[-] ", length( qc[[ filters[[i]] ]] ), " row warning(s) with filter '", names(filters)[i], "' [", filter_expr, "]" ))
    }

  }

  # add filter fail/pass
  qc[[ 'fail' ]] <- unique(unlist(qc))
  if(remove && length(qc[['fail']]) > 0) {
    gwas <- gwas[-qc[['fail']], ]
  }

  return(list(gwas=gwas, qc=qc))
}


standardise_alleles <- new_generic('standardise_alleles', 'gwas', function(gwas, verbose=TRUE) { S7_dispatch() })
method(standardise_alleles, new_S3_class('data.table')) <- function(gwas, verbose=TRUE) {

  if(verbose) message("[i] standardising allele coding")

  # make uppercase
  gwas[, ea := toupper(ea)]
  gwas[, oa := toupper(oa)]

  # define valid allele formats - NA allowed as deletion
  valid_alleles <- (((grepl("^([AGCT]+|[D])$", gwas[["ea"]]) | is.na(gwas[["ea"]])) & grepl("^([AGCT]+|[I])$", gwas[["oa"]])) |   # positive assertion - EA valid SNP or deleteion, OA valid SNP or insertion
                    ((grepl("^([AGCT]+|[D])$", gwas[["oa"]]) | is.na(gwas[["oa"]])) & grepl("^([AGCT]+|[I])$", gwas[["ea"]]))) &  # positive assertion - OA valid SNP or deleteion, EA valid SNP or insertion
                    (!gwas[["ea"]]==gwas[["oa"]])                                                                                  # negative assertion - OA can't be the same as EA

  # invalid alleles to NA
  gwas[!valid_alleles, c("ea", "oa") := list(NA_character_, NA_character_)]

  return(gwas)
}



harmonise_reference <- new_generic('harmonise_reference', 'gwas', function(gwas, reference, ref_map, qc, verbose=TRUE) { S7_dispatch() })
method(harmonise_reference, new_S3_class('data.table')) <- function(gwas, reference, ref_map, qc, verbose=TRUE) {

  # if null just return
  if (is.null(reference) || is.null(ref_map)) return(list(gwas=gwas, qc=qc))

  stop("harmonise_reference not implemented yet")

#   # ensure base reference columns present in mapping
#   req_cols <- c("rsid", "chr", "bp", "ea", "oa", "eaf")
#   req_map  <- sapply(ref_map@map, function(x) x@name %in% req_cols)
#   ref_map  <- ref_map[req_map]
#   map_names<- sapply(ref_map@map, function(x) x@name)
#   stopifnot("Reference map must contain `rsid`, `chr`, `bp`, `ea`, `oa`, `eaf`" = all(sapply(req_cols, function(x) x %in% map_names)))
#
#   # read the reference header and the mapping, ensure base cols in data
#   h <- data.table::fread(reference, nrows=1, nThread=parallel::detectCores())
#   n <- raw_names(ref_map, names(h))
#   n <- n[names(n) %in% req_cols & !is.na(n)]
#   stopifnot("Reference data must contain columns `rsid`, `chr`, `bp`, `ea`, `oa`, `eaf`" = all(sapply(req_cols, function(x) x %in% names(n))))
#
#   # read in all of the reference data
#   r <- data.table::fread(reference, select=unname(n), nThread=parallel::detectCores())
#   r <- load_gwas(dat=r, map=ref_map, drop=TRUE, fill=FALSE, verbose=verbose)
#
#   # enforce types
#   r <- type_columns(r, ref_map, verbose=verbose)
#
#   # check chromosome coding
#   if (any(!r$chr %in% as.character(1:25))) {
#     warning("Non-standard chromosome coding (character '1'-'25') in reference file, please change.")
#   }
#
#   # make all of the non-rsid IDs chr:pos
#   gwas[grepl("^rs[0-9]+$", rsid), rsid := paste0(chr, ":", bp)]
#   r[grepl("^rs[0-9]+$", rsid), rsid := paste0(chr, ":", bp)]
#   r[, c("beta", "se", "p") := list(-9, -9, -9)]
#
#   # format for TSMR
#   r_formatted <- TwoSampleMR::format_data(
#     dat               = r[1:10,] |> as.data.frame(),
#     type              = "outcome",
#     snp_col           = "rsid",
#     chr_col           = "chr",
#     pos_col           = "bp",
#     effect_allele_col = "ea",
#     other_allele_col  = "oa",
#     beta_col          = "beta",
#     eaf_col           = "eaf",
#     se_col            = "se",
#     pval_col          = "p") |> data.table::as.data.table()
#
#   g_formatted <- TwoSampleMR::format_data(
#     dat               = gwas[1:10,] |> as.data.frame(),
#     type              = "exposure",
#     snp_col           = "rsid",
#     chr_col           = "chr",
#     pos_col           = "bp",
#     effect_allele_col = "ea",
#     other_allele_col  = "oa",
#     eaf_col           = "eaf",
#     beta_col          = "beta",
#     se_col            = "se") |> data.table::as.data.table()
#
#   harm <- TwoSampleMR::harmonise_data(g_formatted, r_formatted)


  return(list(gwas=gwas, qc=qc))
}


populate_rsid <- new_generic('populate_rsid', 'gwas', function(gwas, fill_rsid, missing_rsid, parallel_cores, verbose=TRUE) { S7_dispatch() })
method(populate_rsid, new_S3_class('data.table')) <- function(gwas, fill_rsid, missing_rsid, parallel_cores, verbose=TRUE) {

  if(verbose) message("[i] checking rsid coding")

  # parse existing RSID if there
  if("rsid" %in% names(gwas)) {

    gwas[, rsid := tolower(sub("^(rs[0-9]+).*", "\\1", rsid))]
    invalid <- which(!grepl("^(rs[0-9]+)$", gwas$rsid))

  } else {

    gwas[, rsid := NA_character_]
    invalid <- 1:nrow(gwas)

  }

  # some invalid RSIDs and want to populate from dbSNP
  if(length(invalid) > 0 && fill_rsid!=FALSE) {

    if(verbose) {
      message("\t[?] ", length(invalid), " RSIDs could not be parsed, attempting to fetch from dbSNP.")
    }

    # get the RSIDs
    rsid_dat <- chrpos_to_rsid(gwas[invalid, list(chr,bp,ea,oa)],
                               chr_col   = "chr",
                               pos_col   = "bp",
                               ea_col    = "ea",
                               nea_col   = "oa",
                               build     = fill_rsid,
                               flip      = "allow",
                               alt_rsids = FALSE,
                               verbose   = verbose,
                               parallel_cores=parallel_cores)
    # add back
    gwas[rsid_dat, rsid := i.RSID, on = c("chr", "bp", "ea", "oa")]
    rm(rsid_dat)
  }

  # what to do with missing RSIDs;
  switch(missing_rsid,
         'leave'                  = { gwas },
         'none'                   = { still_invalid <- which(!grepl("^(rs[0-9]+)$", gwas$rsid))
                                      if(verbose && length(still_invalid)>0) message("\t[?] ", length(still_invalid), " rsids could not be parsed, setting to NA.")
                                      gwas[invalid      , rsid := NA_character_                    ] },
         'fill_CHR:BP'            = { still_invalid <- which(!grepl("^(rs[0-9]+)$", gwas$rsid))
                                      if(verbose && length(still_invalid)>0) message("\t[?] ", length(still_invalid), " rsids could not be parsed, setting to chr:bp")
                                      gwas[still_invalid, rsid := paste0(chr,':',bp)               ] },
         'fill_CHR:BP_OA_EA'      = { still_invalid <- which(!grepl("^(rs[0-9]+)$", gwas$rsid))
                                      if(verbose && length(still_invalid)>0) message("\t[?] ", length(still_invalid), " rsids could not be parsed, setting to chr:bp_oa_ea")
                                      gwas[still_invalid, rsid := paste0(chr,':',bp,'_',oa,'_',ea) ] },
         'overwrite_CHR:BP'       = { if(verbose) message("\t[?] setting all rsids to chr:bp coding")
                                      gwas[             , rsid := paste0(chr,':',bp)               ] },
         'overwrite_CHR:BP:OA:EA' = { if(verbose) message("\t[?] setting all rsids to chr:bp_oa_ea coding")
                                      gwas[             , rsid := paste0(chr,':',bp,'_',oa,'_',ea) ] })

  # RSIDs at the front
  data.table::setcolorder(gwas, "rsid")

  # return
  return(gwas)
}


#' @title as.data.table
#' @param object GWAS object to covert to data.table
#' @param ... argument for data.table generic, ignored in this implementation
#' @export
as.data.table <- new_generic('as.data.table', 'object')
method(as.data.table, GWAS) <- function(object, ...) {

  data.table::data.table(
    rsid  = object@rsid,
    chr   = object@chr,
    bp    = object@bp,
    ea    = object@ea,
    oa    = object@oa,
    eaf   = object@eaf,
    beta  = object@beta,
    se    = object@se,
    p     = object@p,
    strand  = if(length(object@strand)==0)  { rep(NA_character_, length(object@rsid)) } else if(length(object@strand)==1)  { rep(object@strand, length(object@rsid)) }  else { object@strand },
    imputed = if(length(object@imputed)==0) { rep(NA, length(object@rsid)) }            else if(length(object@imputed)==1) { rep(object@imputed, length(object@rsid)) } else { object@imputed },
    info    = if(length(object@info)==0)    { rep(NA_real_, length(object@rsid)) }      else if(length(object@info)==1)    { rep(object@info, length(object@rsid)) }    else { object@info },
    n       = if(length(object@n)==0)       { rep(NA_integer_, length(object@rsid)) }   else if(length(object@n)==1)       { rep(object@n, length(object@rsid)) }       else { object@n },
    ncase   = if(length(object@ncase)==0)   { rep(NA_integer_, length(object@rsid)) }   else if(length(object@ncase)==1)   { rep(object@ncase, length(object@rsid)) }   else { object@ncase },
    trait   = if(length(object@trait)==0)   { rep(NA_character_, length(object@rsid)) } else if(length(object@trait)==1)   { rep(object@trait, length(object@rsid)) }   else { object@trait },
    id      = if(length(object@id)==0)      { rep(NA_character_, length(object@rsid)) } else if(length(object@id)==1)      { rep(object@id, length(object@rsid)) }      else { object@id }
  )

}


as.twosample.mr <- new_generic('as.twosample.mr', 'x')
method(as.twosample.mr, GWAS) <- function(x, type) { return( as.twosample.mr(list(x), type) ) }
method(as.twosample.mr, class_list) <- function(x, type) {

  # checks
  type <- match.arg(type, choices=c("exposure","outcome"))
  stopifnot("Exposure input list must all be of GWAS class" = all(sapply(x, function(x0) inherits(x0, "genepi.utils::GWAS"))))

  # ensure phenotype/trait col will be unique - (the hack here to avoid de-duplication of multiple exposures is
  # that format de-duplicates by phenotype_col)
  stopifnot("Duplicate exposure ids found" = !any(duplicated(sapply(x, function(x0) x0@id))))
  if(any(duplicated(sapply(x, function(x0) x0@id)))) {
    x <- lapply(x, function(x0) x0@trait <- paste0(x0@trait, " - ", x0@id))
  }

  # exposures to large dt with unique `id`s
  dat <- lapply(x, function(x0) as.data.table(x0)) |> data.table::rbindlist()

  # recode as sorted alleles (to allow matching of unharmonised mutliallelic)
  dat[, rsid := paste0(c(rsid, sort(c(oa,ea))), collapse="_"), by=1:nrow(dat)]

  # 2SMR format
  formatted <- TwoSampleMR::format_data(
    dat               = dat |> as.data.frame(),
    type              = type,
    snp_col           = "rsid",
    chr_col           = "chr",
    pos_col           = "bp",
    effect_allele_col = "ea",
    other_allele_col  = "oa",
    eaf_col           = "eaf",
    beta_col          = "beta",
    se_col            = "se",
    pval_col          = "p",
    samplesize_col    = "n",
    ncase_col         = "ncase",
    id_col            = "id",
    phenotype_col     = "trait") |> data.table::as.data.table()

  # 2SMR converts to all lower-case - put the alleles back
  formatted[, SNP := sub("^(rs[0-9]+|[0-9X]+:[0-9]+)(.*)","\\1\\U\\2", SNP, perl=TRUE)]

  # return
  return(formatted)

}

