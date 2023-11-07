# silence R CMD checks for data.table columns
BETA = BP = CHR = EA = OA = OR = OR_LB = OR_SE = OR_UB = P = SE = SNP = NULL

#' @title Standardise GWAS format
#' @description
#' Standardises a GWAS with appropriate column names and data checks.
#' @param gwas a string file path, or data.frame like object
#' @param input_format a string or list(); if string, must be a valid format - these can be obtained with `input_formats(ColumnMapping())`
#' @param drop a logical, whether to drop additional columns not in the input_format mapping
#' @param fill a logical, whether to add (NAs) missing columns present in the mapping but not present in the GWAS
#' @param build a string, one of GRCh37 or GRCh38
#' @param populate_rsid either FALSE or a valid argument for chrpos_to_rsid `build`, e.g. b37_dbsnp156
#' @return a standardised data.table
#' @export
#'
standardise_gwas <- function(gwas, input_format="default", drop=FALSE, fill=FALSE, build="GRCh37", populate_rsid=FALSE) {
#
#   load_all()
#   gwas="/Users/xx20081/Documents/local_data/hermes_progression/paradigm_hf/raw/paradigm.females.comp2.gz"
#   gwas <- data.table::fread(gwas, nThread=parallel::detectCores())
#
#   # overwrite "D"/"I" coding if present
#   if(any(grepl("^(D|I)$", gwas[[mapping[["EA"]]]]))) {
#     gwas[, mapping[["EA"]] := ORI_EFFECT_ALLELE]
#     gwas[, mapping[["OA"]] := ORI_OTHER_ALLELE]
#   }
#
#   drop=T
#   fill=T
#   build="GRCh37"
#   populate_rsid=FALSE
#   input_format=list(
#     "SNP"=       "MARKER",
#       "CHR"=       "CHR",
#       "BP"=        "POS",
#       "EA"=        "EFFECT_ALLELE",
#       "OA"=        "OTHER_ALLELE",
#       "BETA"=      "BETA",
#       "SE"=        "SE",
#       "P"=         "P",
#       "EAF"=       "EAF",
#       "STRAND"=    "STRAND",
#       "N"=         "N_SAMPLE",
#       "N_CASE"=    "N_EVENTS",
#       "IMPUTED"=   "Imputed",
#       "INFO"=      "QUAL_SCORE"
#   )

  # validate inputs
  build <- match.arg(build, c("GRCh37", "GRCh38"))
  if(populate_rsid!=FALSE) {
    populate_rsid <- match.arg(populate_rsid, names(genepi.utils::which_dbsnp_builds()))
  }

  # get the mapping
  mapping <- get_gwas_mapping(input_format)

  # apply the standardisation procedures
  gwas <- get_input_gwas(gwas) |>
    change_column_names(mapping, drop, fill) |>
    standardise_columns() |>
    filter_incomplete_rows() |>
    filter_invalid_values() |>
    standardise_alleles(build) |>
    health_check() |>
    populate_rsid(populate_rsid)

  # return the standardised gwas
  return(gwas)
}


#' @title Change column names
#' @description
#' Change_column_names: private function that takes a named list of column names
#' and changes the supplied data frame's column names accordingly.
#' @param: gwas of gwas to standardise column names
#' @param: columns named list for
#' @param: opposite_mapping logical flag on if we are mapping from key to value or vice verca
#' @return a data.table
#' @noRd
#'
change_column_names <- function(gwas, columns = list(), drop=FALSE, fill=FALSE, opposite_mapping = FALSE) {

  if(fill) {

    gwas[, unlist(unname(columns))[!unlist(unname(columns)) %in% names(gwas)] := NA]

  } else {

    columns <- columns[unlist(unname(columns)) %in% names(gwas)]

  }

  if (!opposite_mapping) {

    data.table::setnames(gwas, unlist(unname(columns)), names(columns))

  } else {

    data.table::setnames(gwas, names(columns), unlist(unname(columns)))

  }

  if(drop) {

    gwas[, names(gwas)[!names(gwas) %in% names(columns)] := NULL]

  }



  return(gwas)
}


#' @title Standardise column data
#' @description
#' Adjust columns to meeting minimum column specification, including conversion of odds ratio to
#' beta and standard error
#' @param gwas a data.table
#' @return a standardised gwas
#' @noRd
#'
standardise_columns <- function(gwas) {

  # ensure CHR and BP columns present
  if (!all(c("CHR", "BP") %in% names(gwas))) {

    valid_pattern <- grepl("(?i)(\\d+|X|Y)[:_]\\d+", gwas[["SNP"]])

    if(all(valid_pattern)) {

      gwas[, CHR := sub("(?i)(\\d+|X|Y)[:_]\\d+(?:.*)", "\\1", SNP)]
      gwas[, BP  := sub("(?i)(?:\\d+|X|Y)[:_](\\d+)(?:.*)", "\\1", SNP)]

    } else {

      example_failure <- gwas[[which(!valid_pattern)[[1]], "SNP"]]
      stop(paste0("problem standardising CHR and POS columns. No CHR BP columns provided and cannot parse SNP column with regex \\d[:_]\\d and enteries such as: ", example_failure))

    }
  }

  # convert odds ratio to betas and standard error, if BETA, SE not provided
  if(!all(c("BETA", "SE") %in% names(gwas))) {

    if (all(c("OR", "OR_LB", "OR_UB") %in% names(gwas))) {

      gwas <- convert_or_to_beta(gwas)

    } else if (all(c("OR", "OR_SE") %in% names(gwas))) {

      gwas <- convert_or_to_beta(gwas)

    } else {

      stop("problem standardising/creating BETA and SE columns. No BETA SE columns provided and no OR OR_LB OR_UB or OR OR_SE columns provided to compute them.")

    }
  }

  # convert Z-score to Pvalue if P not provided but Z score is
  if(!"P" %in% names(gwas)) {

    if("Z" %in% names(gwas)) {

      gwas <- convert_z_to_p(gwas)

    } else {

      stop("problem standardising/creating P value column. No P value column found and no Z column provided to compute it.")

    }

  }

  # correct types
  gwas[, CHR := as.character(CHR)]
  gwas[, BP  := as.integer(BP)]
  gwas[, BETA:= as.numeric(BETA)]
  gwas[, SE  := as.numeric(SE)]
  gwas[, P   := as.numeric(P)]
  gwas[P==0, P := .Machine$double.xmin]

  return(gwas)
}


#' @title Convert Z score to p-value
#' @description
#' Given a Z-score, calculates the p-value.
#' @param gwas a data.table with the following columns: Z
#' @return a data.table with new column P
#' @noRd
#'
convert_z_to_p <- function(gwas) {

  gwas[, P := 2 * pnorm(-abs(Z))]

  return(gwas)
}


#' @title Convert odds ratio to beta
#' @description
#' Given an OR and lower and upper bounds, calculates the BETA, and SE.
#' Based on this answer: https://stats.stackexchange.com/a/327684
#' @param gwas a data.table with the following columns: OR, LB (lower bound), UB (upper bound)
#' @return a data.table with new columns BETA and SE
#' @noRd
#'
convert_or_to_beta <- function(gwas) {

  gwas[, BETA := log(OR)]

  if ("OR_SE" %in% names(gwas)) {

    gwas[, OR_UB := OR + OR_SE * 1.96]
    gwas[, OR_UB := OR - OR_SE * 1.96]

  }

  gwas[, SE := (OR_UB - OR_LB) / (2 * 1.95996)]

  return(gwas)
}


#' @title Filter incomplete rows
#' @description
#' Remove rows where values are missing or invalid.
#' @param gwas a data.table
#' @return a filtered data.table
#' @noRd
#'
filter_incomplete_rows <- function(gwas) {

  check_na <- c("CHR","BP","EAF","BETA","SE","P")

  missing <- gwas[, rowSums(is.na(.SD)) > 0, .SDcols=check_na]

  if (any(missing)) {
    warning(paste0("Filtering out ", sum(missing), " rows due to NAs one or more of columns: ", paste0(check_na, collapse=", ")))
  }

  return(gwas[!missing, ])
}


#' @title Filter invalid data values
#' @description
#' Remove rows where values do not make sense.
#' @param gwas a data.table
#' @return a filtered data.table
#' @noRd
#'
filter_invalid_values <- function(gwas) {

  invalid_flags <- list(
    "BETA" = is.infinite(gwas[["BETA"]]),
    "P"    = is.infinite(gwas[["P"]]),
    "SE"   = is.infinite(gwas[["SE"]])
  )

  report_cols <- c("CHR","BP","BETA","SE","P")

  any_invalid <- rep(FALSE, nrow(gwas))

  for(i in seq_along(invalid_flags)) {

    col  <- names(invalid_flags)[i]

    invalid <- invalid_flags[[i]]

    if (sum(invalid, na.rm=T) > 0) {
      warning(paste0("Filtering out ", sum(invalid, na.rm=T), " rows due to invalid data in the ", col, " column."))
    }

    any_invalid <- any_invalid | invalid
  }

  return(gwas[!any_invalid, ])
}


#' @title Standardise alleles
#' @description
#' Standardise alleles to upper
#' @param gwas a data.table
#' @return a data.table
#' @importFrom future plan multisession
#' @importFrom progressr with_progress
#' @noRd
#'
standardise_alleles <- function(gwas, build) {

  # make uppercase
  gwas[, EA := toupper(EA)]
  gwas[, OA := toupper(OA)]

  # define valid allele formats (could include 'D' and 'I' here, chosen not to for now)
  valid_alleles <- ((grepl("^[AGCT]+$", gwas[["EA"]]) | is.na(gwas[["EA"]])) & grepl("^[AGCT]+$", gwas[["OA"]])) |
                   ((grepl("^[AGCT]+$", gwas[["OA"]]) | is.na(gwas[["OA"]])) & grepl("^[AGCT]+$", gwas[["EA"]]))

  # report invalid
  if(any(!valid_alleles)) {

    example_invalid <- paste0(gwas[[which(!valid_alleles)[[1]], "EA"]], "/", gwas[[which(!valid_alleles)[[1]], "EA"]])
    warning(paste0("filtering out ", sum(!valid_alleles), " invalid alleles, for example: `", example_invalid, "`"))

  }

  # filter out invalid
  gwas <- gwas[valid_alleles, ]

  # generate the SNP string; plink2 format e.g. "1:100000[b37]REF,ALT"
  if(build == "GRCh37") {

    gwas[, SNP := paste0(CHR,":",BP,"[b37]",OA,",",EA)]

  } else if (build == "GRCh38") {

    gwas[, SNP := paste0(CHR,":",BP,"[b38]",OA,",",EA)]

  }

  # nice column ordering
  cols <- c("SNP","CHR","BP","EA","OA","EAF","BETA","SE","P")
  data.table::setcolorder(gwas, c(cols, names(gwas)[!names(gwas) %in% cols]))

  return(gwas)
}


#' @title Populate RSID column
#' @description
#' Populates an RSID field from CHR, POS, EA and OA information.
#' @param gwas a data.table
#' @param populate_rsid either FALSE (doesn't run), or a valid genome and dbSNP `build` argument for chrpos_to_rsid()
#' @return a data.table
#' @noRd
#'
populate_rsid <- function(gwas, populate_rsid=FALSE) {

  if (populate_rsid == FALSE) return(gwas)

  if("RSID" %in% names(gwas)) {
    warning("GWAS already has an RSID field, this will not be overwritten")
    return(gwas)
  }

  # get the RSIDs
  future::plan(future::multisession, workers=parallel::detectCores())
  progressr::with_progress({
    gwas <- chrpos_to_rsid(gwas,
                           chr_col   = "CHR",
                           pos_col   = "BP",
                           ea_col    = "EA",
                           nea_col   = "OA",
                           build     = populate_rsid,
                           flip      = "report",
                           alt_rsids = FALSE,
                           verbose   = TRUE)
  })

  # RSIDs at the front
  cols <- c("RSID", "rsid_flip_match")
  data.table::setcolorder(gwas, c(cols, names(gwas)[!names(gwas) %in% cols]))

  return(gwas)
}


#' @title GWAS data health check
#' @description
#' Perform several health checks on the data
#' @param gwas a data.table
#' @return a data.table
#' @noRd
#'
health_check <- function(gwas) {

  invalid_pvals <- gwas[["P"]] <= 0 | gwas[["P"]] > 1 | is.na(gwas[["P"]])

  if (sum(invalid_pvals, na.rm=TRUE) > 0 | any(is.na(invalid_pvals))) {
    cols <- c("SNP","CHR","BP","EA","OA","EAF","BETA","SE","P")
    warning(paste0("GWAS has ", sum(invalid_pvals), " P value(s) outside accepted range, for example:"))
    print(gwas[which(invalid_pvals)[1], cols, with=FALSE])
    stop("invalid GWAS P value(s)")
  }

  invalid_eaf <- gwas[["EAF"]] < 0 | gwas[["EAF"]] > 1 | is.na(gwas[["EAF"]])

  if (sum(invalid_eaf, na.rm=T) > 0 | any(is.na(invalid_eaf))) {
    cols <- c("SNP","CHR","BP","EA","OA","EAF","BETA","SE","P")
    warning(paste0("GWAS has ", sum(invalid_eaf), " EAF value(s) outside accepted range, for example:"))
    print(gwas[which(invalid_eaf)[1], cols, with=FALSE])
    stop("invalid GWAS EAF value(s)")
  }

  return(gwas)
}


#' @title Harmonise multiple GWAS
#' @description
#' This takes a list of gwases, gets the SNPs in common across all datasets
#' and arranges them to be in the same order.
#' @param ... ellipses of gwas data.tables
#' @return list of harmonised gwases
#' @export
#'
harmonise_gwases <- function(...) {

  # to a list
  gwases <- list(...)

  # intersection of SNP ids
  snpids <- Reduce(intersect, lapply(gwases, function(gwas) gwas[["SNP"]]))

  # report
  message(paste("Number of shared SNPs after harmonisation:", length(snpids)))

  # apply filtering and order
  gwases <- lapply(gwases, function(gwas) {
    gwas <- gwas[SNP %in% snpids, ]
    data.table::setkey(gwas, "SNP")
  })

  return(gwases)
}


#' @title Get / format the input GWAS
#' @param gwas a file path, or data.frame like object
#' @return a data.table
#' @noRd
#'
get_input_gwas <- function(gwas) {

  if(is.character(gwas)) {

    stopifnot("`gwas` must be a valid file path is supplied as a character" = file.exists(gwas))
    gwas <- data.table::fread(gwas, nThread=parallel::detectCores())

  } else if(inherits(gwas, "data.frame")) {

    gwas <- data.table::as.data.table(gwas)

  } else {

    stop("`gwas` must be a valid file path or data.frame like object")

  }

  return(gwas)
}


#' @title Get the GWAS mapping
#' @param input_format a file path, or data.frame like object
#' @return a data.table
#' @noRd
#'
get_gwas_mapping <- function(input_format) {

  if(is.character(input_format)) {

    mapping <- get_map(ColumnMapping(), input_format)

  } else if(is.list(input_format)) {

    mapping <- input_format

  } else {

    valid_formats <- input_formats(ColumnMapping())
    stop(paste0("`input_format` must be a defined format [", paste0(valid_formats, collapse=", "), "], or a custom list() with at least named elements c(SNP=..., BETA=..., SE=..., EA=..., OA=..., EAF=..., P=..."))

  }

  return(mapping)
}


#' @title Format a GWAS file
#' @description
#' This formats a GWAS file column names and writes out to `output_file`
#' @param file_gwas a string file path, or data.frame like object
#' @param output_file a string, a file path
#' @param output_format a string, a valid format - these can be obtained with `input_formats(ColumnMapping())`
#' @return a file written to output_file
#' @export
#'
format_gwas_output <- function(file_gwas, output_file, output_format="default") {

  mapping <- get_gwas_mapping(output_format)

  gwas <- get_input_gwas(gwas) |>
    change_column_names(mapping, drop=FALSE, opposite_mapping=TRUE)

  data.table::fwrite(gwas, output_file, sep="\t", nThread=parallel::detectCores())
}




# TODO:
# populate_snp_from_rsid <- function(gwas) {
#   start <- Sys.time()
#   marker_to_rsid_file <- paste0(thousand_genomes_dir, "marker_to_rsid_full.tsv.gz")
#   marker_to_rsid <- vroom::vroom(marker_to_rsid_file, col_select=c("HG37", "RSID"))
#   print(paste("loaded file: ", Sys.time()-start))
#
#   matching <- match(gwas$RSID, marker_to_rsid$RSID)
#   gwas$CHRBP<- marker_to_rsid$HG37[matching]
#   gwas <- tidyr::separate(data = gwas, col = "CHRBP", into = c("CHR", "BP"), sep = ":")
#   print(paste("mapped and returned: ", Sys.time()-start))
#
#   return(gwas)
# }




