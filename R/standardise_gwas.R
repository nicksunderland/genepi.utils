# silence R CMD checks for data.table columns
BETA = BP = CHR = EA = OA = OR = OR_LB = OR_SE = OR_UB = P = SE = SNP = NULL

#' @title Standardise GWAS format
#' @description
#' Standardises a GWAS with appropriate column names and data checks.
#' @param gwas a string file path, or data.frame like object
#' @param input_format a string or list()
#' @param build a string, one of GRCh37 or GRCh38
#' @param populate_rsid either FALSE or a valid argument for chrpos_to_rsid `build`, e.g. b37_dbsnp156
#' @return a standardised data.table
#' @export
#'
standardise_gwas <- function(gwas, input_format="default", build="GRCh37", populate_rsid=FALSE) {

  # validate inputs
  build <- match.arg(build, c("GRCh37", "GRCh38"))
  if(populate_rsid!=FALSE) {
    populate_rsid <- match.arg(populate_rsid, names(genepi.utils::which_dbsnp_builds()))
  }

  # define the column mapping
  if(is.character(input_format)) {

    mapping <- get_map(ColumnMapping(), input_format)

  } else if(is.list(input_format)) {

    mapping <- input_format

  } else {

    valid_formats <- input_formats(ColumnMapping())
    stop(paste0("`input_format` must be a defined format [", paste0(valid_formats, collapse=", "), "], or a custom list() with at least named elements c(SNP=..., BETA=..., SE=..., EA=..., OA=..., EAF=..., P=..."))

  }

  # read the gwas if a file path, or ensure is a data.frame like object
  if(is.character(gwas)) {

    stopifnot("`gwas` must be a valid file path is supplied as a character" = file.exists(gwas))
    gwas <- data.table::fread(gwas, nThread=parallel::detectCores())

  } else if(inherits(gwas, "data.frame")) {

    gwas <- data.table::as.data.table(gwas)

  } else {

    stop("`gwas` must be a valid file path or data.frame like object")

  }

  # apply the standardisation procedures
  gwas <- gwas |>
    change_column_names(mapping) |>
    filter_incomplete_rows() |>
    standardise_columns() |>
    standardise_alleles() |>
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
change_column_names <- function(gwas, columns = list(), opposite_mapping = FALSE) {

  if (!opposite_mapping) {

    data.table::setnames(gwas, unname(columns), names(columns))

  } else {

    data.table::setnames(gwas, names(columns), unname(columns))

  }

  return(gwas)
}


#' @title Filter incomplete rows
#' @description
#' Remove rows where "CHR","BP","EA","OA","EAF" are all missing.
#' @param gwas a data.table
#' @return a filtered data.table
#' @noRd
#'
filter_incomplete_rows <- function(gwas) {

  check_na <- c("CHR","BP","EA","OA","EAF")

  empty_row_idx <- gwas[, rowSums(is.na(.SD))==length(check_na), .SDcols=check_na]

  if (length(empty_row_idx) > 0) {
    warning(paste0("Filtering out ", length(empty_row_idx), " rows due to NAs"))
  }

  return(gwas[!empty_row_idx, ])
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

    valid_pattern <- grepl("\\d[:_]\\d", gwas[["SNP"]])

    if(all(valid_pattern)) {

      gwas[, CHR := sub("(\\d)[:_]\\d(?:.*)", "\\1", SNP)]
      gwas[, BP  := sub("\\d[:_](\\d)(?:.*)", "\\1", SNP)]

    } else {

      example_failure <- gwas[[which(!valid_pattern)[[1]], "SNP"]]
      stop(paste0("problem standardising CHR and POS columns. No CHR BP columns provided and cannot parse SNP column with regex \\d[:_]\\d and enteries such as: ", example_failure))

    }
  }

  # convert odds ratio to betas and standard error
  if (all(c("OR", "OR_LB", "OR_UB") %in% names(gwas)) & !all(c("BETA", "SE") %in% names(gwas))) {

    gwas <- convert_or_to_beta(gwas)

  } else if (all(c("OR", "OR_SE") %in% names(gwas)) & !all(c("BETA", "SE") %in% names(gwas))) {

    gwas <- convert_or_to_beta(gwas)

  } else {

    stop("problem standardising BETA and SE columns. No BETA SE columns provided and no OR OR_LB OR_UB or OR OR_SE columns provided to compute them.")

  }

  # correct types
  gwas[, BP := as.integer(BP)]
  gwas[, P  := as.numeric(BP)]
  gwas[P==0, P := .Machine$double.xmin]

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


#' @title Standardise alleles
#' @description
#' Standardise alleles to upper
#' @param gwas a data.table
#' @return a data.table
#' @noRd
#'
standardise_alleles <- function(gwas, build) {

  # make uppercase
  gwas[, EA := toupper(EA)]
  gwas[, OA := toupper(OA)]

  # define valid allele formats (could include 'D' and 'I' here, chosen not to for now)
  valid_alleles <- ((grepl("^[AGCT]+$", gwas[["EA"]]) | is.na(EA)) & grepl("^[AGCT]+$", gwas[["OA"]])) |
                   ((grepl("^[AGCT]+$", gwas[["OA"]]) | is.na(OA)) & grepl("^[AGCT]+$", gwas[["EA"]]))

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
  gwas <- chrpos_to_rsid(gwas,
                         chr_col   = "CHR",
                         pos_col   = "BP",
                         ea_col    = "EA",
                         nea_col   = "OA",
                         build     = populate_rsid,
                         flip      = "report",
                         alt_rsids = FALSE,
                         verbose   = TRUE)

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

  invalid_pvals <- gwas[["P"]] <= 0 | gwas[["P"]] > 1
  if (sum(invalid_pvals) > 0) {
    cols <- c("SNP","CHR","BP","EA","OA","EAF","BETA","SE","P")
    warning(paste0("GWAS has ", sum(invalid_pvals), " P value(s) outside accepted range, for example:"))
    print(gwas[which(invalid_pvals)[1], cols, with=FALSE])
    stop("invalid GWAS P value(s)")
  }

  invalid_eaf <- gwas[["EAF"]] < 0 | gwas[["EAF"]] > 1
  if (sum(invalid_eaf) > 0) {
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

# TODO:
# format_gwas_output <- function(file_gwas, output_file, output_format="default") {
#   gwas <- vroom::vroom(file_gwas) %>%
#     change_column_names(column_map[[output_format]], opposite_mapping = T)
#
#   vroom:vroom_write(gwas, output_file, delim="\t")
# }

