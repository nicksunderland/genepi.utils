# silence R CMD checks for data.table columns
BETA = BP = CHR = OR = OR_LB = OR_SE = OR_UB = P = SE = SNP = NULL

#' @title standardise_gwas
#' @description
#' Standardises a GWAS with appropriate column names and data checks.
#' @param gwas a string file path, or data.frame like object
#' @param input_format a string or list()
#' @param populate_rsid a logical
#' @return a standardised data.table
#' @export
#'
standardise_gwas <- function(gwas, input_format="default", populate_rsid=FALSE) {

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
    standardise_columns()
  #|>
    # standardise_alleles() |>
    # health_check() |>
    # populate_rsid_from_1000_genomes(populate_rsid)

  # return the standardised gwas
  return(gwas)
}


#' @title change_column_names
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


#' @title filter_incomplete_rows
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


#' @title standardise_columns
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


#' @title convert_or_to_beta
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

#'
#'
#' #' @title standardise_alleles
#' #' @description
#' #' Standardise alleles to upper
#' #' @param gwas a data.table
#' #' @return a data.table
#' #' @noRd
#' #'
#' standardise_alleles <- function(gwas) {
#'
#'
#'   # TODO - ask AE about flips
#'
#'   gwas$EA <- toupper(gwas$EA)
#'   gwas$OA <- toupper(gwas$OA)
#'
#'   to_flip <- gwas$EA > gwas$OA
#'   gwas$EAF[to_flip] <- 1 - gwas$EAF[to_flip]
#'   gwas$BETA[to_flip] <- -1 * gwas$BETA[to_flip]
#'
#'   temp <- gwas$OA[to_flip]
#'   gwas$OA[to_flip] <- gwas$EA[to_flip]
#'   gwas$EA[to_flip] <- temp
#'
#'   gwas$SNP <- toupper(paste0(gwas$CHR, ":", format(gwas$BP, scientific = F, trim = T), "_", gwas$EA, "_", gwas$OA))
#'   gwas %>% dplyr::select(SNP, CHR, BP, EA, OA, EAF, BETA, SE, P, everything())
#'
#'   return(gwas)
#' }
#'
#'
#'
#'
#'
#'
#' format_gwas_output <- function(file_gwas, output_file, output_format="default") {
#'   gwas <- vroom::vroom(file_gwas) %>%
#'     change_column_names(column_map[[output_format]], opposite_mapping = T)
#'
#'   vroom:vroom_write(gwas, output_file, delim="\t")
#' }
#'
#'
#'
#'
#'
#'
#'
#'
#' health_check <- function(gwas) {
#'   if (nrow(gwas[gwas$P <= 0 | gwas$P > 1, ]) > 0) {
#'     stop("GWAS has some P values outside accepted range")
#'   }
#'   if (nrow(gwas[gwas$EAF < 0 | gwas$EAF > 1, ]) > 0) {
#'     stop("GWAS has some EAF values outside accepted range")
#'   }
#'   #if ("OR" %in% colnames(gwas) & nrow(gwas[gwas$OR < 0, ]) > 0) {
#'   #  stop("GWAS has some OR values outside accepted range")
#'   #}
#'   return(gwas)
#' }
#'
#'
#'
#'
#'
#'
#' #' harmonise_gwases: takes a list of gwases, get the SNPs in common
#' #' across all datasets arranged to be in the same order
#' #'
#' #' @param: elipses of gwases
#' #' @return: list of harmonised gwases
#' #'
#' harmonise_gwases <- function(...) {
#'   gwases <- list(...)
#'
#'   snpids <- Reduce(intersect, lapply(gwases, function(gwas) gwas$SNP))
#'   print(paste("Number of shared SNPs after harmonisation:", length(snpids)))
#'
#'   gwases <- lapply(gwases, function(gwas) {
#'     gwas %>%
#'       dplyr::filter(SNP %in% snpids) %>%
#'       dplyr::arrange(SNP)
#'   })
#'
#'   return(gwases)
#' }
#'
#' populate_snp_from_rsid <- function(gwas) {
#'   start <- Sys.time()
#'   marker_to_rsid_file <- paste0(thousand_genomes_dir, "marker_to_rsid_full.tsv.gz")
#'   marker_to_rsid <- vroom::vroom(marker_to_rsid_file, col_select=c("HG37", "RSID"))
#'   print(paste("loaded file: ", Sys.time()-start))
#'
#'   matching <- match(gwas$RSID, marker_to_rsid$RSID)
#'   gwas$CHRBP<- marker_to_rsid$HG37[matching]
#'   gwas <- tidyr::separate(data = gwas, col = "CHRBP", into = c("CHR", "BP"), sep = ":")
#'   print(paste("mapped and returned: ", Sys.time()-start))
#'
#'   return(gwas)
#' }
#'
#' populate_rsid_from_1000_genomes <- function(gwas, populate_rsid=F) {
#'   if (populate_rsid == F) return(gwas)
#'
#'   #if (populate_rsid == "full") {
#'   #  marker_to_rsid_file <- paste0(thousand_genomes_dir, "marker_to_rsid_full.tsv.gz")
#'   #}
#'   #else {
#'   #  marker_to_rsid_file <- paste0(thousand_genomes_dir, "marker_to_rsid.tsv.gz")
#'   #}
#'
#'   if(column_map$default$RSID %in% colnames(gwas)) {
#'     print("GWAS already has an RSID field, will not overwrite")
#'     return(gwas)
#'   }
#'   print("populating RSID...")
#'   marker_to_rsid_file <- paste0(genome_data_dir, "marker_to_rsid.tsv.gz")
#'   chrpos_to_rsid <- vroom::vroom(marker_to_rsid_file, col_select=c("HG37", "RSID"))
#'   gwas$RSID <- chrpos_to_rsid$RSID[match(gwas$SNP, chrpos_to_rsid$HG37)]
#'
#'   return(gwas)
#' }
