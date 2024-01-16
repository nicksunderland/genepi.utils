# silence R CMD checks for data.table columns
BETA = BP = CHR = EA = OA = OR = OR_LB = OR_SE = OR_UB = P = Z = SE = SNP = NULL

#' @title Standardise GWAS format
#' @description
#' Standardises a GWAS with appropriate column names and data checks.
#' @param gwas a string file path, or data.frame like object
#' @param input_format a string or list(); if string, must be a valid format - these can be obtained with `input_formats(ColumnMapping())`
#' @param filters a named list of conditions (as strings) to be evaluated as data.table row filters (quoting column names not necessary).
#' These must be positive assertions, i.e. resolve to TRUE to keep the row. If filters==NULL, no filtering occurs. Example:
#' `filters=list('abs(BETA)<10', 'EAF>0.01 & EAF<0.99')` ...will filter for rows with BETA value between -10 and +10 and EAF between
#' 1% and 99%.
#' @param drop a logical, whether to drop additional columns not in the input_format mapping
#' @param fill a logical, whether to add (NAs) missing columns present in the mapping but not present in the GWAS
#' @param build a string, one of GRCh37 or GRCh38
#' @param populate_rsid either FALSE or a valid argument for chrpos_to_rsid `build`, e.g. b37_dbsnp156
#' @param missing_rsid one of c("fill_CHR:BP","fill_CHR:BP_OA_EA","overwrite_CHR:BP","overwrite_CHR:BP:OA:EA","none")
#' @return a standardised data.table
#' @export
#'
standardise_gwas <- function(gwas,
                             input_format="default",
                             drop=FALSE,
                             fill=FALSE,
                             build="GRCh37",
                             populate_rsid=FALSE,
                             missing_rsid="fill_CHR:BP",
                             filters = list(
                               BETA_invalid = "!is.infinite(BETA) & abs(BETA) < 20",
                               P_invalid    = "!is.infinite(P)",
                               SE_invalid   = "!is.infinite(SE)",
                               CHR_missing  = "!is.na(CHR)",
                               BP_missing   = "!is.na(BP)",
                               BETA_missing = "!is.na(BETA)",
                               SE_missing   = "!is.na(SE)",
                               P_missing    = "!is.na(P)",
                               EAF_missing  = "!is.na(EAF)"
                             )) {

  message("Standardising GWAS ...")

  # validate inputs
  build <- match.arg(build, c("GRCh37", "GRCh38"))
  missing_rsid <- match.arg(missing_rsid, c("fill_CHR:BP","fill_CHR:BP_OA_EA","overwrite_CHR:BP","overwrite_CHR:BP:OA:EA","none"))
  if(populate_rsid!=FALSE) {
    populate_rsid <- match.arg(populate_rsid, names(genepi.utils::which_dbsnp_builds()))
  }

  # get the mapping
  mapping <- get_gwas_mapping(input_format)

  # get the gwas from file or object
  gwas <- import_table(gwas)

  # set attribute for starting number of rows
  data.table::setattr(gwas, "input rows", nrow(gwas))

  # apply the standardisation procedures
  gwas <- change_column_names(gwas, mapping, drop, fill)
  gwas <- standardise_columns(gwas)
  gwas <- filter_invalid_values(gwas, filters)
  gwas <- standardise_alleles(gwas, build)
  gwas <- health_check(gwas)
  gwas <- populate_rsid(gwas, populate_rsid, missing_rsid)

  # set attribute for final number of valid rows
  data.table::setattr(gwas, "valid rows", nrow(gwas))

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

    missing_cols <- unlist(columns)[!unlist(unname(columns)) %in% names(gwas)]

    if(length(missing_cols) > 0) {

      message(paste0("[+] adding missing column(s) [", paste0(names(missing_cols),"=",missing_cols, collapse=", "), "] specified in mapping but not found in gwas data  (see `fill` option)"))
      gwas[, unname(missing_cols) := NA]
      data.table::setattr(gwas, "added columns", missing_cols)

    }

  } else {

    columns <- columns[unlist(unname(columns)) %in% names(gwas)]

  }

  if (!opposite_mapping) {

    data.table::setnames(gwas, unlist(unname(columns)), names(columns))

  } else {

    data.table::setnames(gwas, names(columns), unlist(unname(columns)))

  }

  if(drop) {

    unwanted <- names(gwas)[!names(gwas) %in% names(columns)]

    if(length(unwanted) > 0) {

      message(paste0("[-] dropping unwanted column(s) [", paste0(unwanted, collapse=", "), "] found in gwas data but not specified in mapping (see `drop` option)"))
      gwas[, (unwanted) := NULL]
      data.table::setattr(gwas, "dropped columns", unwanted)

    }

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
#' @importFrom stats pnorm
#' @noRd
#'
convert_z_to_p <- function(gwas) {

  message("[+] converting Z to P")

  gwas[, P := 2 * stats::pnorm(-abs(Z))]

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

  message("[+] converting OR to BETA")

  gwas[, BETA := log(OR)]

  if ("OR_SE" %in% names(gwas)) {

    gwas[, OR_UB := OR + OR_SE * 1.96]
    gwas[, OR_UB := OR - OR_SE * 1.96]

  }

  gwas[, SE := (OR_UB - OR_LB) / (2 * 1.95996)]

  return(gwas)
}


#' #' @title Filter incomplete rows
#' #' @description
#' #' Remove rows where values are missing or invalid.
#' #' @param gwas a data.table
#' #' @return a filtered data.table
#' #' @noRd
#' #'
#' filter_incomplete_rows <- function(gwas) {
#'
#'   check_na <- c("CHR","BP","EAF","BETA","SE","P")
#'
#'   missing <- gwas[, rowSums(is.na(.SD)) > 0, .SDcols=check_na]
#'
#'   if (any(missing)) {
#'     warning(paste0("Filtering out ", sum(missing), " rows due to NAs one or more of columns: ", paste0(check_na, collapse=", ")))
#'   }
#'
#'   return(gwas[!missing, ])
#' }


#' @title Filter invalid data values
#' @description
#' Remove rows where values do not make sense.
#' @param gwas a data.table
#' @return a filtered data.table
#' @noRd
#'
filter_invalid_values <- function(gwas, filters=NULL) {

  # apply some filters if provided
  if(!is.null(filters)) {

    # order the filters so that removal of NAs is first, in order to minimise
    # problems later calling conditions that might throw errors when called
    # with NA values
    filters <- c(grep("is\\.na", filters, value=TRUE), setdiff(filters, grep("is\\.na", filters, value=TRUE)))

    # assess each filter
    for(i in seq_along(filters)) {

      # the name of the filter, or assign a number
      filter_names <- names(filters)
      if(is.null(filter_names)) {
        filter_name <- i
      } else {
        filter_name <- filter_names[i]
      }

      # parse the string / filter expression
      filter_expr <- parse(text = filters[[i]])

      # apply filter / create row index logical mask
      filter_true <- gwas[, eval(filter_expr)]

      # if the filter has generated NAs, report which filter errored and stop
      if(any(is.na(filter_true))) {
        na_idx <- which(is.na(filter_true))
        msg <- paste0("filter ", filter_name, " [", filter_expr, "] is invalid as it produced ", length(na_idx), " NAs. Please use an 'is.na()' filter, adjust your current filter, and/or adjust your input file. Example error rows:")
        warning(msg)
        print(head(gwas[na_idx, ]))
        stop("Filter error")
      }

      # add the filter to attributes and report number of rows filtered
      data.table::setattr(gwas, paste0(filter_name, " [", filter_expr, "] rows removed"), sum(!filter_true))

      # apply the filter to the data.table only if any() fail the filter;
      # otherwise leave the table alone (i.e. don't copy it in memory / reassign with `<-`)
      if(any(!filter_true)) {

        gwas <- gwas[filter_true, ]
        message(paste0("[-] ", sum(!filter_true), " rows removed with filter '", filter_name, "' [", filter_expr, "]" ))

      }

    }

  }

  # return the GWAS
  return(gwas)
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

  # define valid allele formats - NA allowed as deletion
  valid_alleles <- (((grepl("^([AGCT]+|[D])$", gwas[["EA"]]) | is.na(gwas[["EA"]])) & grepl("^([AGCT]+|[I])$", gwas[["OA"]])) |   # positive assertion - EA valid SNP or deleteion, OA valid SNP or insertion
                    ((grepl("^([AGCT]+|[D])$", gwas[["OA"]]) | is.na(gwas[["OA"]])) & grepl("^([AGCT]+|[I])$", gwas[["EA"]]))) &  # positive assertion - OA valid SNP or deleteion, EA valid SNP or insertion
                   (!gwas[["EA"]]==gwas[["OA"]])                                                                                  # negative assertion - OA can't be the same as EA

  # report invalid
  if(any(!valid_alleles)) {

    example_invalid <- paste0(gwas[[which(!valid_alleles)[[1]], "EA"]], "/", gwas[[which(!valid_alleles)[[1]], "OA"]])
    message(paste0("[-] ", sum(!valid_alleles), " rows removed due to invalid alleles. Example: '", example_invalid, "'"))

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
populate_rsid <- function(gwas, populate_rsid=FALSE, missing_rsid="none") {

  # rsid regex
  valid_regex <- "^(rs[0-9]+|[0-9XY]+:[0-9]+).*"

  # parse existing RSID if there
  if("RSID" %in% names(gwas)) {

    gwas[, RSID := tolower(sub(valid_regex, "\\1", RSID))]
    invalid <- which(!grepl(valid_regex, gwas$RSID))

  } else {

    invalid <- 1:nrow(gwas)

  }

  # some invalid RSIDs and want to populate from dbSNP
  if(length(invalid) > 0 && populate_rsid!=FALSE) {
    warning("[-] ", length(invalid), " RSIDs could not be parsed, attempting to fetch from dbSNP.")

    # get the RSIDs
    future::plan(future::multisession, workers=parallel::detectCores())
    progressr::with_progress({
      rsid_dat <- chrpos_to_rsid(gwas[invalid, list(CHR,BP,EA,OA)],
                                 chr_col   = "CHR",
                                 pos_col   = "BP",
                                 ea_col    = "EA",
                                 nea_col   = "OA",
                                 build     = populate_rsid,
                                 flip      = "allow",
                                 alt_rsids = FALSE,
                                 verbose   = TRUE)
    })

    # add back
    gwas[invalid, RSID := rsid_dat$RSID]

  # NA out invalid RSIDs
  } else if(length(invalid) > 0 && populate_rsid==FALSE) {
    warning("[-] ", length(invalid), " RSIDs could not be parsed, returning NAs here.")

    gwas[invalid, RSID := NA_character_]

  # everything was fine to begin with
  } # else pass

  # RSIDs at the front
  data.table::setcolorder(gwas, c("RSID", names(gwas)[names(gwas) != "RSID"]))

  # what to do with missing RSIDs; if ignoring, return; else get indices of invalid
  if(missing_rsid=='none') {

    return(gwas)

  } else {

    still_invalid <- which(!grepl(valid_regex, gwas$RSID))

  }

  # fill
  switch(missing_rsid,
         'fill_CHR:BP'            = { gwas[still_invalid, RSID := paste0(CHR,':',BP)               ] },
         'fill_CHR:BP_OA_EA'      = { gwas[still_invalid, RSID := paste0(CHR,':',BP,'_',OA,'_',EA) ] },
         'overwrite_CHR:BP'       = { gwas[             , RSID := paste0(CHR,':',BP)               ] },
         'overwrite_CHR:BP:OA:EA' = { gwas[             , RSID := paste0(CHR,':',BP,'_',OA,'_',EA) ] })

  # return
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

  # TODO: other column health checks

  # check P values; must have P value data
  if("P" %in% names(gwas)) {

    invalid_pvals <- gwas[["P"]] <= 0 | gwas[["P"]] > 1 | is.na(gwas[["P"]])

    if (sum(invalid_pvals, na.rm=TRUE) > 0 | any(is.na(invalid_pvals))) {
      cols <- c("SNP","CHR","BP","EA","OA","EAF","BETA","SE","P")
      warning(paste0("GWAS has ", sum(invalid_pvals), " P value(s) outside accepted range, please check your filter settings."))
    }

  } else {

    # shouldn't get here
    stop("GWAS does not have P value data even after running standardisation procedure")

  }

  # check EAF; don't absolutely have to have EAF data
  if("EAF" %in% names(gwas)) {

    invalid_eaf <- gwas[["EAF"]] < 0 | gwas[["EAF"]] > 1 | is.na(gwas[["EAF"]])

    if (sum(invalid_eaf, na.rm=T) > 0 | any(is.na(invalid_eaf))) {
      cols <- c("SNP","CHR","BP","EA","OA","EAF","BETA","SE","P")
      warning(paste0("GWAS has ", sum(invalid_eaf), " EAF value(s) outside accepted range, please check your filter settings."))
    }

  } else {

    warning("GWAS does not have EAF data")

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

  gwas <- import_table(gwas) |>
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




