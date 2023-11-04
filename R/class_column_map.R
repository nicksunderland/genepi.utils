#' @title ColumnMapping class
#' @return an S3 object of class ColumnMapping
#' @export
#'
ColumnMapping <- function() {

  # possible maps
  maps <- list(
    default = list(SNP="SNP", CHR="CHR", BP="BP", EA="EA", OA="OA", EAF="EAF", P="P", BETA="BETA", SE="SE", OR="OR", OR_SE="OR_SE", OR_LB="OR_LB", OR_UB="OR_UB", RSID="RSID"),
    metal   = list(SNP="MarkerName", EA="Allele1", OA="Allele2", EAF="Freq1", P="P-value", BETA="Effect", SE="StdErr"),
    ieu_ukb = list(SNP="SNP", BETA="BETA", SE="SE", EA="ALLELE1", OA="ALLELE0", EAF="A1FREQ", P="P_BOLT_LMM_INF"),
    ns_map  = list(SNP="MARKER", CHR="CHR", BP="POS", EA="A1", OA="A2", BETA="BETA", SE="SE", EAF="EAF", P="P")
  )

  # return the requested map as S3 object
  structure(
    .Data = maps,
    class = "ColumnMapping"
  )
}


#' @title Check if valid input format
#' @param x a ColumnMapping object
#' @param input_format a string, must be a valid input_format
#' @return a logical
#' @export
is.input_format <- function(x, input_format) UseMethod("is.input_format")
#' @rdname is.input_format
#' @export
is.input_format.ColumnMapping <- function(x, input_format) {
  input_format %in% names(x)
}


#' @title Get valid input formats
#' @param x a ColumnMapping object
#' @return a character vector of valid input formats
#' @export
input_formats <- function(x) UseMethod("input_formats")
#' @rdname input_formats
#' @export
input_formats.ColumnMapping <- function(x) {
  names(x)
}


#' @title Get column mapping
#' @param x a ColumnMapping object
#' @param input_format a string, must be a valid input_format
#' @return a named list (standardised_name: original_name)
#' @export
get_map <- function(x, input_format) UseMethod("get_map")
#' @rdname get_map
#' @export
get_map.ColumnMapping <- function(x, input_format) {

  # checks
  if(!is.input_format(x, input_format)) {
    stop(paste0("ColumnMapping() - `input_format` must be a defined map. Options: ", paste0(names(x), collapse=", ")))
  }

  # return the map
  return(x[[input_format]])
}



