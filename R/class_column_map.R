ColumnMapping <- function() {

  # possible maps
  maps <- list(
    default = list(SNP="SNP", CHR="CHR", BP="BP", EA="EA", OA="OA", EAF="EAF", P="P", BETA="BETA", SE="SE", OR="OR", OR_SE="OR_SE", OR_LB="OR_LB", OR_UB="OR_UB", RSID="RSID"),
    metal   = list(SNP="MarkerName", EA="Allele1", OA="Allele2", EAF="Freq1", P="P-value", BETA="Effect", SE="StdErr"),
    ieu_ukb = list(SNP="SNP", BETA="BETA", SE="SE", EA="ALLELE1", OA="ALLELE0", EAF="A1FREQ", P="P_BOLT_LMM_INF")
  )

  # return the requested map as S3 object
  structure(
    .Data = maps,
    class = "ColumnMapping"
  )
}

is.input_format <- function(x, input_format) UseMethod("is.input_format")
is.input_format.ColumnMapping <- function(x, input_format) {
  input_format %in% names(x)
}

input_formats <- function(x) UseMethod("input_formats")
input_formats.ColumnMapping <- function(x) {
  names(x)
}

get_map <- function(x, input_format) UseMethod("get_map")
get_map.ColumnMapping <- function(x, input_format) {

  # checks
  if(!is.input_format(x, input_format)) {
    stop(paste0("ColumnMapping() - `input_format` must be a defined map. Options: ", paste0(names(x), collapse=", ")))
  }

  # return the map
  return(x[[input_format]])
}



