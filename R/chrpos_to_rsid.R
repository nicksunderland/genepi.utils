#' chrpos_to_rsid
#'
#' @param dt a data.table with at least columns (chrom<chr>, pos<int>, ea<chr>, nea<chr>)
#' @param chr_col a string column name; chromosome position
#' @param pos_col a string column name; base position
#' @param ea_col a string column name; effect allele
#' @param nea_col a string column name; non effect allele
#' @param build a string, options: "GRCh37"
#' @param flip a string, options: "report", "allow", "no_flip"
#'
#' @return a data.table with an RSID column
#' @export
#' @importFrom data.table setkey setkeyv
#' @importFrom progressr progressor
#' @importFrom furrr future_map
#' @importFrom fst fst read_fst
#'
chrpos_to_rsid <- function(dt, chr_col, pos_col, ea_col=NULL, nea_col=NULL, build="GRCh37", flip="report") {

  # checks
  stopifnot("`dt` must be a data.table" = inherits(dt, "data.table"))
  stopifnot("`chr_col` must be a column name in `dt`" = chr_col %in% colnames(dt))
  stopifnot("`pos_col` must be a column name in `dt`" = pos_col %in% colnames(dt))
  stopifnot("`ea_col` must be a column name in `dt`" = ea_col %in% colnames(dt))
  stopifnot("`nea_col` must be a column name in `dt`" = nea_col %in% colnames(dt))
  stopifnot("`chr_col` column must be character type" = is.character(dt[,get(chr_col)]))
  stopifnot("`pos_col` must be numeric/integer type" = is.numeric(dt[,get(pos_col)]))
  stopifnot("`ea_col` column must be character type" = is.character(dt[,get(ea_col)]))
  stopifnot("`nea_col` column must be character type" = is.character(dt[,get(nea_col)]))

  build <- match.arg(build, c("GRCh37"))
  build <- match.arg(flip, c("report", "allow", "no_flip"))

  # data adjustment
  if(any(grepl("(?i)^chr", dt[,get(chr_col)]))) {

    # (?i) case-insensitive
    # ^ at the beginning of the string
    # (?:chr)? non-capture group, matches "chr" 0-1 times
    # ([0-9XYMT]+) capture group //1, matches 1 or more occurances of group [0-9XYMT]
    # \\1 replace with string found in capture group 1
    dt[, "chr_col_orig" := get(chr_col)]
    dt[, (chr_col) := sub("(?i)^(?:chr)?([0-9XYMT]+)", "\\1", get(chr_col))]

  }

  # split the dt into chromosomes
  dts <- split(dt, by=chr_col)

  # set up progress bar
  p <- progressr::progressor(steps = 6*length(dts)) # going to update 5 times in the function

  # set up parallel
  furrr::future_map(dts, ~{
    # the datatable is accessed as `.x`

    # the chromosome to process
    chrom <- .x[[1, chr_col]]

    # the dbSNP file for this chromosome
    # system.file(...., package=Utils)
    dbSNP_path <- paste0("/Users/xx20081/Documents/local_data/genome_reference/rsid_b37_dbsnp155/chr", chrom, ".fst")

    # read just the dbSNP position data (reduce amount of data read in)
    dbSNP_pos <- fst::read_fst(dbSNP_path, "BP", as.data.table=TRUE)

    # increment progress bar #1
    p()

    # set as keys
    data.table::setkey(dbSNP_pos, "BP")
    data.table::setkeyv(.x, pos_col)

    # get the indices of the positions that are needed / also found in the input data
    dbSNP_pos[.x, "found" := get(pos_col)]
    row_idxs <- which(!is.na(dbSNP_pos[["found"]]))
    rm(dbSNP_pos)

    # increment progress bar #2
    p()

    # create a fst object which allows row access without reading the whole file
    dbSNP_fst <- fst::fst(dbSNP_path)

    # read the needed rows
    dbSNP_data <- dbSNP_fst[row_idxs, ] |> data.table::as.data.table()

    # increment progress bar #3
    p()

    # split the ALT column which can be a comma separate vector of alternative alleles
    dbSNP_data[, ALT := lapply(.SD, strsplit, split=','), .SDcols="ALT"]

    # make data.table longer, one allele combination per row
    dbSNP_data <- dbSNP_data[, lapply(.SD, unlist), by=1:nrow(dbSNP_data)]
    dbSNP_data[, nrow := NULL]

    # increment progress bar #4
    p()

    # get the alt rsids (some rsids code for the same position, chr, alt, and ref...)
    alt_rsids  <- dbSNP_data[duplicated(dbSNP_data[, c("CHR","BP","REF","ALT")]),]

    # take unique (first rsid occurance)
    dbSNP_data <- dbSNP_data[!duplicated(dbSNP_data[, c("CHR","BP","REF","ALT")]),]

    # set the keys to match expected way round
    data.table::setkeyv(dbSNP_data, c("CHR","BP","REF","ALT"))
    data.table::setkeyv(.x, c(chr_col, pos_col, nea_col, ea_col))
    .x[dbSNP_data, RSID := i.RSID]

    # increment progress bar # 5
    p()

    # flip the alleles
    if(flip!="no_flip") {

      # set the keys to flipped alleles
      data.table::setkey(dbSNP_data,"CHR","BP","ALT","REF")

      if(flip=="report") {

        # add in flipped matches and a logical flag
        .x[dbSNP_data, c("RSID", "flip_match") := .(i.RSID, TRUE)]

      } else if(flip=="allow") {

        # add in flipped matches
        .x[dbSNP_data, RSID := i.RSID]

      } else {
        stop("internal code problem with `flip` argument")
      }

    }

    # increment progress bar # 6
    p()

    # return
    return(.x)

  }, future.seed=100) |> data.table::rbindlist()

}
