#' #' @title Lift genome positions across genome builds
#' #' @inheritParams chrpos_to_rsid
#' #' @param from a string, genome build code for input `dt` RSID; options: 'b37', 'b38'
#' #' @param to a string, genome build to convert to; options: 'b37', 'b38'
#' #' @param rsid_col a string, RSID column name
#' #' @return an updated data.table with RSID, CHR, and BP for the updated builds
#' #' @noRd
#' #'
#' lift <- function(dt,
#'                  from     = "b37",
#'                  to       = "b38",
#'                  rsid_col = "RSID",
#'                  chr_col  = "CHR",
#'                  pos_col  = "BP",
#'                  dbsnp_dir= genepi.utils::which_dbsnp_directory(),
#'                  verbose = TRUE) {
#'
#'   # import / convert input
#'   dt <- import_table(dt)
#'
#'   stopifnot("`rsid_col` must be a column name in `dt`"= rsid_col %in% colnames(dt))
#'   stopifnot("`chr_col` must be a column name in `dt`" = chr_col %in% colnames(dt))
#'   stopifnot("`pos_col` must be a column name in `dt`" = pos_col %in% colnames(dt))
#'   stopifnot("`pos_col` must be numeric/integer type" = is.numeric(dt[,get(pos_col)]))
#'   stopifnot("`rsid_col` must be character type" = is.character(dt[,get(rsid_col)]))
#'   stopifnot("`chr_col` column must be character, integer, or numeric type" = is.character(dt[,get(chr_col)]) |
#'               is.integer(dt[,get(chr_col)]) |
#'               is.numeric(dt[,get(chr_col)]))
#'   if(is.numeric(dt[,get(chr_col)])) {
#'     stopifnot("`chr_col`, if numeric, must be whole numbers" = all(dt[,get(chr_col)] %% 1 == 0))
#'   }
#'   stopifnot("`dbsnp_dir` must be a valid directory path" = dir.exists(dbsnp_dir))
#'   from <- match.arg(from, c("b37","b38"))
#'   to   <- match.arg(to, c("b37","b38"))
#'   stopifnot("`from` should not be the same as `to`" = from!=to)
#'
#'
#'   #-- Processing ---------------------------------------------------------
#'
#'   # data adjustment - remove chromosome 'Chr6' --> '6' and recode X/Y
#'   # (?i) case-insensitive
#'   # ^ at the beginning of the string
#'   # (?:chr)? non-capture group, matches "chr" 0-1 times
#'   # ([0-9XYMT]+) capture group //1, matches 1 or more occurances of group [0-9XYMT]
#'   # \\1 replace with string found in capture group 1
#'   chr_col_orig = "chr_col_orig"
#'   dt[, (chr_col_orig) := get(chr_col)]
#'   dt[, (chr_col) := sub("(?i)^(?:chr)?([0-9XY]+)", toupper("\\1"), toupper(get(chr_col)))]
#'   # replace X[23] and Y[24]
#'   dt[, (chr_col) := ifelse(get(chr_col)=="X","23",ifelse(get(chr_col)=="Y","24",get(chr_col)))]
#'
#'   # split the dt into a list of data.tables, by chromosome
#'   dts <- split(dt, by=chr_col)
#'
#'   if(verbose) {
#'     message(paste("RSID mapping...\nChromosomes to process: ", paste0(names(dts), collapse=", ")))
#'   }
#'
#'   # set up progress bar
#'   p <- progressr::progressor(steps = 7*length(dts)) # going to update 7 times in the function
#'
#'   # set up parallel processing and process each chromosome independently
#'   out <- furrr::future_map(.x = dts,
#'                            .f = rsid_to_rsid,
#'                            from=from, to=to, rsid_col=rsid_col, chr_col=chr_col, pos_col=pos_col, dbsnp_dir=dbsnp_dir, p=p,
#'                            .options=furrr::furrr_options(seed=TRUE))
#'
#'   out_data <- data.table::rbindlist(out)
#'
#'   # put back the original chromosome column and remove the temporary one
#'   out_data[, (chr_col) := chr_col_orig]
#'   out_data[, chr_col_orig := NULL]
#'
#'   # # report
#'   # msg <- paste0("RSID coverage ", round(100*sum(!is.na(out_data[["RSID"]]))/nrow(out_data), digits=2), "% (", sum(!is.na(out_data[["RSID"]])), "/", nrow(out_data), ")")
#'   # if(flip!="no_flip") {
#'   #   msg <- paste0(msg, ", of which ", round(100*sum(out_data[["rsid_flip_match"]], na.rm=TRUE)/nrow(out_data), digits=2),"% (", sum(out_data[["rsid_flip_match"]], na.rm=TRUE), "/", nrow(out_data), ") where found flipping alleles")
#'   # }
#'   # message(msg)
#'
#'   out_data <- revert_table_type(out_data)
#'   return(out_data)
#' }
#'
#'
#'
#'
#'
#' rsid_to_rsid <- function(chrom_dt, from, to, rsid_col, chr_col, pos_col, dbsnp_dir) {
#'
#'   # silence RMDcheck warning
#'   # RSID = i.RSID = baseRSID = rsid_flip_match = REF = ALT = NULL
#'
#'   # increment progress bar #1
#'   p()
#'
#'   # the chromosome to process, a string value e.g. "1"
#'   chrom <- chrom_dt[[1, chr_col]]
#'
#'   # set the to and from paths
#'   if(from=="b37" & to=="b38") {
#'     from_dbSNP_build_dir  <- file.path(dbsnp_dir, "b37_dbsnp156")
#'     to_dbSNP_build_dir    <- file.path(dbsnp_dir, "b38_dbsnp156")
#'   } else if(from=="b38" & to=="b37") {
#'     from_dbSNP_build_dir  <- file.path(dbsnp_dir, "b38_dbsnp156")
#'     to_dbSNP_build_dir    <- file.path(dbsnp_dir, "b37_dbsnp156")
#'   } else {
#'     stop("Only b37 and b38 builds supported currently")
#'   }
#'
#'   # the dbSNP `.fst` files for this chromosome
#'   to_dbSNP_path   <- file.path(to_dbSNP_build_dir,   paste0("chr", chrom, ".fst"))
#'   stopifnot("Missing .fst reference file [`to`] - is the dbSNP directory set correctly?" = file.exists(to_dbSNP_path))
#'
#'   # read just the RSID data (reduce amount of data read in)
#'   to_dbSNP_rsid <- fst::read_fst(to_dbSNP_path, "RSID", as.data.table=TRUE)
#'
#'   # increment progress bar #2
#'   p()
#'
#'   # set as keys
#'   data.table::setkey(to_dbSNP_rsid, "RSID")
#'   data.table::setkeyv(chrom_dt, rsid_col)
#'
#'   # get the indices of the positions that are needed / also found in the input data
#'   to_dbSNP_rsid[chrom_dt, "found" := get(rsid_col)]
#'   to_row_idxs <- which(!is.na(to_dbSNP_rsid[["found"]]))
#'   rm(to_dbSNP_rsid)
#'
#'   # if matches found, get the from variants
#'   if(length(to_row_idxs)>0) {
#'
#'     # increment progress bar #3
#'     p()
#'
#'     # create a fst object which allows row access without reading the whole file
#'     to_dbSNP_fst <- fst::fst(to_dbSNP_path)
#'
#'     # read the needed rows
#'     to_dbSNP_data <- to_dbSNP_fst[to_row_idxs, c("RSID","CHR","BP")] |> data.table::as.data.table()
#'
#'     # increment progress bar #4
#'     p()
#'
#'     # set the keys on the rsid
#'     data.table::setkey(to_dbSNP_data, "RSID")
#'
#'     # join
#'     chrom_dt[to_dbSNP_data, paste0(c("RSID","CHR","BP"),"_",to) := list(i.RSID,i.CHR,i.BP)]
#'
#'     # no matches found. skip all the processing but add the columns
#'   } else {
#'     p()
#'     p()
#'     p()
#'     p()
#'     chrom_dt[, paste0(c("RSID","CHR","BP"),"_",to) := list(NA_character_,NA_character_,NA_character_)]
#'   }
#'
#'   return(chrom_dt)
#' }



