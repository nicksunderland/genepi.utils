# Silence R CMD check
globalVariables(c("ind", "seqnames", "start"),
                package = "genepi.utils")

#' @title Liftover GWAS positions
#' @description
#' Determine GWAS build and liftover to required build. This is the same function from the
#' GwasDataImport package, the only difference being that you can specify the build rather
#' than it trying to guess the build (which fails if you are trying to liftover
#' small segments of the genome).
#' @references https://github.com/MRCIEU/GwasDataImport
#' @param gwas a data.table, or file path, chr, pos, snp name, effect allele, non-effect allele columns
#' @param from which build to lift from, one of c("Hg18", "Hg19", "Hg38")
#' @param to which build to lift over to, one of c("Hg18", "Hg19", "Hg38")
#' @param chr_col Name of chromosome column name. Required
#' @param pos_col Name of position column name. Required
#' @param snp_col Name of SNP column name. Optional. Uses less certain method of matching if not available
#' @param ea_col Name of effect allele column name. Optional. Might lead to duplicated rows if not presented
#' @param oa_col Name of other allele column name. Optional. Might lead to duplicated rows if not presented
#' @param remove_duplicates a logical, whether to remove duplicate IDs
#' @return data.table with updated position columns
#' @export
#' @importFrom rtracklayer import.chain liftOver
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#'
lift <- function(gwas,
                 from    = "Hg19",
                 to      = "Hg38",
                 snp_col = "snp",
                 chr_col = "chr",
                 pos_col = "pos",
                 ea_col  = "ea",
                 oa_col  = "oa",
                 remove_duplicates = TRUE)
{

  #TODO: try to remove dependencies GenomicRanges IRanges

  # checks
  builds <- c("Hg18", "Hg19", "Hg38")
  from <- match.arg(from, builds)
  to   <- match.arg(to, builds)
  stopifnot("`to` must not be the same as `from`" = to!=from)
  stopifnot("Column name(s) not found in `gwas`" = all(c(chr_col,pos_col) %in% names(gwas)))

  # import and convert
  gwas <- import_table(gwas)

  # start
  message("Lifting build: ", from, " to ", to)

  # the chain file
  chain_file <- system.file("extdata", paste0(tolower(from), "To", to, ".over.chain"), package="genepi.utils")
  stopifnot("Unable to find chain file" = file.exists(chain_file))

  message("Loading chainfile")
  ch <- rtracklayer::import.chain(chain_file)

  message("Converting chromosome codings")
  if(!grepl("chr", gwas[[1, chr_col]]))
  {
    gwas[, (chr_col) := paste0("chr", get(chr_col))]
  }
  gwas[(chr_col)=="chr23", (chr_col) := "chrX"]
  gwas[(chr_col)=="chr24", (chr_col) := "chrY"]
  gwas[(chr_col)=="chr25", (chr_col) := "chrXY"]
  gwas[(chr_col)=="chr26", (chr_col) := "chrM"]
  gwas[(chr_col)=="chrMT", (chr_col) := "chrM"]


  message("Organising")
  datg <- GenomicRanges::GRanges(
    seqnames = gwas[[chr_col]],
    ranges   = IRanges::IRanges(start=gwas[[pos_col]], end=gwas[[pos_col]]),
    ind      = 1:nrow(gwas)
  )

  message("Lifting")
  d19 <- rtracklayer::liftOver(datg, ch) |> unlist() |> data.table::as.data.table()

  message("Organising again")
  gwas <- gwas[d19$ind, ]
  gwas[, (chr_col) :=  d19$seqnames]
  gwas[, (pos_col) :=  d19$start]
  gwas[, (chr_col) :=  gsub("chr", "", get(chr_col))]

  message("Reordering")
  data.table::setkeyv(gwas, c(chr_col, pos_col))

  if(!is.null(ea_col) & !is.null(oa_col) & remove_duplicates)
  {
    message("Removing duplicates")
    gwas <- unique(gwas, by=c(chr_col, pos_col, ea_col, oa_col))
  }

  message("Done")
  return(gwas)
}

