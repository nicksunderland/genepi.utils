library(testthat)
library(data.table)

test_that("harmonise_gwas works", {

  # data files
  gwas_file <- test_path("test_data", "harm_test_input.tsv")
  ref_file  <- test_path("test_data", "harm_test_input_ref.tsv")
  out_file  <- test_path("test_data", "harm_test_output.tsv")

  # read test data
  gwas <- fread(gwas_file)
  ref  <- fread(ref_file)
  out  <- fread(out_file)[, c("chr", "chr_ref") := lapply(.SD, as.character), .SDcols = c("chr", "chr_ref")]
  setkey(out, chr, bp)

  # run harmonise function
  h <- harmonise_gwas(gwas, ref, join = "chr:bp", action = 2, rmap = c(chr = "#CHROM", bp = "POS", ea = "ALT", oa = "REF"))

  # compare
  expect_equal(h, out)
})
