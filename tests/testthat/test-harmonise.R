# test_that("harmonise", {
#
#   gwas1 <- data.table::data.table(
#     chr = c("1","2","3","4","5","6","7"),
#     bp  = c(1,2,3,4,5,6,7),
#     ea  = c("A",  "C",  "ACC",  "D",  "D",  "AT", "A" ),
#     oa  = c("T",  "A",  "A",    "I",  "I",  "AT", "A" ),
#     eaf = c(0.1,  0.9,  0.2,    0.1,  0.9,  0.1,  0.1 ),
#     beta= c(0.5, -0.2,  0.1,    0.1,  0.1,  0.1,  0.1 ),
#     p   = c(1e-3, 1e-4, 3e-4,   3e-3, 3e-3, 3e-3, 3e-3)
#   )
#   gwas1[, rsid := paste0(chr,":",bp,"[b37]",oa,",",ea)]
#
#   gwas2 <- data.table::data.table(
#     chr = gwas1$chr,
#     bp  = gwas1$bp,
#     ea  = c("A",  "A",  "ACC",  "A",  "AC",  "AT",  "A"),
#     oa  = c("T",  "C",  "A",    "AC", "A",   "AT",  "A"),
#     eaf = c(0.2,  0.1,  0.2,    0.1,  0.1,   0.1,   0.1),
#     beta= c(0.1,  0.2,  0.1,    0.1, -0.2,   0.1,   0.1),
#     p   = c(1e-2, 1e-4, 3e-4,   3e-3, 3e-3,  3e-3,  3e-3)
#   )
#   gwas2[, rsid := paste0(chr,":",bp,"[b37]",oa,",",ea)]
#
#   h <- harmonise(gwas1, gwas2, gwas1_trait="in", gwas2_trait="pr", merge=c("chr"="chr", "bp"="bp"))
#
#   # expected output; scenarios
#   # 1 - no changes
#   # 2 - flip SNP alleles, beta, and EAF
#   # 3 - no change, but for nchar indel
#   # 4 - recode gwas1 D/I coding to alleles of gwas2
#   # 5 - recode gwas1 D/I coding to alleles and also flip gwas2 alleles
#   # 6 - incorrect indel allele, will trying flipping then mark for deletion
#   # 7 - incorrect SNP allele, will trying flipping then mark for deletion
#   h_expected <- data.table::data.table(
#     chr_in = gwas1$chr,
#     chr_pr = gwas2$chr,
#     bp_in  = gwas1$bp,
#     bp_pr  = gwas2$bp,
#     ea_in  = c("A","C","ACC","A","A","AT","A"),
#     ea_pr  = c("A","C","ACC","A","A","AT","A"),
#     oa_in  = c("T","A","A","AC","AC","AT","A"),
#     oa_pr  = c("T","A","A","AC","AC","AT","A"),
#     eaf_in = c(0.1,0.9,0.2,0.1,0.9,0.1,0.1),
#     eaf_pr = c(0.2,0.9,0.2,0.1,0.9,0.1,0.1),
#     beta_in= c(0.5,-0.2,0.1,0.1,0.1,0.1,0.1),
#     beta_pr= c(0.1,-0.2,0.1,0.1,0.2,0.1,0.1),
#     p_in   = c(1e-3,1e-4,3e-4,3e-3,3e-3,3e-3,3e-3),
#     p_pr   = c(1e-2,1e-4,3e-4,3e-3,3e-3,3e-3,3e-3),
#     palindromic = c(T,F,F,F,F,F,F),
#     keep   = c(T,T,T,T,T,F,F)
#   )
#   # need to adjust SNP coding if flipped allele
#   h_expected[, rsid_in := paste0(chr_in,":",bp_in,"[b37]",oa_in,",",ea_in)]
#   h_expected[, rsid_pr := paste0(chr_pr,":",bp_pr,"[b37]",oa_pr,",",ea_pr)]
#   data.table::setcolorder(h_expected, c("rsid_in", "rsid_pr"))
#
#   # don't assess the attributes
#   attr(h, "info") <- NULL
#   attr(h, "sorted") <- NULL
#
#   # test
#   expect_equal(h, h_expected)
#
# })
#
#
#
#
# test_that("flip_alleles", {
#
#   # with a vector
#   input_alleles <- c("A","C","T","G","D","ACCT")
#   output_alleles <- flip_alleles(input_alleles)
#   expect_equal(output_alleles, c("T","G","A","C","D","ACCT"))
#
#   # in a data.table
#   input_alleles_dt <- data.table::data.table(ea = c("A","C","T","G","I","ACCT"))
#   input_alleles_dt[, ea := flip_alleles(ea)]
#
#   # test
#   expect_equal(input_alleles_dt, data.table::data.table(ea = c("T","G","A","C","I","ACCT")))
#
# })
#
#
# test_that("is.palindromic", {
#
#   # single
#   expect_true(is.palindromic("A","T"))
#   expect_true(is.palindromic("T","A"))
#   expect_true(is.palindromic("C","G"))
#   expect_true(is.palindromic("G","C"))
#   expect_false(is.palindromic("G","A"))
#   expect_false(is.palindromic("GACC","A"))
#   expect_false(is.palindromic("D","I"))
#
#   # as data.table
#   input_alleles_dt <- data.table::data.table(ea = c("A","C","T","G","GACC","D"), oa = c("T","G","G","A","A","I"))
#   input_alleles_dt[, palindromic := is.palindromic(ea, oa)]
#   expect_equal(input_alleles_dt$palindromic, c(T,T,F,F,F,F))
#
# })
#
#
# test_that("recode_indels", {
#
#   # input data.table
#   h <- data.table::data.table(ea_1 = c("A","C","D","GT","D","GACC","ACC"),
#                               oa_1 = c("T","G","I","A","I","A","ACC"),
#                               ea_2 = c("A","C","T","I","GACC","D","ACC"),
#                               oa_2 = c("T","G","GA","D","G","I","ACC"),
#                               keep = c(T,T,T,T,T,T,T))
#
#   # recode (n.b. this isn't flipping)
#   h <- recode_indels(h)
#
#   # expected
#   h_expected <- data.table::data.table(ea_1 = c("A","C","T","GT","G","GACC","ACC"),
#                                        oa_1 = c("T","G","GA","A","GACC","A","ACC"),
#                                        ea_2 = c("A","C","T","GT","GACC","A","ACC"),
#                                        oa_2 = c("T","G","GA","A","G","GACC","ACC"),
#                                        keep = c(T,T,T,T,T,T,F))
#
#   # test
#   expect_equal(h, h_expected)
# })
#
#
