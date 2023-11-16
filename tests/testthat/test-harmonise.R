test_that("harmonise", {

  gwas1 <- data.table::data.table(
    CHR = c("1","2","3","4","5","6","7"),
    BP  = c(1,2,3,4,5,6,7),
    EA  = c("A",  "C",  "ACC",  "D",  "D",  "AT", "A" ),
    OA  = c("T",  "A",  "A",    "I",  "I",  "AT", "A" ),
    EAF = c(0.1,  0.9,  0.2,    0.1,  0.9,  0.1,  0.1 ),
    BETA= c(0.5, -0.2,  0.1,    0.1,  0.1,  0.1,  0.1 ),
    P   = c(1e-3, 1e-4, 3e-4,   3e-3, 3e-3, 3e-3, 3e-3)
  )
  gwas1[, SNP := paste0(CHR,":",BP,"[b37]",OA,",",EA)]

  gwas2 <- data.table::data.table(
    CHR = gwas1$CHR,
    BP  = gwas1$BP,
    EA  = c("A",  "A",  "ACC",  "A",  "AC",  "AT",  "A"),
    OA  = c("T",  "C",  "A",    "AC", "A",   "AT",  "A"),
    EAF = c(0.2,  0.1,  0.2,    0.1,  0.1,   0.1,   0.1),
    BETA= c(0.1,  0.2,  0.1,    0.1, -0.2,   0.1,   0.1),
    P   = c(1e-2, 1e-4, 3e-4,   3e-3, 3e-3,  3e-3,  3e-3)
  )
  gwas2[, SNP := paste0(CHR,":",BP,"[b37]",OA,",",EA)]

  h <- harmonise(gwas1, gwas2, gwas1_trait="in", gwas2_trait="pr", merge="chrpos")

  # expected output; scenarios
  # 1 - no changes
  # 2 - flip SNP alleles, beta, and EAF
  # 3 - no change, but for nchar indel
  # 4 - recode gwas1 D/I coding to alleles of gwas2
  # 5 - recode gwas1 D/I coding to alleles and also flip gwas2 alleles
  # 6 - incorrect indel allele, will trying flipping then mark for deletion
  # 7 - incorrect SNP allele, will trying flipping then mark for deletion
  h_expected <- data.table::data.table(
    CHR_in = gwas1$CHR,
    CHR_pr = gwas2$CHR,
    BP_in  = gwas1$BP,
    BP_pr  = gwas2$BP,
    EA_in  = c("A","C","ACC","A","A","AT","A"),
    EA_pr  = c("A","C","ACC","A","A","AT","A"),
    OA_in  = c("T","A","A","AC","AC","AT","A"),
    OA_pr  = c("T","A","A","AC","AC","AT","A"),
    EAF_in = c(0.1,0.9,0.2,0.1,0.9,0.1,0.1),
    EAF_pr = c(0.2,0.9,0.2,0.1,0.9,0.1,0.1),
    BETA_in= c(0.5,-0.2,0.1,0.1,0.1,0.1,0.1),
    BETA_pr= c(0.1,-0.2,0.1,0.1,0.2,0.1,0.1),
    P_in   = c(1e-3,1e-4,3e-4,3e-3,3e-3,3e-3,3e-3),
    P_pr   = c(1e-2,1e-4,3e-4,3e-3,3e-3,3e-3,3e-3),
    palindromic = c(T,F,F,F,F,F,F),
    keep   = c(F,T,T,T,T,F,F)
  )
  # need to adjust SNP coding if flipped allele
  h_expected[, SNP_in := paste0(CHR_in,":",BP_in,"[b37]",OA_in,",",EA_in)]
  h_expected[, SNP_pr := paste0(CHR_pr,":",BP_pr,"[b37]",OA_pr,",",EA_pr)]
  data.table::setcolorder(h_expected, c("SNP_in", "SNP_pr"))

  # don't assess the attributes
  attr(h, "info") <- NULL
  attr(h, "sorted") <- NULL

  # test
  expect_equal(h, h_expected)

})




test_that("flip_alleles", {

  # with a vector
  input_alleles <- c("A","C","T","G","D","ACCT")
  output_alleles <- flip_alleles(input_alleles)
  expect_equal(output_alleles, c("T","G","A","C","D","ACCT"))

  # in a data.table
  input_alleles_dt <- data.table::data.table(EA = c("A","C","T","G","I","ACCT"))
  input_alleles_dt[, EA := flip_alleles(EA)]

  # test
  expect_equal(input_alleles_dt, data.table::data.table(EA = c("T","G","A","C","I","ACCT")))

})


test_that("is.palindromic", {

  # single
  expect_true(is.palindromic("A","T"))
  expect_true(is.palindromic("T","A"))
  expect_true(is.palindromic("C","G"))
  expect_true(is.palindromic("G","C"))
  expect_false(is.palindromic("G","A"))
  expect_false(is.palindromic("GACC","A"))
  expect_false(is.palindromic("D","I"))

  # as data.table
  input_alleles_dt <- data.table::data.table(EA = c("A","C","T","G","GACC","D"), OA = c("T","G","G","A","A","I"))
  input_alleles_dt[, palindromic := is.palindromic(EA, OA)]
  expect_equal(input_alleles_dt$palindromic, c(T,T,F,F,F,F))

})


test_that("recode_indels", {

  # input data.table
  h <- data.table::data.table(EA_1 = c("A","C","D","GT","D","GACC","ACC"),
                              OA_1 = c("T","G","I","A","I","A","ACC"),
                              EA_2 = c("A","C","T","I","GACC","D","ACC"),
                              OA_2 = c("T","G","GA","D","G","I","ACC"),
                              keep = c(T,T,T,T,T,T,T))

  # recode (n.b. this isn't flipping)
  h <- recode_indels(h)

  # expected
  h_expected <- data.table::data.table(EA_1 = c("A","C","T","GT","G","GACC","ACC"),
                                       OA_1 = c("T","G","GA","A","GACC","A","ACC"),
                                       EA_2 = c("A","C","T","GT","GACC","A","ACC"),
                                       OA_2 = c("T","G","GA","A","G","GACC","ACC"),
                                       keep = c(T,T,T,T,T,T,F))

  # test
  expect_equal(h, h_expected)
})


