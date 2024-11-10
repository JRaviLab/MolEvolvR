# Load testthat and the script to test
library(testthat)

# Test for reverseOperonSeq function
test_that("reverseOperonSeq reverses directions correctly", {
  prot <- data.frame(GenContext = c("A>B", "C<D", "E=F*G", "H>I"))
  reversed_prot <- reverseOperonSeq(prot)

  # Expected output
  expected <- data.frame(GenContext = c("A<-B", "C->D", "F*G=E", "I<-H"))

  expect_equal(reversed_prot$GenContext, expected$GenContext)
})

# Test for straightenOperonSeq function
test_that("straightenOperonSeq handles equal signs correctly", {
  genomic_context <- c("A", "B", "*", "C", "D", "=", "E", "F")
  result <- straightenOperonSeq(genomic_context)

  # Expected output after processing
  expected <- c("A->", "B->", "*", "<-C", "<-D", "=", "E->", "F->")

  expect_equal(result, expected)
})
