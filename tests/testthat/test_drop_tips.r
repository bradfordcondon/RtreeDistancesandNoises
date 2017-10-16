source("R/auxillary_functions.R")
library("ape")

test_that("Tips drop", {
start_length <- 10
  x <- rmtree(2, start_length)
  x[[2]]$tip.label[1] <- c("mismatch")  #tree 1 and tree 2 no longer match
  trimmed <- trimTreeToTree(x[[1]], x[[2]])
  expect_equal( class(trimmed), "phylo")
  expect_equal( length(trimmed), (start_length - 1) )
})
