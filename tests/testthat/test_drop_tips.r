library("ape")
test_that("Tips dropped from reference tree", {
start_length <- 10
  x <- rmtree(2, start_length)
  x[[2]]$tip.label[1] <- c("mismatch")  #tree 1 and tree 2 no longer match
  trimmed <- trimTreeToTree(x[[1]], x[[2]])
  ## Did we return a tree?
  expect_equal( class(trimmed), "phylo")
  # Did the tree get smaller by 1?
  expect_equal( length(trimmed$tip.label), (start_length - 1) )
  ###trimTreeFromList
  expect_equal(length(trimTreeFromList(x[[2]], c("mismatch"))$tip.label), 1)

})
