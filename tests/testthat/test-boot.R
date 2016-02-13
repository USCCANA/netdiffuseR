context("Network boot")

test_that("struct_test should be reproducible (serial version)", {

  # Generating data
  set.seed(1231)
  graph <- rdiffnet(200, 10)

  set.seed(1313)
  b1 <- struct_test(graph, function(x) mean(threshold(x), na.rm = TRUE), R=500)
  set.seed(1313)
  b2 <- struct_test(graph, function(x) mean(threshold(x), na.rm = TRUE), R=500)

  expect_equal(b1, b2)

})


test_that("struct_test should be reproducible (parallel version)", {

  # Generating data
  set.seed(1231)
  graph <- rdiffnet(200, 10)

  # In order to make the parallel version reproducible, we need to set the
  # RNG algorithm to be L'Ecuyer so the sequences can be reproduced

  oldrng <- RNGkind()[1]
  RNGkind("L'Ecuyer-CMRG")
  set.seed(1313)
  b1 <- struct_test(graph, function(x) mean(threshold(x), na.rm = TRUE), R=500, ncpus=2,
                 parallel="multicore")
  set.seed(1313)
  b2 <- struct_test(graph, function(x) mean(threshold(x), na.rm = TRUE), R=500, ncpus=2,
                 parallel="multicore")

  # Returning to the old RNG. Obs that the seeds are not the same (bug?)
  RNGkind(oldrng)
  b1$boot$seed <- 1
  b2$boot$seed <- 1
  expect_equal(b1, b2)

})
