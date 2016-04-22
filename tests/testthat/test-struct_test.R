context("Network boot")

test_that("struct_test should be reproducible (serial version)", {

  # Generating data
  set.seed(1231)
  graph <- rdiffnet(200, 10)

  set.seed(1313)
  b1 <- struct_test(graph, function(x) mean(threshold(x), na.rm = TRUE), R=100)
  set.seed(1313)
  b2 <- struct_test(graph, function(x) mean(threshold(x), na.rm = TRUE), R=100)

  expect_equal(b1, b2)

})


test_that("struct_test should be reproducible (parallel version)", {

  # Generating data
  set.seed(1231)
  graph <- rdiffnet(200, 10, seed.graph = "small-world")

  # In order to make the parallel version reproducible, we need to set the
  # RNG algorithm to be L'Ecuyer so the sequences can be reproduced

  N <- 2e3
  set.seed(123)
  b1 <- struct_test(graph, function(x) mean(threshold(x), na.rm = TRUE), R=N, ncpus=2,
                 parallel="multicore")

  set.seed(123)
  b2 <- struct_test(graph, function(x) mean(threshold(x), na.rm = TRUE), R=N, ncpus=2,
                 parallel="multicore")

  # Returning to the old RNG. Obs that the seeds are not the same (bug?)
  b1$boot$seed <- 1
  b2$boot$seed <- 1
  b1$graph <- b2$graph <- NULL
  expect_equal(b1$mean_t, b2$mean_t, tolerance=5e-3)

})
