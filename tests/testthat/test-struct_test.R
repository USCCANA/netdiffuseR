context("Structural test")

test_that("struct_test should be reproducible (serial version)", {

  # Generating data
  set.seed(1231)
  graph <- rdiffnet(100, 10)

  set.seed(1313)

  b1 <- struct_test(graph, function(x) mean(threshold(x), na.rm = TRUE), R=100)
  set.seed(1313)
  b2 <- struct_test(graph, function(x) mean(threshold(x), na.rm = TRUE), R=100)

  expect_equal(b1, b2)

})


test_that("struct_test should be reproducible (parallel version)", {

  # I don't want to run this test
  skip_on_cran()
  skip_on_appveyor()

  # Generating data
  set.seed(1231)
  graph <- rdiffnet(150, 5, seed.graph = "small-world")

  # In order to make the parallel version reproducible, we need to set the
  # RNG algorithm to be L'Ecuyer so the sequences can be reproduced

  N <- 100
  set.seed(123, "L'Ecuyer")
  b1 <- struct_test(graph, function(x) mean(threshold(x), na.rm = TRUE), R=N, ncpus=2,
                 parallel="multicore")

  set.seed(123, "L'Ecuyer")
  b2 <- struct_test(graph, function(x) mean(threshold(x), na.rm = TRUE), R=N, ncpus=2,
                 parallel="multicore")

  # Returning to the old RNG. Obs that the seeds are not the same (bug?)
  b1$boot$seed <- 1
  b2$boot$seed <- 1
  b1$graph <- b2$graph <- NULL
  expect_equal(b1$mean_t, b2$mean_t, tolerance=1e-2)

})

# ------------------------------------------------------------------------------
test_that("Methods of struct test", {
  set.seed(1122)
  x <- rdiffnet(100, 4, "central")
  diffnet.toa(x) <- sample(x$toa, nnodes(x))

  ans1 <- struct_test(x, function(g) mean(threshold(g), na.rm = TRUE), 100)
  ans2 <- struct_test(x, function(g) mean(threshold(g), na.rm = TRUE), 100)

  # pvalues of concatenated test should be the same no matter the order
  expect_equal(c(ans1,ans2)$p.value,c(ans2,ans1)$p.value)

  # The number of repetitions should increase
  expect_output(print(c(ans2,ans1,ans1)), "300")

  # Just printing
  expect_silent(hist(ans1, ask=FALSE))
})

# ------------------------------------------------------------------------------
# There is still a lot of work to do in this
test_that("Asymptotic approximations", {
  data(medInnovationsDiffNet)
  g <- medInnovationsDiffNet
  expect_output(print(struct_test_asymp(g$graph[[1]], g$toa, statistic_name = ">=")), "Simulations.+0\n")
})

# ------------------------------------------------------------------------------
test_that("ego_variance", {
  set.seed(9981)
  n <- 100
  g <- rgraph_er(n=100)
  Y <- runif(100)

  # R version
  regovar <- function(graph, Y, funname, all=FALSE) {
    n <- nnodes(graph)


    F_ <- vertex_covariate_compare(graph, cbind(Y), funname)
    F_mean <- mean(outer(Y,Y,function(x,y) abs(x-y)))

    ans <- vector("numeric", n)
    for (i in seq_len(n)) {
      a_ij <- graph[i,,drop=FALSE]@x
      f_ij <- F_[i,,drop=FALSE]@x
      f_i  <- ifelse(all, F_mean, mean(F_[i,,drop=FALSE]@x))

      ans[i] <- sum((f_ij - f_i)^2)/(sum(a_ij) + 1e-15)
    }

    ans
  }

  ans0 <- ego_variance(g, Y, "distance")
  ans1 <- regovar(g,Y,"distance")

  expect_equal(ans0, ans1, tolerance=1e-5)

  ans0 <- ego_variance(g, Y, "distance", TRUE)
  ans1 <- regovar(g,Y,"distance", TRUE)

  expect_equal(ans0, ans1, tolerance=1e-5)

  # microbenchmark(
  #   cpp= ego_variance(g, Y, "distance", TRUE),
  #   r = regovar(g,Y,"distance", TRUE), times=1000, unit="relative"
  # )

})
