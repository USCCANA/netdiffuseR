context("Data-management on diffnet objects")

test_that("Subsetting slices", {
  set.seed(13131)
  diffnet <- rdiffnet(100, 20)

  sum1 <- summary(diffnet, no.print = TRUE)
  sum2 <- summary(diffnet.subset.slices(diffnet, c(7:13)), no.print = TRUE)

  # Number of adopters should hold (adoption rate)
  expect_equal(sum1[nrow(sum1),c("cum_adopt")],sum2[nrow(sum2),c("cum_adopt")])
})

test_that("Error messages", {
  set.seed(009)
  diffnet <- rdiffnet(80, 20, seed.nodes = "random", seed.p.adopt = .1)

  # Invalid attr.class
  expect_error(diffnet.attrs(diffnet, "vertex", attr.class = "opa"), "should only have")
  expect_error(diffnet.attrs(diffnet, "svertex", attr.class = "dyn"), "should only have")

  # Already a diffnet
  expect_message(as_diffnet(diffnet), "already.+diffnet")

  # Must be data frame or matrix
  static_attr <- diffnet[["real_threshold"]]
  expect_error(
    as_diffnet(diffnet$graph, toa=diffnet.toa(diffnet),
               vertex.static.attrs = static_attr), "must be either a data\\.frame")

  # Coercing to data.frame
  static_attr <- matrix(static_attr, ncol=1)
  expect_warning(
    as_diffnet(diffnet$graph, toa=diffnet.toa(diffnet),
               t0=1, t1=20,
               vertex.static.attrs = static_attr), "coerced to a data\\.frame")
})

test_that("Setting attributes", {
  set.seed(909)
  diffnet <- rdiffnet(80, 20, seed.nodes = "random", seed.p.adopt = .1)

  expect_error(diffnet.attrs(diffnet) <- as.matrix(sample(1:4, 20, TRUE)), "different lengths")
})

test_that("Changing toa", {
  set.seed(18231)
  diffnet <- rdiffnet(100, 10)

  # All to the first time period
  diffnet.toa(diffnet) <- 1
  expect_output(print(diffnet), "Final prevalence\\s+[:] 1\\.0")

  # No adopters... what!?
  diffnet.toa(diffnet) <- NA
  expect_output(print(diffnet), "Final prevalence\\s+[:] 0\\.0")
})

# Checking different input classes ---------------------------------------------
test_that("as_diffnet with different graph classes", {
  # Random diffnet
  set.seed(881)
  graph0 <- rdiffnet(100, 10)

  # Getting the bits
  g_array <- as.array(graph0)
  g_dgcMa <- graph0$graph
  toa     <- diffnet.toa(graph0)
  s_attrb <- data.frame(x=graph0[["real_threshold"]])

  # Using different inputs
  dn_arr <- as_diffnet(g_array, toa=toa, vertex.static.attrs = s_attrb, t0=1, t1=10)
  dn_dgc <- as_diffnet(g_dgcMa, toa=toa, vertex.static.attrs = s_attrb, t0=1, t1=10)

  # Comparing adjmats
  test <- sapply(seq_len(10), function(x) identical(as.matrix(dn_arr$graph[[x]]), as.matrix(dn_dgc$graph[[x]])))
  expect_true(all(test))

  # Comparing the rest
  dn_arr$graph      <- dn_dgc$graph      <- NULL
  dn_arr$meta$class <- dn_dgc$meta$class <- NULL
  expect_equal(dn_arr, dn_dgc)

})
