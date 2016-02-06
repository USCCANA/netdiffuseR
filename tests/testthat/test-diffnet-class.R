context("Data-management on diffnet objects")

test_that("Subsetting slices", {
  set.seed(13131)
  diffnet <- rdiffnet(100, 20)

  sum1 <- summary(diffnet, no.print = TRUE)
  sum2 <- summary(diffnet.subset.slices(diffnet, c(7:13)), no.print = TRUE)

  # Number of adopters should hold (adoption rate)
  expect_equal(sum1[nrow(sum1),c("cum_adopt")],sum2[nrow(sum2),c("cum_adopt")])
})

test_that("Extracting attributes", {
  set.seed(009)
  diffnet <- rdiffnet(80, 20, seed.nodes = "random", seed.p.adopt = .1)

  # Testing error messages
  expect_error(diffnet.attrs(diffnet, "vertex", attr.class = "opa"), "should only have")
  expect_error(diffnet.attrs(diffnet, "svertex", attr.class = "dyn"), "should only have")
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
  expect_output(diffnet, "Final prevalence\\s+[:] 1\\.0")

  # No adopters... what!?
  diffnet.toa(diffnet) <- NA
  expect_output(diffnet, "Final prevalence\\s+[:] 0\\.0")
})
