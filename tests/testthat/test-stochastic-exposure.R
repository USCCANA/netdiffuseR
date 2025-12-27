context("Stochastic Exposure")
library(netdiffuseR)

test_that("Stochastic vs Deterministic Exposure", {
  # Create a small random graph
  set.seed(123)
  n <- 10
  g <- rgraph_er(n, t=1, p=0.5)
  # Make it valued with probabilities
  g@x <- runif(length(g@x))

  # Wrap in list to make it "dynamic" (1 time step)
  g_list <- list(g)

  # Create dummy adoption
  cumadopt <- matrix(0, nrow=n, ncol=1)
  cumadopt[1:5, 1] <- 1

  # exposure expects a list for dynamic graphs if not diffnet
  # MUST set valued=TRUE so weights are used as probabilities
  e_det <- exposure(g_list, cumadopt, mode="deterministic", valued=TRUE)
  e_stoch <- exposure(g_list, cumadopt, mode="stochastic", valued=TRUE)

  # They should be different (with high probability)
  expect_false(isTRUE(all.equal(e_det, e_stoch)))

  # Check that stochastic exposure is non-negative and <= 1
  # With the new logic (preserving weights), it should be <= 1
  expect_true(all(e_stoch >= 0))
  expect_true(all(e_stoch <= 1 + 1e-10))
})

test_that("rdiffnet wrapper with stochastic exposure", {
  # We need a dynamic graph for rdiffnet usually, or t > 1
  set.seed(1231)
  diffnet <- rdiffnet(n=20, t=5, seed.graph="small-world", 
                      exposure.mode="stochastic",
                      seed.p.adopt = 0.2,
                      stop.no.diff = FALSE) # Prevent error if no diffusion occurs

  expect_s3_class(diffnet, "diffnet")
  expect_true(all(diffnet$exposure >= 0))
})

test_that("rdiffnet wrapper with stochastic exposure (multiple behaviors)", {
  set.seed(1231)
  # 2 behaviors
  diffnet <- rdiffnet(n=20, t=5, seed.graph="small-world", 
                      exposure.mode="stochastic",
                      seed.p.adopt = list(0.2, 0.2),
                      stop.no.diff = FALSE)

  expect_s3_class(diffnet, "diffnet")
  
  # Calculate exposure to check dimensions
  # Note: We must use the same mode to get consistent dimensions/behavior
  expo <- exposure(diffnet, mode="stochastic")
  
  # Check dimensions of exposure: n x t x 2
  expect_equal(dim(expo), c(20, 5, 2))
  expect_true(all(expo >= 0))
  expect_true(all(expo <= 1 + 1e-10))
})
