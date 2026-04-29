context("rdiffnet disadoption mechanisms")
library(netdiffuseR)

# ----------------------------------------------------------------------------
# disadoptmech_random
# ----------------------------------------------------------------------------

test_that("disadoptmech_random returns a function with the right contract", {
  f <- disadoptmech_random(prob = 0.20)
  expect_type(f, "closure")
  expect_named(formals(f), c("expo", "cumadopt", "time"))

  # Hand-call: 100 nodes, all currently adopted, prob = 0.5.
  # -expo- shape mirrors what rdiffnet passes: n x 1 x Q (current slice).
  set.seed(2026)
  cumadopt <- array(1L, dim = c(100, 5, 1))
  expo     <- array(0,  dim = c(100, 1, 1))
  res <- disadoptmech_random(prob = 0.5)(expo, cumadopt, time = 3L)
  expect_length(res, 1L)
  expect_type(res[[1]], "integer")
  # ~50% of 100 nodes should disadopt
  expect_true(abs(length(res[[1]]) - 50L) < 15L)
})

test_that("disadoptmech_random validates -prob-", {
  expect_error(disadoptmech_random(),               "requires -prob-")
  expect_error(disadoptmech_random(prob = -0.1),    "in \\[0, 1\\]")
  expect_error(disadoptmech_random(prob = 1.1),     "in \\[0, 1\\]")
  expect_error(disadoptmech_random(prob = c(0.1, 0.2)), "single number")
  expect_error(disadoptmech_random(prob = NA_real_), "in \\[0, 1\\]")
})

test_that("disadoptmech_random integrates with rdiffnet", {
  set.seed(2026)
  dn <- rdiffnet(n = 50, t = 10, seed.graph = "small-world",
                 seed.p.adopt = 0.10, stop.no.diff = FALSE,
                 disadopt = disadoptmech_random(prob = 0.15))
  expect_s3_class(dn, "diffnet")
})

# ----------------------------------------------------------------------------
# disadoptmech_bithreshold
# ----------------------------------------------------------------------------

test_that("disadoptmech_bithreshold disadopts only currently-adopted with high expo", {
  f <- disadoptmech_bithreshold(threshold_dis = 0.5)

  # rdiffnet hands disadopt() an -expo- of shape n x 1 x Q (current slice).
  # 5 nodes, behaviour 1, time 3.
  cumadopt <- array(0L, dim = c(5, 4, 1))
  cumadopt[c(1, 2, 4), 3, 1] <- 1L                 # nodes 1, 2, 4 currently adopted
  expo     <- array(c(0.7, 0.3, 0.9, 0.6, 0.8), dim = c(5, 1, 1))
  # node 1 and 4 are currently adopted AND cross threshold;
  # node 2 is adopted but expo = 0.3 < 0.5;
  # node 3 crosses but isn't adopted; node 5 isn't adopted.

  res <- f(expo, cumadopt, time = 3L)
  expect_equal(sort(res[[1]]), c(1L, 4L))
})

test_that("disadoptmech_bithreshold validates -threshold_dis-", {
  expect_error(disadoptmech_bithreshold(),                  "requires -threshold_dis-")
  expect_error(disadoptmech_bithreshold(threshold_dis = NA),"NA")
  expect_error(disadoptmech_bithreshold(threshold_dis = "x"),
               "numeric scalar or vector")
})

test_that("disadoptmech_bithreshold integrates with rdiffnet", {
  set.seed(2026)
  dn <- rdiffnet(n = 50, t = 10, seed.graph = "small-world",
                 seed.p.adopt   = 0.10, stop.no.diff = FALSE,
                 threshold.dist = 0.30,
                 disadopt = disadoptmech_bithreshold(threshold_dis = 0.70))
  expect_s3_class(dn, "diffnet")
})

# ----------------------------------------------------------------------------
# disadoptmech_logit / disadoptmech_probit
# ----------------------------------------------------------------------------

test_that("disadoptmech_logit validates parameters", {
  expect_error(disadoptmech_logit(),                "requires both -beta0- and -beta_expo-")
  expect_error(disadoptmech_logit(beta0 = 0),       "requires both")
  expect_error(disadoptmech_logit(beta_expo = 1),   "requires both")
})

test_that("disadoptmech_logit saturates correctly", {
  # Big positive beta0 -> plogis(...) ~ 1 -> all currently-adopted disadopt
  set.seed(2026)
  f <- disadoptmech_logit(beta0 = 50, beta_expo = 0)
  cumadopt <- array(0L, dim = c(20, 3, 1))
  cumadopt[1:10, 2, 1] <- 1L
  expo <- array(0, dim = c(20, 1, 1))           # current-slice shape
  res <- f(expo, cumadopt, time = 2L)
  expect_equal(sort(res[[1]]), 1:10)

  # Big negative beta0 -> plogis(...) ~ 0 -> no disadoption
  g <- disadoptmech_logit(beta0 = -50, beta_expo = 0)
  res2 <- g(expo, cumadopt, time = 2L)
  expect_length(res2[[1]], 0L)
})

test_that("disadoptmech_probit validates and saturates", {
  expect_error(disadoptmech_probit(beta0 = 0),       "requires both")

  set.seed(2026)
  f <- disadoptmech_probit(beta0 = 8, beta_expo = 0)   # pnorm(8) ~ 1
  cumadopt <- array(0L, dim = c(20, 3, 1))
  cumadopt[1:10, 2, 1] <- 1L
  expo <- array(0, dim = c(20, 1, 1))           # current-slice shape
  res <- f(expo, cumadopt, time = 2L)
  expect_equal(sort(res[[1]]), 1:10)
})

test_that("disadoptmech_logit integrates with rdiffnet", {
  set.seed(2026)
  dn <- rdiffnet(n = 50, t = 10, seed.graph = "small-world",
                 seed.p.adopt = 0.10, stop.no.diff = FALSE,
                 disadopt = disadoptmech_logit(beta0 = -1, beta_expo = -2))
  expect_s3_class(dn, "diffnet")
})

# ----------------------------------------------------------------------------
# Bug fix: error messages reference -adoption_pars- (the user-facing arg)
# ----------------------------------------------------------------------------

test_that("adoptmech_logit error message references -adoption_pars-", {
  expect_error(
    rdiffnet(n = 20, t = 4, seed.graph = "small-world",
             adoption_mechanism = adoptmech_logit,
             adoption_pars = list(beta0 = 0),
             stop.no.diff = FALSE),
    "-adoption_pars-"
  )
})

test_that("adoptmech_probit error message references -adoption_pars-", {
  expect_error(
    rdiffnet(n = 20, t = 4, seed.graph = "small-world",
             adoption_mechanism = adoptmech_probit,
             adoption_pars = list(beta_expo = 1),
             stop.no.diff = FALSE),
    "-adoption_pars-"
  )
})

# ----------------------------------------------------------------------------
# Composition: adoption + disadoption mechanisms together
# ----------------------------------------------------------------------------

test_that("adoptmech_logit composes with disadoptmech_random in rdiffnet", {
  set.seed(2026)
  dn <- rdiffnet(n = 50, t = 12, seed.graph = "small-world",
                 seed.p.adopt = 0.05, stop.no.diff = FALSE,
                 adoption_mechanism = adoptmech_logit,
                 adoption_pars      = list(beta0 = -2, beta_expo = 6),
                 disadopt = disadoptmech_random(prob = 0.10))
  expect_s3_class(dn, "diffnet")
  expect_equal(length(dn$toa), 50)
})
