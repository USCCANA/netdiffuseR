context("rdiffnet logit adoption model")
library(netdiffuseR)

test_that("threshold (default) path is unchanged", {
  # Compare a threshold run under fixed seed against itself on a second
  # call with adoption_model explicitly set to "threshold".
  set.seed(2026)
  dn_default <- rdiffnet(n = 25, t = 6, seed.graph = "small-world",
                         seed.p.adopt = 0.1, stop.no.diff = FALSE)

  set.seed(2026)
  dn_threshold <- rdiffnet(n = 25, t = 6, seed.graph = "small-world",
                           seed.p.adopt = 0.1, stop.no.diff = FALSE,
                           adoption_model = "threshold")

  expect_identical(dn_default$toa, dn_threshold$toa)
})

test_that("logit mode runs and returns a diffnet", {
  set.seed(2026)
  dn <- rdiffnet(n = 40, t = 6, seed.graph = "small-world",
                 seed.p.adopt = 0.05, stop.no.diff = FALSE,
                 adoption_model = "logit",
                 adoption_pars = list(beta0 = -2, beta_expo = 6))

  expect_s3_class(dn, "diffnet")
  expect_equal(dim(dn$cumadopt), c(40, 6))
  expect_equal(length(dn$toa), 40)
})

test_that("logit mode requires both beta0 and beta_expo", {
  expect_error(
    rdiffnet(n = 20, t = 4, seed.graph = "small-world",
             adoption_model = "logit", adoption_pars = list(beta0 = 0),
             stop.no.diff = FALSE),
    "beta0.*beta_expo|beta_expo"
  )
  expect_error(
    rdiffnet(n = 20, t = 4, seed.graph = "small-world",
             adoption_model = "logit", adoption_pars = list(beta_expo = 1),
             stop.no.diff = FALSE),
    "beta0"
  )
  expect_error(
    rdiffnet(n = 20, t = 4, seed.graph = "small-world",
             adoption_model = "logit",
             stop.no.diff = FALSE),
    "beta0"
  )
})

test_that("saturating beta0 drives near-universal adoption", {
  # Very large intercept -> plogis(beta0 + ...) ~ 1 for all exposures.
  set.seed(99)
  dn <- rdiffnet(n = 60, t = 8, seed.graph = "small-world",
                 seed.p.adopt = 0.05, stop.no.diff = FALSE,
                 adoption_model = "logit",
                 adoption_pars = list(beta0 = 50, beta_expo = 0))
  # Everyone has adopted by some t <= T
  expect_true(all(!is.na(dn$toa)))
})

test_that("very negative beta0 + beta_expo = 0 suppresses diffusion", {
  set.seed(99)
  expect_warning(
    dn <- rdiffnet(n = 30, t = 6, seed.graph = "small-world",
                   seed.p.adopt = 0.05, stop.no.diff = FALSE,
                   adoption_model = "logit",
                   adoption_pars = list(beta0 = -50, beta_expo = 0)),
    "No diffusion"
  )
  # All adopters come from the seed round (toa == 1)
  expect_true(all(dn$toa[!is.na(dn$toa)] == 1L))
})

test_that("logit works with multiple behaviors", {
  set.seed(2026)
  dn <- rdiffnet(n = 40, t = 6, seed.graph = "small-world",
                 seed.p.adopt = list(0.05, 0.05),
                 stop.no.diff = FALSE,
                 adoption_model = "logit",
                 adoption_pars = list(beta0 = -1, beta_expo = 4))

  expect_s3_class(dn, "diffnet")
  # Multi-behavior diffnet stores cumadopt as a list of length num_behaviors
  expect_type(dn$cumadopt, "list")
  expect_length(dn$cumadopt, 2L)
  expect_equal(dim(dn$toa), c(40, 2))
})

test_that("unknown adoption_model raises match.arg error", {
  expect_error(
    rdiffnet(n = 20, t = 4, seed.graph = "small-world",
             adoption_model = "probit", stop.no.diff = FALSE),
    "should be one of"
  )
})
