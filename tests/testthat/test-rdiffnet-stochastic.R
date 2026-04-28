context("rdiffnet stochastic adoption mechanism")
library(netdiffuseR)

test_that("default (threshold mechanism) path is unchanged", {
  # A default rdiffnet() call (no adoption_mechanism) and an explicit
  # adoption_mechanism = adoptmech_threshold call with the same seed
  # must produce identical $toa.
  set.seed(2026)
  dn_default <- rdiffnet(n = 25, t = 6, seed.graph = "small-world",
                         seed.p.adopt = 0.1, stop.no.diff = FALSE)

  set.seed(2026)
  dn_thr <- rdiffnet(n = 25, t = 6, seed.graph = "small-world",
                     seed.p.adopt = 0.1, stop.no.diff = FALSE,
                     adoption_mechanism = adoptmech_threshold)

  expect_identical(dn_default$toa, dn_thr$toa)
})

test_that("adoptmech_logit runs and returns a diffnet", {
  set.seed(2026)
  dn <- rdiffnet(n = 40, t = 6, seed.graph = "small-world",
                 seed.p.adopt = 0.05, stop.no.diff = FALSE,
                 adoption_mechanism = adoptmech_logit,
                 adoption_pars = list(beta0 = -2, beta_expo = 6))

  expect_s3_class(dn, "diffnet")
  expect_equal(dim(dn$cumadopt), c(40, 6))
  expect_equal(length(dn$toa), 40)
})

test_that("adoptmech_logit requires both beta0 and beta_expo", {
  expect_error(
    rdiffnet(n = 20, t = 4, seed.graph = "small-world",
             adoption_mechanism = adoptmech_logit,
             adoption_pars = list(beta0 = 0),
             stop.no.diff = FALSE),
    "beta0.*beta_expo|beta_expo"
  )
  expect_error(
    rdiffnet(n = 20, t = 4, seed.graph = "small-world",
             adoption_mechanism = adoptmech_logit,
             adoption_pars = list(beta_expo = 1),
             stop.no.diff = FALSE),
    "beta0"
  )
  expect_error(
    rdiffnet(n = 20, t = 4, seed.graph = "small-world",
             adoption_mechanism = adoptmech_logit,
             stop.no.diff = FALSE),
    "beta0"
  )
})

test_that("adoptmech_probit runs and validates pars", {
  set.seed(2026)
  dn <- rdiffnet(n = 30, t = 5, seed.graph = "small-world",
                 seed.p.adopt = 0.05, stop.no.diff = FALSE,
                 adoption_mechanism = adoptmech_probit,
                 adoption_pars = list(beta0 = -1, beta_expo = 3))
  expect_s3_class(dn, "diffnet")

  expect_error(
    rdiffnet(n = 20, t = 4, seed.graph = "small-world",
             adoption_mechanism = adoptmech_probit,
             stop.no.diff = FALSE),
    "beta0"
  )
})

test_that("saturating beta0 drives near-universal adoption (logit)", {
  # Very large intercept -> plogis(beta0 + ...) ~ 1 for all exposures.
  set.seed(99)
  dn <- rdiffnet(n = 60, t = 8, seed.graph = "small-world",
                 seed.p.adopt = 0.05, stop.no.diff = FALSE,
                 adoption_mechanism = adoptmech_logit,
                 adoption_pars = list(beta0 = 50, beta_expo = 0))
  # Everyone has adopted by some t <= T
  expect_true(all(!is.na(dn$toa)))
})

test_that("very negative beta0 + beta_expo = 0 suppresses diffusion (logit)", {
  set.seed(99)
  expect_warning(
    dn <- rdiffnet(n = 30, t = 6, seed.graph = "small-world",
                   seed.p.adopt = 0.05, stop.no.diff = FALSE,
                   adoption_mechanism = adoptmech_logit,
                   adoption_pars = list(beta0 = -50, beta_expo = 0)),
    "No diffusion"
  )
  # All adopters come from the seed round (toa == 1)
  expect_true(all(dn$toa[!is.na(dn$toa)] == 1L))
})

test_that("logit mechanism works with multiple behaviors", {
  set.seed(2026)
  dn <- rdiffnet(n = 40, t = 6, seed.graph = "small-world",
                 seed.p.adopt = list(0.05, 0.05),
                 stop.no.diff = FALSE,
                 adoption_mechanism = adoptmech_logit,
                 adoption_pars = list(beta0 = -1, beta_expo = 4))

  expect_s3_class(dn, "diffnet")
  # Multi-behavior diffnet stores cumadopt as a list of length num_behaviors
  expect_type(dn$cumadopt, "list")
  expect_length(dn$cumadopt, 2L)
  expect_equal(dim(dn$toa), c(40, 2))
})

test_that("non-function adoption_mechanism raises a clear error", {
  expect_error(
    rdiffnet(n = 20, t = 4, seed.graph = "small-world",
             adoption_mechanism = "logit", stop.no.diff = FALSE),
    "must be a function"
  )
  expect_error(
    rdiffnet(n = 20, t = 4, seed.graph = "small-world",
             adoption_mechanism = 42, stop.no.diff = FALSE),
    "must be a function"
  )
})

test_that("user-defined mechanism can be plugged in", {
  # A "always-adopt" mechanism: ignore exposure, adopt every available node.
  always_adopt <- function(expo, thresholds, not_adopted, time, pars) {
    which(not_adopted)
  }
  set.seed(2026)
  dn <- rdiffnet(n = 20, t = 4, seed.graph = "small-world",
                 seed.p.adopt = 0.1, stop.no.diff = FALSE,
                 adoption_mechanism = always_adopt)
  # By t = 2 every node should have adopted
  expect_true(all(!is.na(dn$toa)))
  expect_true(all(dn$toa <= 2L))
})
