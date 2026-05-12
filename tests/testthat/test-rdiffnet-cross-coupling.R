context("rdiffnet per-behaviour dispatch and cross-state visibility (Toy B)")
library(netdiffuseR)

# ----------------------------------------------------------------------------
# Contract extension: rdiffnet only passes -behavior- / -expo_all- to the
# adoption_mechanism if the function declares them (or -...-). M6 kernels
# (no extras in their formals) continue to work unchanged.
# ----------------------------------------------------------------------------

test_that("M6 kernels with no -behavior- / -expo_all- formals still work", {
  set.seed(2026)
  dn <- rdiffnet(n = 30, t = 5, seed.graph = "small-world",
                 seed.p.adopt = 0.10, stop.no.diff = FALSE,
                 adoption_mechanism = adoptmech_logit,
                 adoption_pars      = list(beta0 = -2, beta_expo = 5))
  expect_s3_class(dn, "diffnet")
})

# ----------------------------------------------------------------------------
# Per-behaviour dispatch: the mechanism sees -behavior- when it declares it
# ----------------------------------------------------------------------------

test_that("adoption_mechanism with -behavior- formal receives the current q", {
  seen_behaviors <- integer(0)

  per_behavior_mech <- function(expo, thresholds, not_adopted, time, pars,
                                behavior) {
    # Record which behaviour was passed for this invocation
    seen_behaviors <<- c(seen_behaviors, behavior)
    # And branch trivially on it (use threshold for q=1, logit for q=2)
    if (behavior == 1L) {
      which(expo >= thresholds & not_adopted)
    } else {
      p <- stats::plogis(pars$beta0 + pars$beta_expo * expo)
      which(stats::runif(length(p)) < p & not_adopted)
    }
  }

  set.seed(2026)
  suppressMessages(
    dn <- rdiffnet(n = 30, t = 5, seed.graph = "small-world",
                   seed.p.adopt = list(0.10, 0.10),
                   stop.no.diff = FALSE,
                   adoption_mechanism = per_behavior_mech,
                   adoption_pars      = list(beta0 = -2, beta_expo = 4))
  )

  expect_s3_class(dn, "diffnet")
  expect_setequal(unique(seen_behaviors), c(1L, 2L))
  expect_true(length(seen_behaviors) > 0L)
})

# ----------------------------------------------------------------------------
# Cross-state visibility: -expo_all- carries the n x 1 x Q slice
# ----------------------------------------------------------------------------

test_that("adoption_mechanism with -expo_all- formal receives n x 1 x Q array", {
  observed_dim <- NULL

  cross_mech <- function(expo, thresholds, not_adopted, time, pars,
                         behavior, expo_all) {
    if (is.null(observed_dim)) observed_dim <<- dim(expo_all)
    # Standard threshold logic per behaviour
    which(expo >= thresholds & not_adopted)
  }

  set.seed(2026)
  suppressMessages(
    dn <- rdiffnet(n = 30, t = 5, seed.graph = "small-world",
                   seed.p.adopt = list(0.10, 0.10),
                   stop.no.diff = FALSE,
                   adoption_mechanism = cross_mech)
  )

  expect_s3_class(dn, "diffnet")
  expect_equal(length(observed_dim), 3L)
  expect_equal(observed_dim[1], 30L)             # n
  expect_equal(observed_dim[2], 1L)              # current slice only
  expect_equal(observed_dim[3], 2L)              # Q behaviours
})

# ----------------------------------------------------------------------------
# Toy B (the central paper demo): coupled disease + mask simulation
# ----------------------------------------------------------------------------

test_that("Toy B - SIR + mask cross-coupled via -behavior- and -expo_all-", {
  # behavior 1 = disease (SIR-style transmission, masking reduces it)
  # behavior 2 = mask (logit on social exposure to maskers + local prevalence)
  both_mechanisms <- function(expo, thresholds, not_adopted, time, pars,
                              behavior, expo_all) {
    if (behavior == 1L) {
      # Disease: transmission rate scaled by neighbours masked
      masked_now        <- expo_all[, 1L, 2L]
      protective_factor <- 1 - pars$disease$mask_efficacy * masked_now
      p <- pmax(pmin(pars$disease$transmission_rate *
                       expo * protective_factor, 1), 0)
      which((stats::runif(length(p)) < p) & not_adopted)
    } else {
      # Mask: logit on mask exposure + local disease prevalence (proxy = expo_disease)
      expo_mask    <- expo
      expo_disease <- expo_all[, 1L, 1L]
      p <- stats::plogis(pars$mask$beta0 +
                         pars$mask$beta_expo    * expo_mask +
                         pars$mask$beta_disease * expo_disease)
      which((stats::runif(length(p)) < p) & not_adopted)
    }
  }

  set.seed(2026)
  suppressMessages(
    dn <- rdiffnet(n = 50, t = 8, seed.graph = "small-world",
                   seed.p.adopt = list(0.10, 0.10),
                   stop.no.diff = FALSE,
                   adoption_mechanism = both_mechanisms,
                   adoption_pars = list(
                     disease = list(transmission_rate = 0.30,
                                    mask_efficacy     = 0.50),
                     mask    = list(beta0       = -2,
                                    beta_expo   = 4,
                                    beta_disease = 3)
                   ),
                   source_attribution = list(source_attribution_uniform, NULL))
  )

  expect_true(is.diffnet_epi(dn))
  expect_equal(dim(dn$toa), c(50L, 2L))

  tr <- transmission_tree(dn)
  # Only disease (virus_id = 1) is lineage-tracked
  expect_true(all(tr$virus_id == 1L))
  # Should contain at least the seed rows for behaviour 1
  expect_true(nrow(tr) >= 1L)
})
