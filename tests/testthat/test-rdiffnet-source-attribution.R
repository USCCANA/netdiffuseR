context("rdiffnet source_attribution callback")
library(netdiffuseR)

# ----------------------------------------------------------------------------
# Default-NULL path: no behaviour change vs pre-M8 rdiffnet()
# ----------------------------------------------------------------------------

test_that("source_attribution = NULL keeps rdiffnet returning a plain diffnet", {
  set.seed(2026)
  dn <- rdiffnet(n = 25, t = 6, seed.graph = "small-world",
                 seed.p.adopt = 0.10, stop.no.diff = FALSE)
  expect_false(is.diffnet_epi(dn))
  expect_false("transmission" %in% names(dn))
  expect_s3_class(dn, "diffnet")
})

# ----------------------------------------------------------------------------
# Single-function broadcast: auto-promotion + tree well-formed
# ----------------------------------------------------------------------------

test_that("source_attribution_uniform produces a diffnet_epi with a valid tree", {
  set.seed(2026)
  dn <- rdiffnet(n = 30, t = 6, seed.graph = "small-world",
                 seed.p.adopt = 0.10, stop.no.diff = FALSE,
                 source_attribution = source_attribution_uniform)

  expect_true(is.diffnet_epi(dn))
  tr <- transmission_tree(dn)
  expect_s3_class(tr, "data.frame")
  expect_setequal(
    names(tr),
    c("date", "source", "target", "source_exposure_date", "virus_id", "virus")
  )

  # Seed rows have NA source and NA source_exposure_date
  seed_rows <- tr[is.na(tr$source), , drop = FALSE]
  expect_true(nrow(seed_rows) >= 1L)
  expect_true(all(is.na(seed_rows$source_exposure_date)))

  # Non-seed rows: source != NA, source_exposure_date < date, source has
  # toa equal to source_exposure_date
  non_seed <- tr[!is.na(tr$source), , drop = FALSE]
  if (nrow(non_seed) > 0L) {
    expect_true(all(non_seed$source_exposure_date < non_seed$date))
    expect_true(all(non_seed$source_exposure_date == dn$toa[non_seed$source]))
  }

  # Every target appears at most once per virus_id (single adoption per
  # behaviour in absorbing simulation; without disadopt no re-adoption)
  by_virus <- split(tr, tr$virus_id)
  for (chunk in by_virus) {
    expect_equal(length(unique(chunk$target)), nrow(chunk))
  }
})

# ----------------------------------------------------------------------------
# source_attribution_earliest: deterministic choice (sorted-by-toa first)
# ----------------------------------------------------------------------------

test_that("source_attribution_earliest is deterministic for a given seed/graph", {
  set.seed(2026)
  dn1 <- rdiffnet(n = 25, t = 6, seed.graph = "small-world",
                  seed.p.adopt = 0.10, stop.no.diff = FALSE,
                  source_attribution = source_attribution_earliest)

  set.seed(2026)
  dn2 <- rdiffnet(n = 25, t = 6, seed.graph = "small-world",
                  seed.p.adopt = 0.10, stop.no.diff = FALSE,
                  source_attribution = source_attribution_earliest)

  expect_identical(transmission_tree(dn1), transmission_tree(dn2))
})

# ----------------------------------------------------------------------------
# source_attribution_weighted: fallback to uniform on unweighted graph
# ----------------------------------------------------------------------------

test_that("source_attribution_weighted runs end-to-end on an unweighted graph", {
  set.seed(2026)
  dn <- rdiffnet(n = 30, t = 6, seed.graph = "small-world",
                 seed.p.adopt = 0.10, stop.no.diff = FALSE,
                 source_attribution = source_attribution_weighted)
  expect_true(is.diffnet_epi(dn))
  expect_true(nrow(transmission_tree(dn)) >= 1L)
})

# ----------------------------------------------------------------------------
# Multi-behaviour: per-behaviour list with one tracked, one not
# ----------------------------------------------------------------------------

test_that("Per-behaviour source_attribution list tracks selected behaviours only", {
  set.seed(2026)
  suppressMessages(
    dn <- rdiffnet(n = 40, t = 6, seed.graph = "small-world",
                   seed.p.adopt = list(0.10, 0.10),
                   stop.no.diff = FALSE,
                   source_attribution = list(source_attribution_uniform, NULL))
  )

  expect_true(is.diffnet_epi(dn))
  tr <- transmission_tree(dn)

  # Only virus_id == 1 appears in the tree
  expect_true(all(tr$virus_id == 1L))
})

test_that("Per-behaviour source_attribution list with both populated mixes virus_ids", {
  set.seed(2026)
  suppressMessages(
    dn <- rdiffnet(n = 40, t = 6, seed.graph = "small-world",
                   seed.p.adopt = list(0.10, 0.10),
                   stop.no.diff = FALSE,
                   source_attribution = list(source_attribution_uniform,
                                             source_attribution_earliest))
  )

  expect_true(is.diffnet_epi(dn))
  tr <- transmission_tree(dn)
  expect_setequal(unique(tr$virus_id), c(1L, 2L))
})

# ----------------------------------------------------------------------------
# Validation
# ----------------------------------------------------------------------------

test_that("Bad -source_attribution- inputs error informatively", {
  set.seed(2026)
  expect_error(
    rdiffnet(n = 20, t = 4, seed.graph = "small-world",
             seed.p.adopt = 0.10, stop.no.diff = FALSE,
             source_attribution = 42),
    "must be NULL, a function, or a length-Q list"
  )

  # Wrong-length list
  expect_error(
    suppressMessages(rdiffnet(n = 20, t = 4, seed.graph = "small-world",
             seed.p.adopt = list(0.10, 0.10), stop.no.diff = FALSE,
             source_attribution = list(source_attribution_uniform))),
    "must have length"
  )

  # List element that is neither NULL nor a function
  expect_error(
    suppressMessages(rdiffnet(n = 20, t = 4, seed.graph = "small-world",
             seed.p.adopt = list(0.10, 0.10), stop.no.diff = FALSE,
             source_attribution = list(source_attribution_uniform, "earliest"))),
    "must be NULL or a function"
  )
})

# ----------------------------------------------------------------------------
# Individual kernels (out-of-loop unit tests)
# ----------------------------------------------------------------------------

test_that("source_attribution_uniform returns NA when no adopted neighbours", {
  expect_equal(
    source_attribution_uniform(target = 5L,
                               adopted_neighbours = integer(0),
                               weights = numeric(0),
                               time = 3L, pars = list()),
    NA_integer_
  )
})

test_that("source_attribution_earliest picks the first (sorted by toa) entry", {
  expect_equal(
    source_attribution_earliest(target = 5L,
                                adopted_neighbours = c(7L, 3L, 12L),
                                weights = c(1, 1, 1),
                                time = 4L, pars = list()),
    7L
  )
})

test_that("source_attribution_weighted falls back to uniform on equal weights", {
  set.seed(2026)
  picks <- replicate(200, source_attribution_weighted(
    target = 1L,
    adopted_neighbours = c(2L, 3L),
    weights = c(1, 1),
    time = 3L, pars = list()
  ))
  # With equal weights and uniform fallback, both candidates should appear
  expect_true(all(picks %in% c(2L, 3L)))
  expect_true(any(picks == 2L))
  expect_true(any(picks == 3L))
})

test_that("source_attribution_weighted concentrates mass on the heavier neighbour", {
  set.seed(2026)
  picks <- replicate(400, source_attribution_weighted(
    target = 1L,
    adopted_neighbours = c(2L, 3L),
    weights = c(99, 1),
    time = 3L, pars = list()
  ))
  expect_gt(mean(picks == 2L), 0.85)
})
