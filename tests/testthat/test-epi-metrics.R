context("Epidemiological metrics for diffnet / diffnet_epi (M10)")
library(netdiffuseR)

mk_absorbing_dn <- function() {
  set.seed(2026)
  rdiffnet(n = 30, t = 6, seed.graph = "small-world",
           seed.p.adopt = 0.10, stop.no.diff = FALSE)
}

mk_recovery_dn <- function() {
  # Hand-built diffnet with explicit recovery via -status-, so toa is preserved
  # for all adopters (rdiffnet's -disadopt- resets toa to NA on disadoption,
  # which is not what we want for the KM curve test).
  set.seed(2026)
  g <- lapply(1:6, function(t) rgraph_ba(t = 5L))
  status <- rbind(
    c(1L, 1L, 1L, 0L, 0L, 0L),    # node 1: adopt t=1, recover t=4
    c(0L, 1L, 1L, 1L, 0L, 0L),    # node 2: adopt t=2, recover t=5
    c(0L, 0L, 0L, 0L, 0L, 0L),    # node 3: never adopts
    c(0L, 0L, 1L, 1L, 1L, 1L),    # node 4: adopt t=3, absorbing (censored)
    c(0L, 0L, 0L, 1L, 1L, 0L),    # node 5: adopt t=4, recover t=6
    c(0L, 1L, 1L, 1L, 1L, 1L)     # node 6: adopt t=2, absorbing
  )
  new_diffnet(g, status = status, t0 = 1L, t1 = 6L)
}

mk_epi_dn <- function() {
  set.seed(2026)
  rdiffnet(n = 30, t = 6, seed.graph = "small-world",
           seed.p.adopt = 0.10, stop.no.diff = FALSE,
           source_attribution = source_attribution_uniform)
}

# ----------------------------------------------------------------------------
# peak_prevalence / peak_time
# ----------------------------------------------------------------------------

test_that("peak_prevalence and peak_time work on plain diffnet", {
  dn <- mk_absorbing_dn()
  pp <- peak_prevalence(dn)
  pt <- peak_time(dn)
  expect_type(pp, "double")
  expect_length(pp, 1L)
  expect_true(pp >= 0 && pp <= 1)
  expect_type(pt, "integer")
  expect_length(pt, 1L)
  expect_true(pt %in% dn$meta$pers)
})

test_that("peak_prevalence is non-decreasing for an absorbing single-virus sim", {
  dn <- mk_absorbing_dn()
  prev_per_t <- colSums(dn$status) / nrow(dn$status)
  expect_identical(peak_prevalence(dn), max(prev_per_t))
  # In absorbing diffusion the peak is at the final period
  expect_equal(peak_time(dn), dn$meta$pers[which.max(prev_per_t)])
})

test_that("peak_prevalence returns named vector for multi-behaviour", {
  set.seed(2026)
  suppressMessages(
    dn <- rdiffnet(n = 30, t = 6, seed.graph = "small-world",
                   seed.p.adopt = list(0.10, 0.10), stop.no.diff = FALSE)
  )
  pp <- peak_prevalence(dn)
  pt <- peak_time(dn)
  expect_length(pp, 2L)
  expect_length(pt, 2L)
  expect_true(!is.null(names(pp)))
})

test_that("peak_prevalence on diffnet_epi delegates correctly", {
  dn_epi <- mk_epi_dn()
  expect_equal(peak_prevalence(dn_epi),
               peak_prevalence(structure(dn_epi, class = "diffnet")))
})

# ----------------------------------------------------------------------------
# survival_curve
# ----------------------------------------------------------------------------

test_that("survival_curve returns a netdiffuseR_survival data.frame", {
  dn <- mk_absorbing_dn()
  s  <- survival_curve(dn)
  expect_s3_class(s, "netdiffuseR_survival")
  expect_s3_class(s, "data.frame")
  expect_setequal(names(s), c("time", "n_at_risk", "n_recovered", "survival"))
})

test_that("survival_curve for absorbing diffusion has survival = 1 throughout", {
  dn <- mk_absorbing_dn()
  s  <- survival_curve(dn)
  expect_true(all(s$survival == 1))
  expect_true(all(s$n_recovered == 0))
})

test_that("survival_curve drops below 1 for a recovery simulation", {
  dn <- mk_recovery_dn()
  s  <- survival_curve(dn)
  # Some recoveries should be present
  expect_true(sum(s$n_recovered) >= 1L)
  expect_true(any(s$survival < 1))
})

test_that("survival_curve print method runs silently", {
  dn  <- mk_absorbing_dn()
  out <- capture.output(print(survival_curve(dn)))
  expect_true(any(grepl("Survival curve", out)))
  expect_true(any(grepl("as.data.frame", out)))
})

# ----------------------------------------------------------------------------
# secondary_attack_rate
# ----------------------------------------------------------------------------

test_that("secondary_attack_rate errors on plain diffnet", {
  dn <- mk_absorbing_dn()
  expect_error(secondary_attack_rate(dn), "diffnet_epi")
})

test_that("secondary_attack_rate returns a netdiffuseR_sar data.frame", {
  dn  <- mk_epi_dn()
  sar <- secondary_attack_rate(dn)
  expect_s3_class(sar, "netdiffuseR_sar")
  expect_s3_class(sar, "data.frame")
  expect_setequal(
    names(sar),
    c("source", "virus_id", "n_secondary", "n_contacts", "sar")
  )
  # All per-source rates are in [0, 1] (or NA when contacts = 0)
  sar_vals <- sar$sar[!is.na(sar$sar)]
  expect_true(all(sar_vals >= 0))
  expect_true(all(sar_vals <= 1))
})

test_that("secondary_attack_rate aggregate matches sum-over-sum formula", {
  dn  <- mk_epi_dn()
  sar <- secondary_attack_rate(dn)
  expected_global <- sum(sar$n_secondary) / sum(sar$n_contacts)
  expect_equal(attr(sar, "global"), expected_global)
})

test_that("secondary_attack_rate print mentions Aggregate and infectors", {
  dn  <- mk_epi_dn()
  out <- capture.output(print(secondary_attack_rate(dn)))
  expect_true(any(grepl("Secondary Attack Rate", out)))
  expect_true(any(grepl("Aggregate", out)))
  expect_true(any(grepl("infector", out)))
})

# ----------------------------------------------------------------------------
# generation_time
# ----------------------------------------------------------------------------

test_that("generation_time errors on plain diffnet", {
  dn <- mk_absorbing_dn()
  expect_error(generation_time(dn), "diffnet_epi")
})

test_that("generation_time returns a netdiffuseR_generation_time data.frame", {
  dn <- mk_epi_dn()
  gt <- generation_time(dn)
  expect_s3_class(gt, "netdiffuseR_generation_time")
  expect_s3_class(gt, "data.frame")
  expect_true("gen_time" %in% names(gt))
  # All gen_times are positive (date - source_exposure_date > 0)
  expect_true(all(gt$gen_time > 0))
})

test_that("generation_time gen_time equals date - source_exposure_date", {
  dn <- mk_epi_dn()
  gt <- generation_time(dn)
  expect_equal(as.integer(gt$gen_time),
               as.integer(gt$date - gt$source_exposure_date))
})

test_that("generation_time print mentions Mean and Median", {
  dn  <- mk_epi_dn()
  out <- capture.output(print(generation_time(dn)))
  expect_true(any(grepl("Generation time", out)))
  expect_true(any(grepl("Mean", out)))
  expect_true(any(grepl("Median", out)))
})

# ----------------------------------------------------------------------------
# repr_number (M12)
# ----------------------------------------------------------------------------

test_that("repr_number errors on plain diffnet (and on arbitrary input)", {
  dn <- mk_absorbing_dn()
  expect_error(repr_number(dn), "diffnet_epi")
  expect_error(repr_number(list(a = 1)), "diffnet_epi")
})

test_that("repr_number returns a netdiffuseR_repr data.frame", {
  dn <- mk_epi_dn()
  R  <- repr_number(dn)
  expect_s3_class(R, "netdiffuseR_repr")
  expect_s3_class(R, "data.frame")
  expect_setequal(names(R), c("node", "virus_id", "n_offspring"))
  expect_type(R$n_offspring, "integer")
  expect_true(all(R$n_offspring >= 0L))
})

test_that("repr_number global = mean(n_offspring) and matches non-seed source count", {
  dn <- mk_epi_dn()
  R  <- repr_number(dn)
  tr <- transmission_tree(dn)

  # Aggregate is mean offspring
  expect_equal(attr(R, "global"), mean(R$n_offspring))

  # Sum of offspring equals number of non-seed source rows in the tree
  expect_equal(sum(R$n_offspring), sum(!is.na(tr$source)))

  # Row count equals number of unique (target, virus_id) cases
  expect_equal(
    nrow(R),
    nrow(unique(tr[, c("target", "virus_id"), drop = FALSE]))
  )
})

test_that("repr_number print mentions Reproduction number and Mean offspring", {
  dn  <- mk_epi_dn()
  out <- capture.output(print(repr_number(dn)))
  expect_true(any(grepl("Reproduction number", out)))
  expect_true(any(grepl("Mean offspring", out)))
  expect_true(any(grepl("as.data.frame", out)))
  # Single-virus path: no aggregate-banner line
  expect_false(any(grepl("Aggregate over", out)))
})

test_that("repr_number print flags aggregate when tree has multiple viruses", {
  set.seed(2026)
  suppressMessages(
    dn <- rdiffnet(n = 40, t = 6, seed.graph = "small-world",
                   seed.p.adopt = list(0.10, 0.05), stop.no.diff = FALSE,
                   source_attribution = source_attribution_uniform)
  )
  R   <- repr_number(dn)
  out <- capture.output(print(R))
  expect_true(any(grepl("Aggregate over 2 behaviours", out)))
  expect_true(any(grepl("Per-virus R", out)))
  # Both virus_ids must show up in the per-virus rollup lines
  expect_true(any(grepl("virus 1:", out)))
  expect_true(any(grepl("virus 2:", out)))
  # Aggregate R must equal mean over all rows (pooled)
  expect_equal(attr(R, "global"), mean(R$n_offspring))
})

# ----------------------------------------------------------------------------
# summary.diffnet_epi
# ----------------------------------------------------------------------------

test_that("summary.diffnet_epi extends summary.diffnet with epi block", {
  dn  <- mk_epi_dn()
  out <- capture.output(summary(dn))
  # Base diffnet summary fields
  expect_true(any(grepl("Diffusion network summary statistics", out)))
  # New epi block fields
  expect_true(any(grepl("Epidemiological metrics", out)))
  expect_true(any(grepl("Peak prevalence", out)))
  expect_true(any(grepl("Secondary Attack Rate", out)))
})

test_that("summary.diffnet_epi returns invisibly (no print on assignment)", {
  dn  <- mk_epi_dn()
  out <- capture.output(invisible(summary(dn)))
  expect_true(length(out) > 0L)
})
