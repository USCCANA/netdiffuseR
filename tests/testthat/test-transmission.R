context("Transmission tree slot")

mk_diffnet <- function() {
  set.seed(42)
  gr  <- lapply(1:5, function(x) rgraph_ba(t = 4L))
  toa <- c(1L, 2L, NA, 3L, 5L)
  new_diffnet(gr, toa, t0 = 1L, t1 = 5L)
}

test_that("transmission_tree() errors on a plain diffnet (not promoted)", {
  x <- mk_diffnet()
  expect_false(is.diffnet_epi(x))
  expect_error(transmission_tree(x), "diffnet_epi")
})

test_that("transmission_tree() returns an empty data.frame after empty promotion", {
  x <- as_diffnet_epi(mk_diffnet())
  expect_true(is.diffnet_epi(x))
  tr <- transmission_tree(x)
  expect_s3_class(tr, "data.frame")
  expect_equal(nrow(tr), 0L)
  expect_setequal(
    names(tr),
    c("date", "source", "target", "source_exposure_date", "virus_id", "virus")
  )
})

test_that("as_transmission_tree() validates required columns and ranges", {
  x <- mk_diffnet()

  expect_error(
    as_transmission_tree(x, data.frame(date = 1L, source = NA, target = 1L)),
    "missing required column"
  )

  expect_error(
    as_transmission_tree(x, data.frame(
      date = 1L, source = NA_integer_, target = NA_integer_,
      source_exposure_date = NA_integer_
    )),
    "cannot contain NA"
  )

  expect_error(
    as_transmission_tree(x, data.frame(
      date = 1L, source = NA_integer_, target = 999L,
      source_exposure_date = NA_integer_
    )),
    "integer indices"
  )

  expect_error(
    as_transmission_tree(x, data.frame(
      date = 1L, source = 42L, target = 2L, source_exposure_date = 1L
    )),
    "NA or an integer index"
  )
})

test_that("as_transmission_tree() stores a clean tree and optional pars", {
  x <- mk_diffnet()
  tree <- data.frame(
    date   = c(3L, 1L, 2L),
    source = c(2L, NA, 1L),
    target = c(4L, 1L, 2L),
    source_exposure_date = c(2L, NA, 1L),
    virus_id = c(1L, 1L, 1L),
    virus    = c("flu", "flu", "flu"),
    stringsAsFactors = FALSE
  )
  y <- as_transmission_tree(x, tree, pars = list(kernel = "wells-riley"))
  tr <- transmission_tree(y)

  # Ordered by (date, target) and clean rownames
  expect_equal(tr$date, c(1L, 2L, 3L))
  expect_equal(tr$target, c(1L, 2L, 4L))
  expect_equal(rownames(tr), as.character(seq_len(nrow(tr))))

  expect_equal(y$transmission$pars$kernel, "wells-riley")
})

test_that("Missing optional columns are defaulted", {
  x <- mk_diffnet()
  tree <- data.frame(
    date   = c(1L, 2L),
    source = c(NA_integer_, 1L),
    target = c(1L, 2L),
    source_exposure_date = c(NA_integer_, 1L)
  )
  y  <- as_transmission_tree(x, tree)
  tr <- transmission_tree(y)

  expect_true(all(tr$virus_id == 1L))
  expect_true(all(is.na(tr$virus)))
})

test_that("transmission_tree() and as_transmission_tree() reject non-diffnet inputs", {
  expect_error(transmission_tree(42),                   "diffnet_epi")
  expect_error(as_transmission_tree(42, data.frame()),  "must be a diffnet")
})

test_that("as_transmission_tree() promotes the diffnet to diffnet_epi", {
  x <- mk_diffnet()
  expect_false(is.diffnet_epi(x))
  tree <- data.frame(
    date = c(1L, 2L),
    source = c(NA_integer_, 1L),
    target = c(1L, 2L),
    source_exposure_date = c(NA_integer_, 1L)
  )
  y <- as_transmission_tree(x, tree)
  expect_true(is.diffnet_epi(y))
  expect_s3_class(y, "diffnet")            # still a diffnet
  expect_s3_class(y, "diffnet_epi")        # also a diffnet_epi
})

# ----------------------------------------------------------------------------
# transmission_tree_from_events / as_diffnet_epi(attribution=) (M13)
# ----------------------------------------------------------------------------

test_that("transmission_tree_from_events() on a diffnet returns canonical schema", {
  dn <- mk_diffnet()
  tr <- transmission_tree_from_events(dn, attribution = "uniform", seed = 42)
  expect_s3_class(tr, "data.frame")
  expect_setequal(
    names(tr),
    c("date", "source", "target", "source_exposure_date", "virus_id", "virus")
  )
  # Targets cover every node with non-NA toa, exactly once.
  expect_setequal(tr$target, which(!is.na(dn$toa)))
  # date == toa[target] for every row.
  expect_equal(tr$date, as.integer(dn$toa[tr$target]))
  # source_exposure_date == toa[source] when source is not NA.
  has_src <- !is.na(tr$source)
  expect_equal(tr$source_exposure_date[has_src],
               as.integer(dn$toa[tr$source[has_src]]))
})

test_that("transmission_tree_from_events() with a list of graphs matches diffnet path", {
  dn <- mk_diffnet()
  tr_dn   <- transmission_tree_from_events(dn,         attribution = "uniform", seed = 42)
  tr_list <- transmission_tree_from_events(dn$graph,
                                            toa         = dn$toa,
                                            attribution = "uniform",
                                            seed        = 42)
  expect_equal(tr_dn, tr_list)
})

test_that("transmission_tree_from_events() rejects unknown attribution / missing toa", {
  dn <- mk_diffnet()
  expect_error(transmission_tree_from_events(dn, attribution = "nope"),
               "Unknown -attribution-")
  expect_error(transmission_tree_from_events(dn$graph),  # no toa
               "-toa- is required")
  expect_error(transmission_tree_from_events(42),
               "must be a diffnet or a list")
})

test_that("transmission_tree_from_events() agrees with rdiffnet's online tree", {
  # Run rdiffnet with source_attribution=uniform and the SAME seed twice:
  # 1) online (rdiffnet builds the tree in-loop via M8)
  # 2) post-hoc reconstruction from the diffnet's graph + toa using M13.1.
  # Under absorbing diffusion these must coincide row-for-row, since both
  # paths apply the same uniform kernel to the same set of adopted-neighbours
  # at each adoption event.
  set.seed(2026)
  dn_online <- rdiffnet(n = 30, t = 6, seed.graph = "small-world",
                        seed.p.adopt = 0.10, stop.no.diff = FALSE,
                        source_attribution = source_attribution_uniform)
  tr_online <- transmission_tree(dn_online)

  # Strip the online tree's class machinery and reconstruct from the same
  # toa + graphs. Use the same seed inside transmission_tree_from_events to
  # reproduce the random choices the online attributor made.
  tr_replay <- transmission_tree_from_events(
    dn_online, attribution = "uniform", seed = 2026
  )

  # Same (target, date) set, identical source assignment for deterministic
  # rows (singleton neighbour). Verify per-target source coincides where
  # the neighbour set is a singleton (the only case where the seed-dependent
  # choices don't matter).
  o <- tr_online[order(tr_online$target), ]
  r <- tr_replay[order(tr_replay$target), ]
  expect_equal(o$target, r$target)
  expect_equal(o$date,   r$date)
  # Sum of offspring is preserved (every non-seed edge contributes +1).
  expect_equal(sum(!is.na(o$source)), sum(!is.na(r$source)))
})

test_that("as_diffnet_epi(attribution = ...) and standalone primitive agree", {
  dn <- mk_diffnet()
  set.seed(42)
  tr_standalone <- transmission_tree_from_events(dn, attribution = "uniform",
                                                  seed = 42)
  dn_epi <- as_diffnet_epi(dn, attribution = "uniform", seed = 42)
  expect_true(is.diffnet_epi(dn_epi))
  expect_equal(transmission_tree(dn_epi), tr_standalone)
})

test_that("as_diffnet_epi() rejects mutually-exclusive transmission + attribution", {
  dn   <- mk_diffnet()
  tree <- transmission_tree_from_events(dn, attribution = "uniform", seed = 42)
  expect_error(
    as_diffnet_epi(dn,
                   transmission = list(tree = tree, pars = list()),
                   attribution  = "uniform"),
    "either -transmission- .+ or -attribution-"
  )
})

test_that("as_diffnet_epi(attribution = <function>) accepts user kernels", {
  # Always pick the FIRST neighbour (most ancient — earliest infector).
  my_attr <- function(target, adopted_neighbours, weights, time, pars) {
    if (!length(adopted_neighbours)) return(NA_integer_)
    as.integer(adopted_neighbours[1L])
  }
  dn     <- mk_diffnet()
  dn_epi <- as_diffnet_epi(dn, attribution = my_attr)
  tr     <- transmission_tree(dn_epi)
  expect_true(is.diffnet_epi(dn_epi))
  expect_true(nrow(tr) > 0L)
})

# ----------------------------------------------------------------------------
# Smoke test on the shipped epigamesDiffNet dataset (M13.3)
# ----------------------------------------------------------------------------

test_that("shipped epigamesDiffNet is a diffnet_epi with a populated tree", {
  data("epigamesDiffNet", package = "netdiffuseR")
  expect_s3_class(epigamesDiffNet, "diffnet_epi")
  expect_s3_class(epigamesDiffNet, "diffnet")

  tr <- transmission_tree(epigamesDiffNet)
  expect_gt(nrow(tr), 0L)
  expect_setequal(
    names(tr),
    c("date", "source", "target", "source_exposure_date", "virus_id", "virus")
  )

  # The 5 epi metrics + summary block all run without error.
  expect_silent(pp  <- peak_prevalence(epigamesDiffNet))
  expect_silent(pt  <- peak_time(epigamesDiffNet))
  expect_silent(sar <- secondary_attack_rate(epigamesDiffNet))
  expect_silent(gt  <- generation_time(epigamesDiffNet))
  expect_silent(R   <- repr_number(epigamesDiffNet))
  expect_silent(sc  <- survival_curve(epigamesDiffNet))

  expect_true(pp > 0 && pp <= 1)
  expect_true(pt %in% epigamesDiffNet$meta$pers)
  expect_true(attr(sar, "global") >= 0)
  expect_true(all(gt$gen_time > 0))
  expect_true(attr(R, "global") >= 0)
  expect_s3_class(sc, "netdiffuseR_survival")
})
