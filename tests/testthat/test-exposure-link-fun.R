context("Exposure link / kernel function")
library(netdiffuseR)

# ---- Helpers --------------------------------------------------------------
# Build a small valued dynamic graph + cumadopt that we can reuse.
mk_fixture <- function(seed = 71) {
  set.seed(seed)
  n <- 8L
  g <- rgraph_er(n, t = 1, p = 0.5)
  g@x <- runif(length(g@x), min = 0.1, max = 2.0)
  cumadopt <- matrix(0, nrow = n, ncol = 1)
  cumadopt[1:3, 1] <- 1
  list(g = list(g), cumadopt = cumadopt, n = n)
}

# ---- Identity (default) ---------------------------------------------------
test_that("link_fun = 'identity' reproduces the default exposure", {
  fx <- mk_fixture()
  e_default <- exposure(fx$g, fx$cumadopt, valued = TRUE)
  e_ident   <- exposure(fx$g, fx$cumadopt, valued = TRUE,
                        link_fun = "identity")
  expect_equal(e_default, e_ident)

  # NULL is treated as identity as well
  e_null <- exposure(fx$g, fx$cumadopt, valued = TRUE, link_fun = NULL)
  expect_equal(e_default, e_null)
})

# ---- Linear ---------------------------------------------------------------
test_that("link_fun = 'linear' applies min(beta * w, 1) to edge weights", {
  fx <- mk_fixture()

  g_lin <- fx$g
  g_lin[[1]]@x <- pmin(0.4 * g_lin[[1]]@x, 1)
  e_manual <- exposure(g_lin, fx$cumadopt, valued = TRUE)

  e_kernel <- exposure(fx$g, fx$cumadopt, valued = TRUE,
                       link_fun = "linear", link_pars = list(beta = 0.4))
  expect_equal(e_kernel, e_manual)

  # Linear with beta large enough to saturate at 1
  g_sat <- fx$g
  g_sat[[1]]@x <- pmin(1000 * g_sat[[1]]@x, 1)
  e_sat_manual <- exposure(g_sat, fx$cumadopt, valued = TRUE)
  e_sat_kernel <- exposure(fx$g, fx$cumadopt, valued = TRUE,
                           link_fun = "linear",
                           link_pars = list(beta = 1000))
  expect_equal(e_sat_kernel, e_sat_manual)
})

# ---- Sigmoid --------------------------------------------------------------
test_that("link_fun = 'sigmoid' applies plogis((w - h)/scale)", {
  fx <- mk_fixture()

  g_sig <- fx$g
  g_sig[[1]]@x <- stats::plogis((g_sig[[1]]@x - 1.0) / 0.5)
  e_manual <- exposure(g_sig, fx$cumadopt, valued = TRUE)

  e_kernel <- exposure(fx$g, fx$cumadopt, valued = TRUE,
                       link_fun = "sigmoid",
                       link_pars = list(h = 1.0, scale = 0.5))
  expect_equal(e_kernel, e_manual)
})

# ---- Wells-Riley ----------------------------------------------------------
test_that("link_fun = 'wells-riley' applies 1 - exp(-beta * w)", {
  fx <- mk_fixture()

  g_wr <- fx$g
  g_wr[[1]]@x <- 1 - exp(-0.7 * g_wr[[1]]@x)
  e_manual <- exposure(g_wr, fx$cumadopt, valued = TRUE)

  e_kernel <- exposure(fx$g, fx$cumadopt, valued = TRUE,
                       link_fun = "wells-riley",
                       link_pars = list(beta = 0.7))
  expect_equal(e_kernel, e_manual)
})

# ---- User-supplied function ----------------------------------------------
test_that("link_fun accepts a single-argument user function with baked pars", {
  fx <- mk_fixture()
  # A fully-specified closure -- no link_pars needed.
  random_func <- function(x) 1 - exp(-0.3 * x)

  g_ref <- fx$g
  g_ref[[1]]@x <- 1 - exp(-0.3 * g_ref[[1]]@x)
  e_manual <- exposure(g_ref, fx$cumadopt, valued = TRUE)

  e_user <- exposure(fx$g, fx$cumadopt, valued = TRUE,
                     link_fun = random_func)
  expect_equal(e_user, e_manual)
})

test_that("user-supplied link_fun must return same-length output", {
  fx <- mk_fixture()
  bad_kernel <- function(w) w[-1L]
  expect_error(
    exposure(fx$g, fx$cumadopt, valued = TRUE, link_fun = bad_kernel),
    "same length"
  )
})

# ---- Parameter validation -------------------------------------------------
test_that("Missing required link_pars raises informative errors", {
  fx <- mk_fixture()

  expect_error(
    exposure(fx$g, fx$cumadopt, valued = TRUE,
             link_fun = "linear", link_pars = list()),
    "beta"
  )
  expect_error(
    exposure(fx$g, fx$cumadopt, valued = TRUE,
             link_fun = "wells-riley", link_pars = list()),
    "beta"
  )
  expect_error(
    exposure(fx$g, fx$cumadopt, valued = TRUE,
             link_fun = "sigmoid", link_pars = list(h = 1)),
    "scale"
  )
  expect_error(
    exposure(fx$g, fx$cumadopt, valued = TRUE,
             link_fun = "sigmoid", link_pars = list(scale = 1)),
    "h"
  )
})

test_that("Unknown link_fun name errors informatively", {
  fx <- mk_fixture()
  expect_error(
    exposure(fx$g, fx$cumadopt, valued = TRUE,
             link_fun = "gompertz", link_pars = list()),
    "Unknown link_fun"
  )
})

# ---- Interaction with `valued` -------------------------------------------
test_that("Non-identity link_fun warns and forces valued = TRUE", {
  fx <- mk_fixture()
  expect_warning(
    e_forced <- exposure(fx$g, fx$cumadopt, valued = FALSE,
                         link_fun = "wells-riley",
                         link_pars = list(beta = 0.5)),
    "valued"
  )
  e_valued <- exposure(fx$g, fx$cumadopt, valued = TRUE,
                       link_fun = "wells-riley",
                       link_pars = list(beta = 0.5))
  expect_equal(e_forced, e_valued)
})

# ---- diffnet dispatch threads link_fun through ----------------------------
test_that("diffnet dispatch accepts link_fun/link_pars", {
  set.seed(9)
  dn <- suppressWarnings(rdiffnet(n = 20, t = 4, seed.graph = "small-world",
                                  seed.p.adopt = 0.1,
                                  stop.no.diff = FALSE))
  # Force a non-trivial valued graph
  for (i in seq_along(dn$graph)) {
    dn$graph[[i]]@x <- runif(length(dn$graph[[i]]@x), 0.1, 2.0)
  }

  e_ident <- exposure(dn, valued = TRUE)
  e_wr    <- exposure(dn, valued = TRUE,
                      link_fun = "wells-riley",
                      link_pars = list(beta = 1.2))

  # Output must have the same shape as identity case
  expect_equal(dim(e_ident), dim(e_wr))

  # Wells-Riley output is bounded in [0, 1]
  expect_true(all(e_wr >= 0))
  expect_true(all(e_wr <= 1 + 1e-10))
})
