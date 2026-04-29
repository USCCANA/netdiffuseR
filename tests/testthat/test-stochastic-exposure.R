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

# ---- Continuous (non-Bernoulli) edge weights -----------------------------
# Fixture: graph with floating-point weights resembling seconds of contact.
mk_seconds_graph <- function(seed = 17, n = 10L) {
  set.seed(seed)
  g <- rgraph_er(n, t = 1, p = 0.6)
  g@x <- runif(length(g@x), min = 30, max = 3600)  # 30s .. 1h
  cumadopt <- matrix(0, nrow = n, ncol = 1)
  cumadopt[1:4, 1] <- 1
  list(g = list(g), cumadopt = cumadopt, n = n)
}

test_that("Each link_fun maps seconds-scale weights into valid probabilities", {
  fx <- mk_seconds_graph()

  # Identity: raw seconds > 1 warns and saturates the sampler; the
  # normalised exposure still ends up in [0, 1].
  expect_warning(
    e_id <- exposure(fx$g, fx$cumadopt, valued = TRUE, mode = "stochastic"),
    "weights in \\[0, 1\\]"
  )
  expect_true(all(e_id >= 0) && all(e_id <= 1 + 1e-10))

  set.seed(101)
  e_wr <- exposure(fx$g, fx$cumadopt, valued = TRUE, mode = "stochastic",
                   link_fun = "wells-riley",
                   link_pars = list(beta = 1 / 1800))
  expect_true(all(e_wr >= 0) && all(e_wr <= 1 + 1e-10))

  set.seed(101)
  e_lin <- exposure(fx$g, fx$cumadopt, valued = TRUE, mode = "stochastic",
                    link_fun = "linear",
                    link_pars = list(beta = 1 / 5000))
  expect_true(all(e_lin >= 0) && all(e_lin <= 1 + 1e-10))

  set.seed(101)
  e_sig <- exposure(fx$g, fx$cumadopt, valued = TRUE, mode = "stochastic",
                    link_fun = "sigmoid",
                    link_pars = list(h = 600, scale = 300))
  expect_true(all(e_sig >= 0) && all(e_sig <= 1 + 1e-10))
})

test_that("Out-of-range warning fires only when raw weights leave [0, 1]", {
  fx <- mk_seconds_graph()

  # Kernel maps into [0, 1]: no warning.
  expect_warning(
    exposure(fx$g, fx$cumadopt, valued = TRUE, mode = "stochastic",
             link_fun = "wells-riley",
             link_pars = list(beta = 1 / 1800)),
    regexp = NA
  )

  # Caller pre-normalises: no warning.
  set.seed(99); n <- fx$n
  g_ok <- rgraph_er(n, t = 1, p = 0.6); g_ok@x <- runif(length(g_ok@x))
  expect_warning(
    exposure(list(g_ok), fx$cumadopt, valued = TRUE, mode = "stochastic"),
    regexp = NA
  )
})

test_that("Stochastic denominator excludes zero-weight stored neighbours", {
  # Node 1 has 4 stored edges to adopters; 3 weigh 0, 1 weighs 1.
  # Old behaviour: exposure[1] = 1/4; corrected: 1/1 = 1.
  n <- 5L
  A <- matrix(0, n, n); A[1, 2:5] <- 1
  G <- methods::as(A, "CsparseMatrix")
  stopifnot(length(G@x) == 4)
  G@x <- c(0, 0, 0, 1)  # three stored zeros, one stored one

  cumadopt <- matrix(0, n, 1); cumadopt[2:5, 1] <- 1

  set.seed(123)
  e <- exposure(list(G), cumadopt, valued = TRUE, mode = "stochastic",
                self = TRUE)

  expect_equal(e[1, 1], 1)
})

test_that("Zero-weight self-loops do not blow up the exposure", {
  set.seed(7); n <- 6L
  g <- rgraph_er(n, t = 1, p = 0.6)
  g[1, 1] <- 1; g[1, 1] <- 0  # force a stored 0 on the diagonal
  cumadopt <- matrix(0, n, 1); cumadopt[1, 1] <- 1

  e <- exposure(list(g), cumadopt, valued = TRUE, mode = "stochastic",
                self = FALSE,
                link_fun = "wells-riley",
                link_pars = list(beta = 0.5))
  expect_false(anyNA(e))
  expect_true(all(is.finite(e)))
  expect_true(all(e >= 0) && all(e <= 1 + 1e-10))
})
