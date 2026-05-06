context("Canonical -status- slot and accessors")

mk_graph <- function(n = 5L, T = 5L, seed = 1L) {
  set.seed(seed)
  lapply(seq_len(T), function(x) rgraph_ba(t = n - 1L))
}

# ----------------------------------------------------------------------------
# Constructor: -toa- only path (legacy, must be bit-identical)
# ----------------------------------------------------------------------------

test_that("toa-only path keeps legacy behaviour and gains a $status slot", {
  gr  <- mk_graph()
  toa <- c(1L, 2L, NA, 3L, 5L)
  x   <- new_diffnet(gr, toa, t0 = 1L, t1 = 5L)

  expect_true(inherits(x, "diffnet"))
  expect_null(x$transmission)
  expect_null(x$tod)                                      # legacy slot retired
  expect_false(is.null(x$status))                         # status slot present
  expect_identical(x$status, x$cumadopt)                  # alias of $cumadopt
})

# ----------------------------------------------------------------------------
# Constructor: -status- only path (new)
# ----------------------------------------------------------------------------

test_that("status-only path derives toa and rebuilds cumadopt", {
  gr <- mk_graph()
  # Node 1: adopts t=1, recovers t=3.   -> 1,1,0,0,0
  # Node 2: adopts t=2, recovers t=4.   -> 0,1,1,0,0
  # Node 3: never adopts.                -> 0,0,0,0,0
  # Node 4: adopts t=3, absorbing.       -> 0,0,1,1,1
  # Node 5: adopts t=5, absorbing.       -> 0,0,0,0,1
  status <- rbind(
    c(1L, 1L, 0L, 0L, 0L),
    c(0L, 1L, 1L, 0L, 0L),
    c(0L, 0L, 0L, 0L, 0L),
    c(0L, 0L, 1L, 1L, 1L),
    c(0L, 0L, 0L, 0L, 1L)
  )
  x <- new_diffnet(gr, status = status, t0 = 1L, t1 = 5L)

  expect_equal(as.integer(x$toa), c(1L, 2L, NA, 3L, 5L))
  expect_identical(x$status, x$cumadopt)
  expect_equal(as.integer(x$cumadopt[1, ]), c(1L, 1L, 0L, 0L, 0L))
  expect_equal(as.integer(x$cumadopt[2, ]), c(0L, 1L, 1L, 0L, 0L))
  expect_equal(as.integer(x$cumadopt[4, ]), c(0L, 0L, 1L, 1L, 1L))
})

test_that("status-only path captures multi-cycle (re-adoption)", {
  gr <- mk_graph(n = 3L, T = 6L, seed = 7L)
  # Node 1: adopt t=1, recover t=3, re-adopt t=5    -> 1,1,0,0,1,1
  # Node 2: never adopts                             -> 0,0,0,0,0,0
  # Node 3: adopt t=2, never recovers (absorbing)    -> 0,1,1,1,1,1
  status <- rbind(
    c(1L, 1L, 0L, 0L, 1L, 1L),
    c(0L, 0L, 0L, 0L, 0L, 0L),
    c(0L, 1L, 1L, 1L, 1L, 1L)
  )
  x <- new_diffnet(gr, status = status, t0 = 1L, t1 = 6L)

  # toa() returns the FIRST adoption time per node
  expect_equal(as.integer(toa(x)), c(1L, NA, 2L))

  # adopt is a "fresh adoption" indicator -> two 1s for node 1 (t=1 and t=5)
  expect_equal(as.integer(rowSums(x$adopt)), c(2L, 0L, 1L))
  expect_equal(as.integer(x$adopt[1, ]), c(1L, 0L, 0L, 0L, 1L, 0L))
})

# ----------------------------------------------------------------------------
# Constructor: -toa- and -status- both passed (warn and prefer status)
# ----------------------------------------------------------------------------

test_that("supplying both toa and status warns and uses status", {
  gr <- mk_graph()
  status <- rbind(
    c(1L, 1L, 0L, 0L, 0L),
    c(0L, 1L, 1L, 0L, 0L),
    c(0L, 0L, 0L, 0L, 0L),
    c(0L, 0L, 1L, 1L, 1L),
    c(0L, 0L, 0L, 0L, 1L)
  )
  toa_consistent   <- c(1L, 2L, NA, 3L, 5L)        # matches status
  toa_inconsistent <- c(1L, 2L, NA, 4L, 5L)        # node 4 disagrees

  expect_warning(
    x1 <- new_diffnet(gr, toa = toa_consistent, status = status,
                      t0 = 1L, t1 = 5L),
    "Both -toa- and -status- supplied"
  )
  expect_equal(as.integer(x1$toa), c(1L, 2L, NA, 3L, 5L))

  expect_warning(
    x2 <- new_diffnet(gr, toa = toa_inconsistent, status = status,
                      t0 = 1L, t1 = 5L),
    "does NOT match"
  )
  # Status wins: derived toa[4] = 3, not 4
  expect_equal(as.integer(x2$toa), c(1L, 2L, NA, 3L, 5L))
})

# ----------------------------------------------------------------------------
# Constructor: validation
# ----------------------------------------------------------------------------

test_that("status validation errors are informative", {
  gr <- mk_graph()

  # Wrong number of rows
  expect_error(
    new_diffnet(gr, status = matrix(0L, nrow = 4L, ncol = 5L),
                t0 = 1L, t1 = 5L),
    "rows"
  )
  # Non-binary entries
  expect_error(
    new_diffnet(gr, status = matrix(c(0L, 0L, 0L, 0L, 0L,
                                      0L, 0L, 0L, 0L, 0L,
                                      0L, 0L, 2L, 0L, 0L,
                                      0L, 0L, 0L, 0L, 0L,
                                      0L, 0L, 0L, 0L, 0L),
                                    nrow = 5L, byrow = TRUE),
                t0 = 1L, t1 = 5L),
    "0 or 1"
  )

  # Neither toa nor status
  expect_error(new_diffnet(gr), "either -toa- or -status-")
})

# ----------------------------------------------------------------------------
# Accessors: toa(x), tod(x), toa_all(x), tod_all(x)
# ----------------------------------------------------------------------------

test_that("toa(x) returns the same vector as x$toa", {
  gr  <- mk_graph()
  toa_v <- c(1L, 2L, NA, 3L, 5L)
  x <- new_diffnet(gr, toa_v, t0 = 1L, t1 = 5L)
  expect_identical(toa(x), x$toa)
})

test_that("tod(x) returns NA for absorbing diffnet", {
  gr <- mk_graph()
  toa_v <- c(1L, 2L, NA, 3L, 5L)
  x <- new_diffnet(gr, toa_v, t0 = 1L, t1 = 5L)
  expect_equal(as.integer(tod(x)), rep(NA_integer_, 5L))
})

test_that("tod(x) returns first recovery time for non-absorbing diffnet", {
  gr <- mk_graph()
  status <- rbind(
    c(1L, 1L, 0L, 0L, 0L),    # node 1: recovers at t=3
    c(0L, 1L, 1L, 0L, 0L),    # node 2: recovers at t=4
    c(0L, 0L, 0L, 0L, 0L),    # node 3: never adopts
    c(0L, 0L, 1L, 1L, 1L),    # node 4: absorbing
    c(0L, 0L, 0L, 0L, 1L)     # node 5: absorbing (last period)
  )
  x <- new_diffnet(gr, status = status, t0 = 1L, t1 = 5L)
  expect_equal(as.integer(tod(x)),
               c(3L, 4L, NA, NA, NA))
})

test_that("tod(x) reports only the FIRST recovery in multi-cycle", {
  gr <- mk_graph(n = 3L, T = 6L, seed = 7L)
  status <- rbind(
    c(1L, 1L, 0L, 0L, 1L, 1L),    # adopt -> recover t=3 -> re-adopt t=5
    c(0L, 0L, 0L, 0L, 0L, 0L),
    c(0L, 1L, 1L, 1L, 1L, 1L)
  )
  x <- new_diffnet(gr, status = status, t0 = 1L, t1 = 6L)
  expect_equal(as.integer(tod(x)), c(3L, NA, NA))
})

test_that("toa_all(x) returns one row per fresh adoption event", {
  gr <- mk_graph(n = 3L, T = 6L, seed = 7L)
  status <- rbind(
    c(1L, 1L, 0L, 0L, 1L, 1L),
    c(0L, 0L, 0L, 0L, 0L, 0L),
    c(0L, 1L, 1L, 1L, 1L, 1L)
  )
  x <- new_diffnet(gr, status = status, t0 = 1L, t1 = 6L)

  ta <- toa_all(x)
  expect_s3_class(ta, "data.frame")
  expect_named(ta, c("node", "behavior", "episode", "time"))
  expect_equal(nrow(ta), 3L)                                 # 2 (node 1) + 1 (node 3)
  expect_equal(ta$node,     c(1L, 1L, 3L))
  expect_equal(ta$episode,  c(1L, 2L, 1L))
  expect_equal(ta$time,     c(1L, 5L, 2L))
})

test_that("tod_all(x) returns one row per recovery event", {
  gr <- mk_graph(n = 3L, T = 6L, seed = 7L)
  status <- rbind(
    c(1L, 1L, 0L, 0L, 1L, 1L),    # one recovery at t=3
    c(0L, 0L, 0L, 0L, 0L, 0L),
    c(0L, 1L, 1L, 1L, 1L, 1L)     # absorbing
  )
  x <- new_diffnet(gr, status = status, t0 = 1L, t1 = 6L)

  td <- tod_all(x)
  expect_s3_class(td, "data.frame")
  expect_named(td, c("node", "behavior", "episode", "time"))
  expect_equal(nrow(td), 1L)
  expect_equal(td$node,    1L)
  expect_equal(td$time,    3L)
})

test_that("tod_all(x) on an absorbing diffnet returns a 0-row data.frame", {
  gr <- mk_graph()
  toa_v <- c(1L, 2L, NA, 3L, 5L)
  x <- new_diffnet(gr, toa_v, t0 = 1L, t1 = 5L)
  td <- tod_all(x)
  expect_s3_class(td, "data.frame")
  expect_equal(nrow(td), 0L)
  expect_named(td, c("node", "behavior", "episode", "time"))
})

# ----------------------------------------------------------------------------
# Multi-behaviour status path
# ----------------------------------------------------------------------------

test_that("multi-behaviour status path produces parallel matrices", {
  gr <- mk_graph()
  s_q1 <- rbind(
    c(1L, 1L, 0L, 0L, 0L),
    c(0L, 1L, 1L, 0L, 0L),
    c(0L, 0L, 0L, 0L, 0L),
    c(0L, 0L, 1L, 1L, 1L),
    c(0L, 0L, 0L, 0L, 1L)
  )
  s_q2 <- rbind(
    c(0L, 0L, 0L, 0L, 1L),
    c(0L, 1L, 1L, 1L, 1L),
    c(1L, 1L, 1L, 1L, 1L),
    c(0L, 0L, 0L, 0L, 0L),
    c(0L, 0L, 1L, 1L, 1L)
  )

  x <- new_diffnet(gr, status = list(s_q1, s_q2), t0 = 1L, t1 = 5L)

  expect_type(x$cumadopt, "list")
  expect_length(x$cumadopt, 2L)
  expect_type(x$status,   "list")
  expect_length(x$status, 2L)
  expect_identical(x$status, x$cumadopt)

  expect_equal(dim(toa(x)), c(5L, 2L))
  expect_equal(unname(toa(x)[, 1]), c(1L, 2L, NA, 3L, 5L))
  expect_equal(unname(toa(x)[, 2]), c(5L, 2L, 1L, NA, 3L))

  # tod() gives the first recovery per (node, behavior); only behavior 1
  # has any recoveries (nodes 1 and 2 in s_q1).
  td <- tod(x)
  expect_equal(dim(td), c(5L, 2L))
  expect_equal(unname(td[, 1]), c(3L, 4L, NA, NA, NA))
  expect_equal(unname(td[, 2]), rep(NA_integer_, 5L))
})

test_that("toa_all(x) and tod_all(x) span behaviours in multi-behaviour", {
  gr <- mk_graph(n = 3L, T = 4L, seed = 13L)
  s_q1 <- rbind(c(1L, 0L, 1L, 0L),    # node 1: adopt/recover/adopt/recover
                c(0L, 0L, 0L, 0L),
                c(0L, 1L, 1L, 1L))
  s_q2 <- rbind(c(0L, 0L, 0L, 0L),
                c(1L, 1L, 1L, 1L),
                c(0L, 0L, 0L, 0L))
  x <- new_diffnet(gr, status = list(s_q1, s_q2), t0 = 1L, t1 = 4L)

  ta <- toa_all(x)
  expect_equal(nrow(ta), 4L)
  expect_equal(sum(ta$behavior == 1L), 3L)
  expect_equal(sum(ta$behavior == 2L), 1L)

  td <- tod_all(x)
  expect_equal(nrow(td), 2L)                          # both for node 1, behavior 1
  expect_equal(td$node,     c(1L, 1L))
  expect_equal(td$behavior, c(1L, 1L))
  expect_equal(td$episode,  c(1L, 2L))
  expect_equal(td$time,     c(2L, 4L))
})
