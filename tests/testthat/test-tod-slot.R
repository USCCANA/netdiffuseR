context("Time of disadoption (tod) slot")

mk_graph <- function(n = 5L, T = 5L, seed = 1L) {
  set.seed(seed)
  lapply(seq_len(T), function(x) rgraph_ba(t = n - 1L))
}

test_that("Default new_diffnet still has tod = NULL and transmission = NULL", {
  gr  <- mk_graph()
  toa <- c(1L, 2L, NA, 3L, 5L)
  x   <- new_diffnet(gr, toa, t0 = 1L, t1 = 5L)

  expect_null(x$tod)
  expect_null(x$transmission)
  expect_true(inherits(x, "diffnet"))
})

test_that("tod rebuilds cumadopt using [toa, tod - 1] intervals", {
  gr  <- mk_graph()
  toa <- c(1L, 2L, NA, 3L, 5L)
  tod <- c(3L, 4L, NA, NA, NA)

  x <- new_diffnet(gr, toa, tod = tod, t0 = 1L, t1 = 5L)

  expect_identical(x$tod, setNames(tod, as.character(seq_len(5L))))

  # Node 1: adopts t=1, disadopts t=3 -> cumadopt c(1,1,0,0,0)
  expect_equal(as.integer(x$cumadopt[1, ]), c(1L, 1L, 0L, 0L, 0L))
  # Node 2: adopts t=2, disadopts t=4 -> cumadopt c(0,1,1,0,0)
  expect_equal(as.integer(x$cumadopt[2, ]), c(0L, 1L, 1L, 0L, 0L))
  # Node 3: never adopts -> all zero
  expect_equal(as.integer(x$cumadopt[3, ]), rep(0L, 5L))
  # Node 4: adopts t=3, absorbing -> c(0,0,1,1,1)
  expect_equal(as.integer(x$cumadopt[4, ]), c(0L, 0L, 1L, 1L, 1L))
  # Node 5: adopts t=5, absorbing -> c(0,0,0,0,1)
  expect_equal(as.integer(x$cumadopt[5, ]), c(0L, 0L, 0L, 0L, 1L))

  # adopt matrix is preserved (single 1 at the adoption period).
  expect_equal(as.integer(rowSums(x$adopt)), c(1L, 1L, 0L, 1L, 1L))
})

test_that("tod validation errors are informative", {
  gr  <- mk_graph()
  toa <- c(1L, 2L, NA, 3L, 5L)

  # tod <= toa somewhere
  expect_error(
    new_diffnet(gr, toa, tod = c(0L, 2L, NA, NA, NA), t0 = 1L, t1 = 5L),
    "strictly greater"
  )

  # length mismatch
  expect_error(
    new_diffnet(gr, toa, tod = c(3L, 4L, NA, NA, NA, NA), t0 = 1L, t1 = 5L),
    "different lengths"
  )

  # tod defined where toa is NA
  expect_error(
    new_diffnet(gr, toa, tod = c(3L, 4L, 5L, NA, NA), t0 = 1L, t1 = 5L),
    "NA wherever"
  )

  # tod must be a vector (matrix not yet supported)
  expect_error(
    new_diffnet(gr, toa,
                tod = matrix(c(3L, 4L, NA, NA, NA), ncol = 1L),
                t0 = 1L, t1 = 5L),
    "vector of the same length"
  )
})

test_that("tod coerces non-integer input with a warning", {
  gr  <- mk_graph()
  toa <- c(1L, 2L, NA, 3L, 5L)

  expect_warning(
    x <- new_diffnet(gr, toa, tod = c(3, 4, NA, NA, NA), t0 = 1L, t1 = 5L),
    "Coercing -tod-"
  )
  expect_true(is.integer(x$tod))
})
