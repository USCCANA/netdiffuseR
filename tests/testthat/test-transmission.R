context("Transmission tree slot")

mk_diffnet <- function() {
  set.seed(42)
  gr  <- lapply(1:5, function(x) rgraph_ba(t = 4L))
  toa <- c(1L, 2L, NA, 3L, 5L)
  new_diffnet(gr, toa, t0 = 1L, t1 = 5L)
}

test_that("get_transmissions() returns an empty data.frame when unset", {
  x <- mk_diffnet()
  tr <- get_transmissions(x)
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
  tr <- get_transmissions(y)

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
  tr <- get_transmissions(y)

  expect_true(all(tr$virus_id == 1L))
  expect_true(all(is.na(tr$virus)))
})

test_that("get_transmissions() requires a diffnet", {
  expect_error(get_transmissions(42), "must be a diffnet")
  expect_error(as_transmission_tree(42, data.frame()), "must be a diffnet")
})
