context("diffnet_epi subclass")

mk_diffnet <- function() {
  set.seed(42)
  gr  <- lapply(1:5, function(x) rgraph_ba(t = 4L))
  new_diffnet(gr, toa = c(1L, 2L, NA, 3L, 5L), t0 = 1L, t1 = 5L)
}

# ----------------------------------------------------------------------------
# Promotion: as_diffnet_epi()
# ----------------------------------------------------------------------------

test_that("as_diffnet_epi() with no tree promotes to subclass and seeds empty tree", {
  x  <- mk_diffnet()
  y  <- as_diffnet_epi(x)
  expect_true(is.diffnet_epi(y))
  expect_s3_class(y, c("diffnet_epi", "diffnet"))
  expect_equal(nrow(transmission_tree(y)), 0L)
  expect_setequal(
    names(transmission_tree(y)),
    c("date", "source", "target", "source_exposure_date", "virus_id", "virus")
  )
})

test_that("as_diffnet_epi() with a pre-built transmission list stores it verbatim", {
  x <- mk_diffnet()
  pre_tree <- data.frame(
    date = 1L, source = NA_integer_, target = 1L,
    source_exposure_date = NA_integer_, virus_id = 1L, virus = NA_character_,
    stringsAsFactors = FALSE
  )
  pre <- list(tree = pre_tree, pars = list(kernel = "wells-riley"))
  y <- as_diffnet_epi(x, transmission = pre)
  expect_true(is.diffnet_epi(y))
  expect_identical(transmission_tree(y), pre_tree)
  expect_equal(y$transmission$pars$kernel, "wells-riley")
})

test_that("as_diffnet_epi() rejects malformed -transmission- inputs", {
  x <- mk_diffnet()
  expect_error(as_diffnet_epi(x, transmission = 42),       "NULL or a list")
  expect_error(as_diffnet_epi(x, transmission = "tree"),   "NULL or a list")
  expect_error(as_diffnet_epi(x, transmission = list(foo = 1)),
               "NULL or a list with -tree- and -pars-")
  # Even a data.frame is rejected here on purpose -- you'd use
  # as_transmission_tree() to attach a data.frame.
  expect_error(as_diffnet_epi(x, transmission = data.frame(a = 1)),
               "NULL or a list")
})

test_that("as_diffnet_epi() refuses non-diffnet input", {
  expect_error(as_diffnet_epi(42), "must be a diffnet")
})

test_that("Promotion is monotone (idempotent)", {
  x  <- mk_diffnet()
  y  <- as_diffnet_epi(x)
  y2 <- as_diffnet_epi(y)
  expect_identical(class(y), class(y2))     # still c("diffnet_epi", "diffnet")
  expect_equal(sum(class(y2) == "diffnet_epi"), 1L)   # not duplicated
})

# ----------------------------------------------------------------------------
# is.diffnet_epi()
# ----------------------------------------------------------------------------

test_that("is.diffnet_epi() distinguishes plain diffnet from the subclass", {
  x <- mk_diffnet()
  y <- as_diffnet_epi(x)
  expect_false(is.diffnet_epi(x))
  expect_true(is.diffnet_epi(y))
  expect_false(is.diffnet_epi("not a diffnet"))
  expect_false(is.diffnet_epi(NULL))
})

# ----------------------------------------------------------------------------
# print.diffnet_epi()
# ----------------------------------------------------------------------------

test_that("print.diffnet_epi extends print.diffnet with a transmission line", {
  x <- as_diffnet_epi(mk_diffnet())
  out <- capture.output(print(x))
  # All the base diffnet print fields are still there
  expect_true(any(grepl("Dynamic network of class -diffnet-", out)))
  expect_true(any(grepl("Final prevalence", out)))
  # And the appended epi-specific line
  expect_true(any(grepl("Transmission tree", out)))
  expect_true(any(grepl("empty", out)))
})

test_that("print.diffnet_epi shows edge/seed/virus counts for a populated tree", {
  x <- mk_diffnet()
  tree <- data.frame(
    date = c(1L, 2L, 2L),
    source = c(NA_integer_, 1L, 1L),
    target = c(1L, 2L, 4L),
    source_exposure_date = c(NA_integer_, 1L, 1L),
    virus_id = c(1L, 1L, 1L),
    virus = c("flu", "flu", "flu"),
    stringsAsFactors = FALSE
  )
  y <- as_transmission_tree(x, tree)
  out <- capture.output(print(y))
  expect_true(any(grepl("Transmission tree", out)))
  expect_true(any(grepl("3 edges", out)))
  expect_true(any(grepl("1 seeds?", out)))
  expect_true(any(grepl("1 virus(es)?", out)))
})

# ----------------------------------------------------------------------------
# Inheritance: existing diffnet methods still work
# ----------------------------------------------------------------------------

test_that("diffnet methods dispatch on diffnet_epi via inheritance", {
  y <- as_diffnet_epi(mk_diffnet())

  # nnodes/nslices come from the base class
  expect_equal(nnodes(y),  5L)
  expect_equal(nslices(y), 5L)

  # toa(), tod() (M5 accessors) work on diffnet_epi as well
  expect_equal(unname(toa(y)), c(1L, 2L, NA, 3L, 5L))
  expect_equal(unname(tod(y)), rep(NA_integer_, 5L))

  # The subsetting still produces a valid object (still inheriting from both)
  y2 <- y[-3]
  expect_true(is.diffnet_epi(y2))
  expect_true(inherits(y2, "diffnet"))
})

# ----------------------------------------------------------------------------
# Base diffnet does NOT carry a $transmission slot
# ----------------------------------------------------------------------------

test_that("a fresh new_diffnet() does not include a $transmission slot", {
  x <- mk_diffnet()
  # The slot is not part of the base structure anymore (M7).
  expect_false("transmission" %in% names(x))
})
