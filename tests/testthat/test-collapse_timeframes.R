################################################################################
# Tests for collapse_timeframes()
################################################################################

context("collapse_timeframes: collapsing longitudinal edgelists")

# Base edgelist used across most tests:
# - 2 directed pairs: (1->2) and (2->3)
# - 4 time points: 1, 2, 3, 4
# - Each pair appears twice per time point
el <- data.frame(
  sender   = c(1, 1, 1, 1, 2, 2, 2, 2),
  receiver = c(2, 2, 2, 2, 3, 3, 3, 3),
  time     = c(1, 2, 3, 4, 1, 2, 3, 4),
  weight   = c(1, 1, 1, 1, 1, 1, 1, 1)
)

################################################################################
# Block 1: Output structure
################################################################################

test_that("collapse_timeframes returns a data.frame", {
  result <- collapse_timeframes(el, ego = "sender", alter = "receiver",
                                timevar = "time", weightvar = "weight",
                                window_size = 1)
  expect_s3_class(result, "data.frame")
})

test_that("collapse_timeframes returns exactly 4 columns", {
  result <- collapse_timeframes(el, ego = "sender", alter = "receiver",
                                timevar = "time", weightvar = "weight",
                                window_size = 1)
  expect_equal(ncol(result), 4L)
})

test_that("output column names match inputs", {
  result <- collapse_timeframes(el, ego = "sender", alter = "receiver",
                                timevar = "time", weightvar = "weight",
                                window_size = 1)
  expect_named(result, c("sender", "receiver", "time", "weight"))
})

test_that("output has fewer or equal rows than input after collapsing", {
  result <- collapse_timeframes(el, ego = "sender", alter = "receiver",
                                timevar = "time", weightvar = "weight",
                                window_size = 2)
  expect_lte(nrow(result), nrow(el))
})

################################################################################
# Block 2: Binning logic (window_size)
################################################################################

test_that("window_size=1 does not merge periods", {
  result <- collapse_timeframes(el, ego = "sender", alter = "receiver",
                                timevar = "time", weightvar = "weight",
                                window_size = 1)
  expect_equal(length(unique(result$time)), 4L)
})

test_that("window_size=2 merges 4 periods into 2 bins", {
  result <- collapse_timeframes(el, ego = "sender", alter = "receiver",
                                timevar = "time", weightvar = "weight",
                                window_size = 2)
  expect_equal(length(unique(result$time)), 2L)
})

test_that("window_size=4 merges all periods into 1 bin", {
  result <- collapse_timeframes(el, ego = "sender", alter = "receiver",
                                timevar = "time", weightvar = "weight",
                                window_size = 4)
  expect_equal(length(unique(result$time)), 1L)
})

test_that("aggregated weight is sum of constituent weights", {
  # Two rows with weight=0.5 in bin 1 should aggregate to 1.0
  el2 <- data.frame(
    sender   = c(1, 1),
    receiver = c(2, 2),
    time     = c(1, 2),
    weight   = c(0.5, 0.5)
  )
  result <- collapse_timeframes(el2, ego = "sender", alter = "receiver",
                                timevar = "time", weightvar = "weight",
                                window_size = 2)
  expect_equal(result$weight, 1.0)
})

################################################################################
# Block 3: relative_time TRUE / FALSE
################################################################################

# Edgelist with a gap: time points 1, 2, 5, 6 (no 3 or 4)
el_gap <- data.frame(
  sender   = c(1, 1, 1, 1),
  receiver = c(2, 2, 2, 2),
  time     = c(1, 2, 5, 6),
  weight   = c(1, 1, 1, 1)
)

test_that("relative_time=TRUE produces a strict 1,2,... sequence", {
  result <- collapse_timeframes(el_gap, ego = "sender", alter = "receiver",
                                timevar = "time", weightvar = "weight",
                                window_size = 1, relative_time = TRUE)
  expect_equal(sort(unique(result$time)), 1:4)
})

test_that("relative_time=FALSE preserves original bin values (may have gaps)", {
  result <- collapse_timeframes(el_gap, ego = "sender", alter = "receiver",
                                timevar = "time", weightvar = "weight",
                                window_size = 1, relative_time = FALSE)
  expect_false(identical(sort(unique(result$time)), 1:4))
})

################################################################################
# Block 4: Time column parsing (integer, POSIXct, character string)
################################################################################

test_that("integer time column is handled", {
  el_int <- el
  el_int$time <- as.integer(el_int$time)
  result <- collapse_timeframes(el_int, ego = "sender", alter = "receiver",
                                timevar = "time", weightvar = "weight",
                                window_size = 2)
  expect_equal(length(unique(result$time)), 2L)
})

test_that("POSIXct time column is handled", {
  origin <- as.POSIXct("2024-01-01 00:00:00", tz = "UTC")
  el_posix <- el
  el_posix$time <- origin + (el$time - 1) * 3600  # 1 hour apart
  result <- collapse_timeframes(el_posix, ego = "sender", alter = "receiver",
                                timevar = "time", weightvar = "weight",
                                window_size = 7200)  # 2-hour windows (in seconds)
  expect_equal(length(unique(result$time)), 2L)
})

test_that("character time with time_format is parsed correctly", {
  el_chr <- el
  el_chr$time <- format(
    as.POSIXct("2024-01-01", tz = "UTC") + (el$time - 1) * 86400,
    "%Y-%m-%d"
  )
  result <- collapse_timeframes(el_chr, ego = "sender", alter = "receiver",
                                timevar = "time", weightvar = "weight",
                                window_size = 2 * 86400,  # 2-day windows
                                time_format = "%Y-%m-%d")
  expect_equal(length(unique(result$time)), 2L)
})

################################################################################
# Block 5: weightvar = NULL (count mode) vs explicit weight column
################################################################################

test_that("weightvar=NULL counts interactions as weight", {
  el_now <- data.frame(
    sender   = c(1, 1, 1),
    receiver = c(2, 2, 2),
    time     = c(1, 1, 1)
  )
  result <- collapse_timeframes(el_now, ego = "sender", alter = "receiver",
                                timevar = "time", weightvar = NULL,
                                window_size = 1)
  # 3 interactions in 1 bin -> weight should be 3
  expect_equal(result$weight, 3)
})

test_that("weightvar=NULL output column is named 'weight'", {
  result <- collapse_timeframes(el[, c("sender","receiver","time")],
                                ego = "sender", alter = "receiver",
                                timevar = "time", weightvar = NULL,
                                window_size = 1)
  expect_true("weight" %in% names(result))
})

test_that("explicit weight column is summed correctly", {
  el_w <- data.frame(
    sender   = c(1, 1),
    receiver = c(2, 2),
    time     = c(1, 1),
    w        = c(3, 7)
  )
  result <- collapse_timeframes(el_w, ego = "sender", alter = "receiver",
                                timevar = "time", weightvar = "w",
                                window_size = 1)
  expect_equal(result$w, 10)
})

################################################################################
# Block 6: Edge cases and error handling
################################################################################

test_that("NAs in time column produce a warning", {
  el_na <- el
  el_na$time[1] <- NA
  expect_warning(
    collapse_timeframes(el_na, ego = "sender", alter = "receiver",
                        timevar = "time", weightvar = "weight",
                        window_size = 1),
    "NA"
  )
})

test_that("minimal input (1 pair, 1 period) works", {
  el_min <- data.frame(sender = 1, receiver = 2, time = 1, weight = 1)
  result <- collapse_timeframes(el_min, ego = "sender", alter = "receiver",
                                timevar = "time", weightvar = "weight",
                                window_size = 1)
  expect_equal(nrow(result), 1L)
  expect_equal(result$time, 1L)
})

test_that("custom ego/alter column names are respected in output", {
  el_custom <- el
  names(el_custom) <- c("from", "to", "period", "weight")
  result <- collapse_timeframes(el_custom, ego = "from", alter = "to",
                                timevar = "period", weightvar = "weight",
                                window_size = 2)
  expect_named(result, c("from", "to", "period", "weight"))
})

test_that("output time starts at 1", {
  result <- collapse_timeframes(el, ego = "sender", alter = "receiver",
                                timevar = "time", weightvar = "weight",
                                window_size = 1)
  expect_equal(min(result$time), 1L)
})
