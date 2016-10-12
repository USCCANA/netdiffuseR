context("Silly tests for imports.R")

test_that("release_questions", {
  expect_equal(length(netdiffuseR:::release_questions()), 4L)
})
