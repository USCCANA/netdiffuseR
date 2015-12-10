context("Misc functions")

test_that("The recode function for data.frame works", {
  edgelist <- data.frame(
    ego   = c(1, 10, 8),
    alter = c(10, 2, 1), stringsAsFactors = FALSE
  )

  expected_edgelist <- data.frame(
    ego   = c(1, 2, 4),
    alter = c(2, 3, 1), stringsAsFactors = FALSE
  )

  expect_equivalent(recode(edgelist), expected_edgelist)
})
