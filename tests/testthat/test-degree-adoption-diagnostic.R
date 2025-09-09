context("Degree and Time of Adoption Diagnostic")

test_that("degree_adoption_diagnostic works with basic inputs", {
  # Create a simple test diffnet
  set.seed(123)
  dn <- rdiffnet(30, 4, seed.p.adopt = 0.2)

  # Basic test
  result <- degree_adoption_diagnostic(dn, bootstrap = FALSE)

  expect_s3_class(result, "degree_adoption_diagnostic")
  expect_named(result, c("correlations", "bootstrap", "call", "degree_strategy", "sample_size"))
  expect_named(result$correlations, c("indegree_toa", "outdegree_toa"))
  expect_type(result$correlations, "double")
  expect_length(result$correlations, 2)
  expect_null(result$bootstrap)
  expect_equal(result$degree_strategy, "mean")
  expect_gt(result$sample_size, 0)
})

test_that("degree_adoption_diagnostic works with different degree strategies", {
  set.seed(456)
  dn <- rdiffnet(25, 3, seed.p.adopt = 0.3)

  result_mean <- degree_adoption_diagnostic(dn, degree_strategy = "mean", bootstrap = FALSE)
  result_first <- degree_adoption_diagnostic(dn, degree_strategy = "first", bootstrap = FALSE)
  result_last <- degree_adoption_diagnostic(dn, degree_strategy = "last", bootstrap = FALSE)

  expect_equal(result_mean$degree_strategy, "mean")
  expect_equal(result_first$degree_strategy, "first")
  expect_equal(result_last$degree_strategy, "last")

  # Different strategies should potentially give different results
  expect_true(
    !all(result_mean$correlations == result_first$correlations) ||
    !all(result_mean$correlations == result_last$correlations)
  )
})

test_that("degree_adoption_diagnostic bootstrap works", {
  set.seed(789)
  dn <- rdiffnet(20, 3, seed.p.adopt = 0.4)

  # Skip if boot package not available
  skip_if_not_installed("boot")

  result <- degree_adoption_diagnostic(dn, bootstrap = TRUE, R = 50, conf.level = 0.9)

  expect_false(is.null(result$bootstrap))
  expect_named(result$bootstrap, c("indegree", "outdegree", "R", "boot_object"))
  expect_equal(result$bootstrap$R, 50)

  # Check bootstrap structure for indegree
  expect_named(result$bootstrap$indegree,
               c("correlation", "bias", "std_error", "conf_int", "conf_level"))
  expect_equal(result$bootstrap$indegree$conf_level, 0.9)
  expect_length(result$bootstrap$indegree$conf_int, 2)

  # Check bootstrap structure for outdegree
  expect_named(result$bootstrap$outdegree,
               c("correlation", "bias", "std_error", "conf_int", "conf_level"))
  expect_equal(result$bootstrap$outdegree$conf_level, 0.9)
  expect_length(result$bootstrap$outdegree$conf_int, 2)
})

test_that("degree_adoption_diagnostic handles edge cases", {
  # Test with diffnet that has no adopters
  set.seed(111)
  dn <- rdiffnet(100, 3, seed.p.adopt = 0.1)
  dn$toa[] <- NA  # Force no adopters

  expect_error(
    degree_adoption_diagnostic(dn, bootstrap = FALSE),
    "At least 3 adopters are required"
  )

  # Test with insufficient adopters
  set.seed(222)
  dn <- rdiffnet(100, 3, seed.p.adopt = 0.1)
  # Force only 2 adopters
  dn$toa[!is.na(dn$toa)][-(1:2)] <- NA

  expect_error(
    degree_adoption_diagnostic(dn, bootstrap = FALSE),
    "At least 3 adopters are required"
  )
})

test_that("degree_adoption_diagnostic input validation", {
  set.seed(333)
  dn <- rdiffnet(15, 3, seed.p.adopt = 0.3)

  # Test invalid toa for non-diffnet input
  expect_error(
    degree_adoption_diagnostic("not a diffnet"),
    "When 'graph' is not a diffnet, argument 'toa' must be provided."
  )

  # Test valid sparse matrix input
  A <- Matrix::rsparsematrix(10,10,0.1); A@x[] <- 1
  toa <- sample(1:3, 10, TRUE)
  x <- degree_adoption_diagnostic(A, toa = toa, bootstrap = FALSE)
  expect_s3_class(x, "degree_adoption_diagnostic")

  # Test invalid degree_strategy
  expect_error(
    degree_adoption_diagnostic(dn, degree_strategy = "invalid"),
    "'arg' should be one of"
  )

  # Test invalid bootstrap
  expect_error(
    degree_adoption_diagnostic(dn, bootstrap = "yes"),
    "'bootstrap' must be a logical scalar"
  )

  # Test invalid R
  expect_error(
    degree_adoption_diagnostic(dn, R = -1),
    "'R' must be a positive integer"
  )

  # Test invalid conf.level
  expect_error(
    degree_adoption_diagnostic(dn, conf.level = 1.5),
    "'conf.level' must be between 0 and 1"
  )
})

test_that("degree_adoption_diagnostic works with Korean Family Planning data", {
  # Test with real data
  data(kfamilyDiffNet)

  # Basic test with real data
  result <- degree_adoption_diagnostic(kfamilyDiffNet, bootstrap = FALSE)

  expect_s3_class(result, "degree_adoption_diagnostic")
  expect_gt(result$sample_size, 100)  # Korean data should have many adopters
  expect_true(all(is.finite(result$correlations)))

  # Test that correlations are within valid range
  expect_true(all(result$correlations >= -1 & result$correlations <= 1))
})

test_that("degree_adoption_diagnostic print method works", {
  set.seed(444)
  dn <- rdiffnet(20, 3, seed.p.adopt = 0.3)

  result <- degree_adoption_diagnostic(dn, bootstrap = FALSE)

  # Test that print doesn't error and returns invisibly
  expect_output(print(result), "Degree and Time of Adoption Diagnostic")
  expect_output(print(result), "Correlations:")
  expect_output(print(result), "In-degree")
  expect_output(print(result), "Out-degree")
  expect_output(print(result), "Interpretation:")

  # Test that it returns the object invisibly
  returned <- capture.output(printed_result <- print(result))
  expect_identical(result, printed_result)
})

test_that("degree_adoption_diagnostic reproduces manual calculation", {
  # Create a simple case where we can manually verify
  set.seed(555)
  dn <- rdiffnet(15, 3, seed.p.adopt = 0.4)

  result <- degree_adoption_diagnostic(dn, degree_strategy = "mean", bootstrap = FALSE)

  # Manual calculation to verify
  adopters <- !is.na(dn$toa)
  toa_adopters <- dn$toa[adopters]

  # Compute degrees manually
  indegree_full <- dgr(dn, cmode = "indegree", valued = FALSE)
  outdegree_full <- dgr(dn, cmode = "outdegree", valued = FALSE)

  indegree_mean <- rowMeans(indegree_full, na.rm = TRUE)
  outdegree_mean <- rowMeans(outdegree_full, na.rm = TRUE)

  indegree_adopters <- indegree_mean[adopters]
  outdegree_adopters <- outdegree_mean[adopters]

  # Manual correlations
  manual_cor_in <- cor(indegree_adopters, toa_adopters, use = "complete.obs")
  manual_cor_out <- cor(outdegree_adopters, toa_adopters, use = "complete.obs")

  # Compare
  expect_equal(result$correlations[["indegree_toa"]], manual_cor_in, tolerance = 1e-10)
  expect_equal(result$correlations[["outdegree_toa"]], manual_cor_out, tolerance = 1e-10)
  expect_equal(result$sample_size, sum(adopters))
})
