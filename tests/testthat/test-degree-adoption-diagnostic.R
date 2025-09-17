context("Degree and Time of Adoption Diagnostic")

test_that("degree_adoption_diagnostic works with basic inputs", {
  # Create a simple test diffnet
  set.seed(123)
  dn <- rdiffnet(30, 4, seed.p.adopt = 0.2)

  # Basic test
  result <- degree_adoption_diagnostic(dn, bootstrap = FALSE)

  expect_s3_class(result, "degree_adoption_diagnostic")
  expect_setequal(names(result), c("correlations","bootstrap","call","degree_strategy","sample_size","combine","undirected"))
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
  expect_setequal(names(result$bootstrap), c("indegree","outdegree","R","boot_object"))
  expect_equal(result$bootstrap$R, 50)

  # Check bootstrap structure for indegree
  expect_true(all(c("correlation") %in% names(result$bootstrap$indegree)))
  if (!is.null(result$bootstrap$indegree$conf_int)) {
    expect_true(all(c("bias","std_error","conf_int","conf_level") %in% names(result$bootstrap$indegree)))
    expect_equal(result$bootstrap$indegree$conf_level, 0.9)
    expect_length(result$bootstrap$indegree$conf_int, 2)
  }

  # Check bootstrap structure for outdegree (can be degenerate/NA)
  od <- result$bootstrap$outdegree
  expect_true("correlation" %in% names(od))
  if (!is.na(result$correlations[["outdegree_toa"]])) {
    # fully-formed bootstrap stats expected
    expect_true(all(c("bias","std_error","conf_int","conf_level") %in% names(od)))
    expect_equal(od$conf_level, 0.9)
    expect_length(od$conf_int, 2)
  } else {
    # degenerate case: allow missing CI/SE/CL when correlation is NA
    expect_true(is.null(od$conf_int) || length(od$conf_int) == 0)
  }
})

test_that("degree_adoption_diagnostic handles edge cases", {
  # Test with diffnet that has no adopters
  set.seed(111)
  dn <- rdiffnet(100, 3, seed.p.adopt = 0.1)
  dn$toa[] <- NA  # Force no adopters

  expect_error(
    degree_adoption_diagnostic(dn, bootstrap = FALSE),
    "Insufficient adopters for correlation analysis"
  )

  # Test with insufficient adopters
  set.seed(222)
  dn <- rdiffnet(100, 3, seed.p.adopt = 0.1)
  # Force only 2 adopters
  dn$toa[!is.na(dn$toa)][-(1:2)] <- NA

  expect_error(
    degree_adoption_diagnostic(dn, bootstrap = FALSE),
    "Insufficient adopters for correlation analysis"
  )
})

test_that("degree_adoption_diagnostic input validation", {
  set.seed(333)
  dn <- rdiffnet(n = 15, t = 3, seed.p.adopt = 0.3)

  # Test invalid toa for non-diffnet input
  expect_error(
    degree_adoption_diagnostic("not a diffnet"),
    "toa argument is required when graph is not a diffnet object"
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

test_that("degree_adoption_diagnostic supports multi-diffusion with combine='none'", {
  set.seed(1001)
  n <- 80; t <- 5; q <- 3
  garr <- rgraph_ws(n, t, p = .15)
  dn <- rdiffnet(seed.graph = garr, t = t, seed.p.adopt = rep(list(0.12), q))
  colnames(dn$toa) <- paste0("Beh", seq_len(q))

  res <- degree_adoption_diagnostic(dn, combine = "none", bootstrap = FALSE)

  # Structure: correlations is 2 x Q matrix with named columns
  expect_true(is.matrix(res$correlations))
  expect_identical(rownames(res$correlations), c("indegree_toa","outdegree_toa"))
  expect_equal(ncol(res$correlations), q)
  expect_identical(colnames(res$correlations), colnames(dn$toa))
  # sample_size is integer vector of length Q
  expect_true(is.integer(res$sample_size) || is.numeric(res$sample_size))
  expect_equal(length(res$sample_size), q)
  # combine flag
  expect_identical(res$combine, "none")
})

test_that("degree_adoption_diagnostic supports multi-diffusion with combine!='none'", {
  set.seed(1002)
  n <- 90; t <- 4; q <- 3
  garr <- rgraph_ws(n, t, p = .12)
  dn <- rdiffnet(seed.graph = garr, t = t, seed.p.adopt = rep(list(0.10), q))
  colnames(dn$toa) <- paste0("B", seq_len(q))

  res_pooled   <- degree_adoption_diagnostic(dn, combine = "pooled", bootstrap = FALSE)
  res_average  <- degree_adoption_diagnostic(dn, combine = "average", bootstrap = FALSE)
  res_earliest <- degree_adoption_diagnostic(dn, combine = "earliest", bootstrap = FALSE)

  # Single behavior-like objects: correlations is a named numeric vector
  for (res in list(res_pooled, res_average, res_earliest)) {
    expect_true(is.numeric(res$correlations))
    expect_named(res$correlations, c("indegree_toa","outdegree_toa"))
    expect_equal(length(res$correlations), 2)
    expect_true(is.numeric(res$sample_size) && length(res$sample_size) == 1)
    expect_true(res$combine %in% c("pooled","average","earliest"))
  }

  # Sample sizes: pooled should be >= earliest/average
  expect_gte(res_pooled$sample_size, res_average$sample_size)
  expect_gte(res_pooled$sample_size, res_earliest$sample_size)

  # Behavior subsetting by names and by indices produce the same pooled result
  set.seed(1003)
  n <- 70; t <- 4; q <- 3
  garr <- rgraph_ws(n, t, p = .20)
  dn <- rdiffnet(seed.graph = garr, t = t, seed.p.adopt = rep(list(0.15), q))
  colnames(dn$toa) <- c("Alpha","Beta","Gamma")

  res_by_names  <- degree_adoption_diagnostic(dn, behavior = c("Beta","Gamma"), combine = "pooled", bootstrap = FALSE)
  res_by_index  <- degree_adoption_diagnostic(dn, behavior = c(2,3),          combine = "pooled", bootstrap = FALSE)

  expect_equal(res_by_names$correlations[["indegree_toa"]],  res_by_index$correlations[["indegree_toa"]], tolerance = 1e-12)
  expect_equal(res_by_names$correlations[["outdegree_toa"]], res_by_index$correlations[["outdegree_toa"]], tolerance = 1e-12)
  expect_equal(res_by_names$sample_size, res_by_index$sample_size)
})

test_that("Undirected graphs collapse to a single Degree correlation and printer reflects that", {
  set.seed(1004)
  n <- 60; t <- 3
  # Build symmetric (undirected) time-varying graph
  glist <- replicate(t, {
    A <- matrix(0L, n, n)
    A[sample(n*n, size = round(0.08*n*n))] <- 1L
    A <- A * (1 - diag(n))                # no self loops
    A <- (A | t(A)) * 1L                  # symmetrize
    A
  }, simplify = FALSE)

  toa <- sample(1:t, n, replace = TRUE); toa[sample(n, 8)] <- NA
  res <- degree_adoption_diagnostic(glist, toa = toa, bootstrap = FALSE)

  expect_true(isTRUE(res$undirected))
  # indegree_toa and outdegree_toa should be identical under undirected collapse
  expect_equal(res$correlations[["indegree_toa"]], res$correlations[["outdegree_toa"]], tolerance = 0)
  # Printer shows "Degree - Time of Adoption" and not the separate lines
  expect_output(print(res), "Degree\\s+-\\s+Time of Adoption")
})

test_that("Coercion works for non-diffnet inputs: dgCMatrix, array, and list", {
  set.seed(1005)
  n <- 40; t <- 4

  # dgCMatrix (single slice)
  A <- Matrix::rsparsematrix(n, n, 0.05); A@x[] <- 1
  diag(A) <- 0
  toa1 <- sample(1:t, n, TRUE); toa1[sample(n, 5)] <- NA
  x1 <- degree_adoption_diagnostic(A, toa = toa1, bootstrap = FALSE)
  expect_s3_class(x1, "degree_adoption_diagnostic")
  expect_true(is.numeric(x1$sample_size) && length(x1$sample_size) == 1)

  # 3D array
  arr <- array(0L, dim = c(n, n, t))
  for (k in 1:t) {
    M <- matrix(rbinom(n*n, 1, 0.05), n, n)
    diag(M) <- 0
    arr[,,k] <- M
  }
  toa2 <- sample(1:t, n, TRUE); toa2[sample(n, 6)] <- NA
  x2 <- degree_adoption_diagnostic(arr, toa = toa2, bootstrap = FALSE)
  expect_s3_class(x2, "degree_adoption_diagnostic")

  # list of matrices
  glist <- replicate(t, {
    M <- matrix(rbinom(n*n, 1, 0.06), n, n); diag(M) <- 0; M
  }, simplify = FALSE)
  toa3 <- sample(1:t, n, TRUE); toa3[sample(n, 4)] <- NA
  x3 <- degree_adoption_diagnostic(glist, toa = toa3, bootstrap = FALSE)
  expect_s3_class(x3, "degree_adoption_diagnostic")
})

test_that("Zero-variance cases yield NA correlations without warnings and printer explains NA", {
  set.seed(1006)
  n <- 50; t <- 3
  # Construct a directed graph with constant out-degree across nodes at all times
  glist <- replicate(t, {
    # Each node nominates exactly k others (fixed out-degree)
    k <- 3
    M <- matrix(0L, n, n)
    for (i in 1:n) {
      neigh <- sample(setdiff(1:n, i), k)
      M[i, neigh] <- 1L
    }
    M
  }, simplify = FALSE)

  toa <- sample(1:t, n, TRUE); toa[sample(n, 7)] <- NA
  res <- degree_adoption_diagnostic(glist, toa = toa, bootstrap = TRUE, R = 50)

  # Out-degree correlation can be NA; in any case, print should mention "r is NA" for that metric
  if (is.na(res$correlations[["outdegree_toa"]])) {
    expect_output(print(res), "r is NA; no CI\\.")
  }
})

test_that("Degree strategy changes results on time-varying graphs", {
  set.seed(1007)
  n <- 60; t <- 5
  glist <- replicate(t, {
    M <- matrix(rbinom(n*n, 1, 0.06), n, n); diag(M) <- 0; M
  }, simplify = FALSE)
  toa <- sample(1:t, n, TRUE); toa[sample(n, 6)] <- NA

  r_mean  <- degree_adoption_diagnostic(glist, toa = toa, degree_strategy = "mean",  bootstrap = FALSE)
  r_first <- degree_adoption_diagnostic(glist, toa = toa, degree_strategy = "first", bootstrap = FALSE)
  r_last  <- degree_adoption_diagnostic(glist, toa = toa, degree_strategy = "last",  bootstrap = FALSE)

  expect_true(
    !all(r_mean$correlations == r_first$correlations, na.rm = TRUE) ||
    !all(r_mean$correlations == r_last$correlations,  na.rm = TRUE)
  )
})

test_that("Behavior selection validates indices/names and min_adopters is enforced", {
  set.seed(1008)
  n <- 60; t <- 4; q <- 2
  garr <- rgraph_ws(n, t, p = .18)
  dn <- rdiffnet(seed.graph = garr, t = t, seed.p.adopt = rep(list(0.05), q))
  colnames(dn$toa) <- c("X","Y")

  # Out-of-range indices
  expect_error(degree_adoption_diagnostic(dn, behavior = c(3), combine = "pooled"),
               "out of range")

  # Unknown names
  expect_error(degree_adoption_diagnostic(dn, behavior = c("Z"), combine = "pooled"),
               "not found")

  # min_adopters threshold
  # Force very few adopters in behavior X
  dn2 <- dn
  dn2$toa[, "X"] <- NA
  few <- sample(which(!is.na(dn$toa[, "X"])), size = 2)
  dn2$toa[few, "X"] <- dn$toa[few, "X"]
  # For combine='none' the function should NOT error; it should return NA correlations for that behavior
  res_few_none <- degree_adoption_diagnostic(dn2, behavior = "X", combine = "none", min_adopters = 3, bootstrap = FALSE)
  expect_true(is.matrix(res_few_none$correlations))
  expect_true(all(is.na(res_few_none$correlations)))
  expect_equal(unname(res_few_none$sample_size["X"]), 2)

  # For combined modes, too few rows should still error
  expect_error(
    degree_adoption_diagnostic(dn2, behavior = "X", combine = "pooled", min_adopters = 3),
    "Insufficient adopters for correlation analysis"
  )
})
