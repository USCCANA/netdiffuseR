context("Sparse matrices indexing in Cpp")

test_that("Are we recovering the right coords?", {
  set.seed(13131311)

  n <- 10L
  N <- 100L
  O <- NULL

  # Empty matrix
  x0 <- methods::new("dgCMatrix", Dim=c(n,n), p=rep(0L,n+1L))
  for (j in 1L:((n-1L)*2)) {
    for (i in 1L:N) {
      x1 <- x0

      # Filling place
      x1[sample.int(n*n,j)] <- 1L

      # Getting the indexes
      index <- netdiffuseR:::sparse_indexes(x1) + 1L

      # if (!(i %% 10)) print(x)

      # Reconstructing
      x2 <- as.matrix(x0)
      x2[index] <- 1L

      O <- c(O, identical(as.matrix(x1), x2))
    }
  }

  # print(x1)
  # print(methods::as(x2, "dgCMatrix"))

  expect_true(all(O))

})
