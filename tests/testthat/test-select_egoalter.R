context("Select ego/alter")

test_that("Simple example", {

  # Creating graph
  g <- matrix(0, ncol=3, nrow=3)
  g[cbind(1,2:3)] <- 1
  g[cbind(2:3,1)] <- 1
  g[2,3] <- 1
  g[3,2] <- 1

  g <- lapply(1:4, function(x) methods::as(g, "dgCMatrix"))

  # Checking error
  expect_error(select_egoalter(g), "provided when")

  # Cumultive adoption matrix
  cumadopt <- matrix(0,ncol=4, nrow=3)
  cumadopt[1,1:3] <- 1
  cumadopt[2,2:4] <- 1
  cumadopt[3,3:4] <- 1


  ans <- adopt_changes(g, cumadopt)
  # class(ans) <- c("diffnet_adoptChanges", "data.frame")

  # Only over unchanged
  ans2 <- summary(ans)$`unchanged (s)`
  expect_equivalent(sapply(ans2, sum), rep(6,3))

  # Matrices calculated manually
  per2 <- matrix(0, ncol=4, nrow=4)
  per2[rbind(c(1,2),c(1,4),c(2,1), c(2,4), c(4,1), c(4,2))] <- 1

  per3 <- matrix(0, ncol=4, nrow=4)
  per3[rbind(c(2,4), c(4,2), c(4,4))] <- 2

  per4 <- matrix(0, ncol=4, nrow=4)
  per4[rbind(c(3,4), c(4,3), c(4,4))] <- 2


  expect_equivalent(ans2, list(per2, per3, per4))


  # Should be equivalent
  g <- lapply(g, as.matrix)
  g_array <- array(unlist(g), dim=c(3,3,4))

  ans <- adopt_changes(g_array, cumadopt)
  ans3 <- summary(ans)$`unchanged (s)`

  expect_equal(ans2, ans3)

})
