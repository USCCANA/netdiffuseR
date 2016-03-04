context("Checking diffnet-indexing")

test_that("Dynamic attributes assignment work", {
  # Preparing the data
  set.seed(1313)
  x <- rdiffnet(20, 5)

  # Calculating exposure
  expoM <- exposure(x)
  expoL <- lapply(seq_len(x$meta$nper), function(x) expoM[,x,drop=FALSE])
  expoD <- do.call(rbind, expoL)

  # Adding data
  x[["expoM"]] <- expoM
  x[["expoL"]] <- expoL
  x[["expoD"]] <- expoD

  # Checking
  eM <- x[["expoM"]]
  eL <- x[["expoL"]]
  eD <- x[["expoD"]]
  expect_equal(eM, eL)
  expect_equal(eL, eD)

})
