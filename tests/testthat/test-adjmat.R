################################################################################
# edgelist_to_adjmat <--> adjmat_to_edgelist
################################################################################
context("Read/Write graphs edgelist to adjmat and vice-versa")

# Directed graph with attributes
edgelist <- cbind(
  c(2:4,1),
  1:4
)
w <- rowMeans(edgelist)
set.seed(123)
tim <- sample(1:4, 4, TRUE)

EL_digraph <- list(`edgelist matrix` = edgelist, `edgelist data.frame` = as.data.frame(edgelist))

# Undirected graph, generating the data
edgelist <- cbind(
  c(2,1,4,3),
  c(1,2,3,4)
)
EL_graph <- list(`edgelist matrix` = edgelist, `edgelist data.frame` = as.data.frame(edgelist))

#-Loop over classes
for (g in names(EL_digraph)) {

  # Arguments length
  test_that(paste0("Length of inputs in edgelist_adjmat should match (detecting error) - ",g), {
    expect_error(edgelist_to_adjmat(EL_digraph[[g]], w = w[-1], undirected=FALSE), "should have the same length")
    expect_error(edgelist_to_adjmat(EL_digraph[[g]], t0 = tim[-1], undirected=FALSE), "should have the same length")
    expect_error(edgelist_to_adjmat(EL_digraph[[g]], t0 = tim, w = w[-1], undirected=FALSE), "should have the same length")
  })

  #-----------------------------------------------------------------------------
  # Comparing back and forth using matrices
  #-----------------------------------------------------------------------------

  # Data
  adjmat <- edgelist_to_adjmat(EL_digraph[[g]], undirected = FALSE)
  edgelist_recovered <- adjmat_to_edgelist(adjmat , undirected = FALSE)

  # Comparing edgelists
  # EL_digraph[[g]] <- EL_digraph[[g]][order(EL_digraph[[g]][,1]),]
  test_that(paste0("Undirected static edgelist-adjmat-edgelist (should hold) - ",g), {
    expect_equivalent(edgelist_recovered[,1:2], as.matrix(EL_digraph[[g]]))
  })

  # Comparing adjmatrices
  test_that(paste0("Undirected static adjmat-edgelist-adjmat (should hold) - ",g), {
    expect_equivalent(edgelist_recovered, adjmat_to_edgelist(as.matrix(adjmat)))
  })

  # ----------------------------------------------------------------------------
  # Dynamic graphs (explicitly): As lists
  # ----------------------------------------------------------------------------

  # Creating the output graph
  edgelist_recovered <- adjmat_to_edgelist(
    edgelist_to_adjmat(EL_digraph[[g]], t0 = tim, undirected = FALSE),
    undirected = FALSE)

  # Collapsing
  edgelist_recovered <- unique(edgelist_recovered[,1:2])
  edgelist_recovered <- edgelist_recovered[order(edgelist_recovered[,1]),]
  EL_digraph[[g]]    <- EL_digraph[[g]][order(EL_digraph[[g]][,1]),]

  test_that(paste0("Undirected dynamic edgelist-adjmat-edgelist (should hold) - ",g), {
    expect_equivalent(edgelist_recovered, as.matrix(EL_digraph[[g]]))
  })

  # ----------------------------------------------------------------------------
  # Dynamic graphs (explicitly): As arrays
  # ----------------------------------------------------------------------------
  array_recovered <-
    lapply(edgelist_to_adjmat(EL_digraph[[g]], t0 = tim, undirected = FALSE),
           as.matrix)

  dn <- list(rownames(array_recovered[[1]]), colnames(array_recovered[[1]]),
             names(array_recovered))

  array_recovered <- array(
    unlist(array_recovered), dim=c(dim(array_recovered[[1]]), length(array_recovered)),
    dimnames = dn)

  edgelist_recovered <- adjmat_to_edgelist(array_recovered, undirected = FALSE)

  times <- edgelist_recovered[,"time"]
  edgelist_recovered <- unique(edgelist_recovered[,1:2])
  edgelist_recovered <- edgelist_recovered[order(edgelist_recovered[,1]),]

  test_that(paste0("Undirected dynamic edgelist-adjmat-edgelist (should hold) - ", g), {
    expect_equivalent(edgelist_recovered, as.matrix(EL_digraph[[g]]))
    # expect_equivalent(edgelist_recovered$times, tim)
  })

}

################################################################################
# Time of adoption
################################################################################
context("Time of Adoption (toa_mat, toa_dif)")

times <- c(2001L, 2004L, 2003L, 2008L)

graph <- lapply(2001:2008, function(x) rgraph_er(4))
diffnet <- as_diffnet(graph, times)

test_that("Should warn about -times- not been integer", {
  expect_warning(toa_mat(as.numeric(times)), "will be coersed to integer")
})

test_that("Dimensions of TOA mat should be ok", {
  toa <- toa_mat(times)
  cumadopt <- apply(toa$adopt, 2, cumsum)
  expect_equal(dim(toa$adopt), c(4, length(min(times):max(times))))
  expect_equal(dim(toa$adopt), dim(toa$cumadopt), info = "adopt and cumadopt are equal dim")
  expect_equal(t(apply(toa$adopt, 1, cumsum)), toa$cumadopt, info = "cumadopt is the cumsum")
})

test_that("Passing labels should work", {
  labs <- letters[1:length(times)]
  toa <- toa_mat(times, labels=labs)
  expect_equal(rownames(toa$adopt), labs)
  expect_equal(rownames(toa$cumadopt), labs)
})

test_that("In toa_diff, its dim should be equal to the input mat", {
  expect_equal(dim(toa_diff(times)), c(4,4))
  expect_equal(dim(toa_diff(as.integer(times))), c(4,4))
  expect_equal(toa_diff(times), toa_diff(as.integer(times)))
})

test_that("Checking toa_mat output", {

  # Manual calc
  mat <- matrix(0, nrow=4, ncol=8)
  dimnames(mat) <- list(1:4, 2001:2008)
  amat <- list(adopt=mat, cumadopt=mat)
  for (i in 1:4) {
    amat$adopt[i,times[i] - 2000] <- 1
    amat$cumadopt[i,] <- cumsum(amat$adopt[i,])
  }

  expect_equal(amat, toa_mat(diffnet))
})

################################################################################
# Isolated
################################################################################
context("Isolated vertices (isolated, drop_isolated)")

# Preparing data ---------------------------------------------------------------
edgelist <- cbind(
  c(2:4,1),
  1:4
)
set.seed(123)
tim <- sample(1:4, 4, TRUE)
adjmat <- edgelist_to_adjmat(edgelist)
dynadjmat <- edgelist_to_adjmat(edgelist, t0=tim)
diffnet <- as_diffnet(dynadjmat, tim)

test_that("Finding isolated nodes", {
  # Static graphs --------------------------------------------------------------

  # Test with dgCMatrix
  iso23<-iso2<-adjmat
  iso2[2,1:4] <- 0
  iso2[1:4,2] <- 0

  iso23[c(2,3),1:4] <- 0
  iso23[1:4,c(2,3)] <- 0

  expect_equal(which(isolated(iso2)==1), 2, info = "only one (dgCMatrix)")
  expect_equal(which(isolated(iso23)==1), c(2,3), info = "two (dgCMatrix)")
  expect_equal(which(isolated(as.matrix(iso2))==1), 2, info = "only one (matrix)")
  expect_equal(which(isolated(as.matrix(iso23))==1), c(2,3), info = "two (matrix)")

  # Test with sparse matrix
  iso2 <- as(iso2, "dgCMatrix")
  iso23 <- as(iso23, "dgCMatrix")

  expect_equal(which(isolated(iso2)==1), 2, info = "only one (array)")
  expect_equal(which(isolated(iso23)==1), c(2,3), info = "two (array)")

  # Dynamic graphs -------------------------------------------------------------

  # Test with lists
  iso2 <- lapply(dynadjmat, "[<-", i=2, j=1:4, value=0)
  iso2 <- lapply(iso2, "[<-", i=1:4, j=2, value=0)

  iso23 <- lapply(dynadjmat, "[<-", i=c(2,3), j=1:4, value=0)
  iso23 <- lapply(iso23, "[<-", i=1:4, j=c(2,3), value=0)

  expect_equal(which(isolated(iso2)$isolated==1), 2, info = "only one (list)")
  expect_equal(which(isolated(iso23)$isolated==1), c(2,3), info = "two (list)")

  # Test with array
  iso2 <- array(unlist(lapply(iso2, as.matrix)), dim=c(4,4,3))
  iso23 <- array(unlist(lapply(iso23, as.matrix)), dim=c(4,4,3))

  expect_equal(which(isolated(iso2)$isolated==1), 2, info = "only one (array)")
  expect_equal(which(isolated(iso23)$isolated==1), c(2,3), info = "two (array)")

})

test_that("Dropping isolated nodes", {
  # Static graphs --------------------------------------------------------------

  # Test with matrix
  iso23<-iso2<-adjmat
  iso2[2,1:4] <- 0
  iso2[1:4,2] <- 0

  iso23[c(2,3),1:4] <- 0
  iso23[1:4,c(2,3)] <- 0

  expect_equal(dim(drop_isolated(iso2)), c(3,3))
  expect_equal(dim(drop_isolated(iso23)), c(2,2))

  # Test with sparse matrix
  iso2 <- as(iso2, "dgCMatrix")
  iso23 <- as(iso23, "dgCMatrix")

  expect_equal(dim(drop_isolated(iso2)), c(3,3), info = "only one (dgCMatrix)")
  expect_equal(dim(drop_isolated(iso23)), c(2,2), info = "two (dgCMatrix)")
  expect_equal(dim(drop_isolated(as.matrix(iso2))), c(3,3), info = "only one (matrix)")
  expect_equal(dim(drop_isolated(as.matrix(iso23))), c(2,2), info = "two (matrix)")

  # Dynamic graphs -------------------------------------------------------------

  # Test with lists
  iso2 <- lapply(dynadjmat, "[<-", i=2, j=1:4, value=0)
  iso2 <- lapply(iso2, "[<-", i=1:4, j=2, value=0)

  iso23 <- lapply(dynadjmat, "[<-", i=c(2,3), j=1:4, value=0)
  iso23 <- lapply(iso23, "[<-", i=1:4, j=c(2,3), value=0)

  expect_equal(dim(drop_isolated(iso2)[[1]]), c(3,3))
  expect_equal(dim(drop_isolated(iso23)[[1]]), c(2,2))

  # Test with array
  iso2 <- array(unlist(lapply(iso2, as.matrix)), dim=c(4,4,3),
                dimnames=list(1:4,1:4, names(iso2)))
  iso23 <- array(unlist(lapply(iso23, as.matrix)), dim=c(4,4,3),
                 dimnames=list(1:4,1:4, names(iso2)))

  dim_list <- function(x) {
    d <- dim(x[[1]])
    c(d[1],d[2], length(x))
  }

  expect_equal(dim_list(drop_isolated(iso2)), c(3,3,3))
  expect_equal(dim_list(drop_isolated(iso23)), c(2,2,3))
})
