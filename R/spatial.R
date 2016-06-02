diag_expand <- function(...) UseMethod("diag_expand")

diag_expand.default <- function(m, times=NULL, is.array=FALSE, auto.lag=TRUE) {

  # Getting the info
  d <- dim(m)

  if (is.array) times <- d[3]

  # Structure
  W <- methods::as(Matrix::sparseMatrix(
    i={}, j={},
    dims=d[1:2]*times,
    giveCsparse = TRUE), "dgCMatrix")

  # Filling
  for (p in 1:times) {
    i <- ((p-1)*d[1]+1):(d[1]*p)
    j <- ((p-1)*d[2]+1):(d[2]*p)

    if (is.array) W[i,j] <- m[,,p]
    else W[i,j] <- m
  }

  # Autolag
  al <- cbind(1:d[1], 1:d[2])

  W
}

diag_expand.diffnet <- function(g, ...) {
  d     <- rep(nnodes(g),2)
  times <- nslices(g)

  # Casket
  W <- methods::as(Matrix::sparseMatrix(i={}, j={}, dims = d*times),
                   "dgCMatrix")

  # Filling
  for (p in 1:times) {
    i <- ((p-1)*d[1]+1):(d[1]*p)
    j <- ((p-1)*d[2]+1):(d[2]*p)
    W[i,j] <- g$graph[[p]]
  }

  # Autolag
  al <- cbind(1:d[1], 1:d[2])

  W
}
