diag_expand <- function(...) UseMethod("diag_expand")

diag_expand.default <- function(m, times, auto.lag=TRUE) {

  # Checking class
  meta <- classify_graph(m)

  # Getting the info
  d <- with(meta, c(n, n, nper))
  if (missing(times)) times <- d[3]

  if (!times)
    stop("It must be a dynamic graph. nslices() = 0.")

  # Structure
  W <- methods::new("dgCMatrix", Dim=d[1:2]*times, p=rep(0L,d[1]*times+1L))

  # Filling
  for (p in 1:times) {
    i <- ((p-1)*d[1]+1):(d[1]*p)
    j <- ((p-1)*d[2]+1):(d[2]*p)

    if (meta$class=="array") W[i,j] <- m[,,p]
    else if (meta$class=="list") W[i,j] <- m[[p]]
    else if (meta$class=="matrix") W[i,j] <- m
  }

  # Autolag
  al <- cbind(1:d[1], 1:d[2])

  W
}

diag_expand.diffnet <- function(g, ...) {
  diag_expand.default(g$graph, g$meta$nper)
}
