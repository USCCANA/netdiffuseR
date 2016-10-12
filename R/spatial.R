#' Creates a square matrix suitable for spatial statistics models.
#'
#' @param g A graph.
#' @param nper Integer scalar. Number of time periods of the graph.
#' @param self Logical scalar. When \code{TRUE} non-zero diagonal elements are used.
#' @param valued Logical scalar. When \code{FALSE} \code{g} is treated as a binary graph.
#' @param ... Further arguments to be passed to the method.
#' @return A square matrix of class \code{\link[Matrix:dgCMatrix-class]{dgCMatrix}} of
#' size \code{(nnode(g)*nper)^2}
#' @export
diag_expand <- function(...) UseMethod("diag_expand")

.diag_expand <- function(g, nper,
                                self=getOption("diffnet.self"),
                                valued=getOption("diffnet.valued")) {

  # Checking class
  meta <- classify_graph(g)

  # Getting the info
  d <- with(meta, c(n, n, nper))
  if (missing(nper)) nper <- d[3]

  if (!nper)
    stop("It must be a dynamic graph. nslices() = 0.")

  # checking options
  if (!self) {
    g <- Map(function(x) {
      Matrix::diag(x) <- rep(0,nnodes(g))
      x
    }, x = g)
  }
  if (!valued)
    g <- Map(function(x) {
    x@x <- rep(1, length(x@x))
    x
    }, x=g)

  # Structure
  W <- methods::new("dgCMatrix", Dim=d[1:2]*nper, p=rep(0L,d[1]*nper+1L))

  # Filling
  for (p in 1:nper) {
    i <- ((p-1)*d[1]+1):(d[1]*p)
    j <- ((p-1)*d[2]+1):(d[2]*p)

    W[i,j] <- g[[p]]
  }

  # Autolag
  al <- cbind(1:d[1], 1:d[2])

  W
}

#' @export
#' @rdname diag_expand
diag_expand.list <- function(g, self=getOption("diffnet.self"),
                                valued=getOption("diffnet.valued"), ...) {
  .diag_expand(g, length(g), self, valued)
}


#' @export
#' @rdname diag_expand
diag_expand.diffnet <- function(g, self=getOption("diffnet.self"),
                                valued=getOption("diffnet.valued"), ...) {
  .diag_expand(g$graph, g$meta$nper, self, valued)
}

#' @export
#' @rdname diag_expand
diag_expand.matrix <- function(g, nper, self=getOption("diffnet.self"),
                            valued=getOption("diffnet.valued"), ...) {

  .diag_expand(list(methods::as(g, "dgCMatrix")), nper, self, valued)
}


#' @export
#' @rdname diag_expand
diag_expand.array <- function(g, self=getOption("diffnet.self"),
                              valued=getOption("diffnet.valued"), ...) {

  g <- apply(g, 3, function(x) methods::as(x, "dgCMatrix"))
  diag_expand(g, nslices(g), self, valued)

}

#' @export
#' @rdname diag_expand
diag_expand.dgCMatrix <- function(g, nper, self=getOption("diffnet.self"),
                               valued=getOption("diffnet.valued"), ...) {

  .diag_expand(list(g), nper, self, valued)
}
