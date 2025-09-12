#' Creates a square matrix suitable for spatial statistics models.
#' @templateVar self TRUE
#' @templateVar valued TRUE
#' @template graph_template
#' @param nper Integer scalar. Number of time periods of the graph.
#' @param ... Further arguments to be passed to the method.
#' @return A square matrix of class \code{\link[Matrix:dgCMatrix-class]{dgCMatrix}} of
#' size \code{(nnode(g)*nper)^2}
#' @examples
#' # Simple example ------------------------------------------------------------
#' set.seed(23)
#' g <- rgraph_er(n=10, p=.5, t=2,undirected=TRUE)
#'
#' # What we've done: A list with 2 bernoulli graphs
#' g
#'
#' # Expanding to a 20*20 matrix with structural zeros on the diagonal
#' # and on cell 'off' adjacency matrix
#' diag_expand(g)
#' @export
diag_expand <- function(...) UseMethod("diag_expand")

.diag_expand <- function(
  graph,
  nper,
  self   = getOption("diffnet.self"),
  valued = getOption("diffnet.valued")
  ) {

  # Checking class
  meta <- classify_graph(graph)

  # Getting the info
  d <- with(meta, c(n, n, nper))
  if (missing(nper)) nper <- d[3]

  if (!nper)
    stop("It must be a dynamic graph. nslices() = 0.")

  # checking options
  if (!self) {
    graph <- Map(function(x) {
      Matrix::diag(x) <- rep(0,nnodes(graph))
      x
    }, x = graph)
  }
  if (!valued)
    graph <- Map(function(x) {
    x@x <- rep(1, length(x@x))
    x
    }, x=graph)

  # Structure
  W <- methods::new("dgCMatrix", Dim=d[1:2]*nper, p=rep(0L,d[1]*nper+1L))

  # Filling
  for (p in 1:nper) {
    i <- ((p-1)*d[1]+1):(d[1]*p)
    j <- ((p-1)*d[2]+1):(d[2]*p)

    W[i,j] <- graph[[p]]
  }

  # Autolag
  al <- cbind(1:d[1], 1:d[2])

  W
}

#' @export
#' @rdname diag_expand
diag_expand.list <- function(
    graph,
    self   = is_self(graph),
    valued = is_valued(graph),
    ...
    ) {
  .diag_expand(graph, length(graph), self, valued)
}


#' @export
#' @rdname diag_expand
diag_expand.diffnet <- function(
    graph,
    self   = is_self(graph),
    valued = is_valued(graph),
    ...
    ) {
  .diag_expand(graph$graph, graph$meta$nper, self, valued)
}

#' @export
#' @rdname diag_expand
diag_expand.matrix <- function(
    graph,
    nper,
    self   = is_self(graph),
    valued = is_valued(graph),
    ...) {

  .diag_expand(list(methods::as(graph, "dgCMatrix")), nper, self, valued)
}


#' @export
#' @rdname diag_expand
diag_expand.array <- function(
    graph,
    self   = is_self(graph),
    valued = is_valued(graph),
    ...
    ) {

  graph <- apply(graph, 3, function(x) methods::as(x, "dgCMatrix"))
  diag_expand(graph, nslices(graph), self, valued)

}

#' @export
#' @rdname diag_expand
diag_expand.dgCMatrix <- function(
    graph,
    nper,
    self   = is_self(graph),
    valued = is_valued(graph),
    ...) {

  .diag_expand(list(graph), nper, self, valued)
}
