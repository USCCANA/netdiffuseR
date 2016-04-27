#' Computes Moran's I correlation index
#' @param x Numeric vector of size \eqn{n}.
#' @param w Numeric matrix of size \eqn{n\times n}{n * n}. Weights.
#' @param normalize.w Logical scalar. When TRUE normalizes rowsums to one (or zero).
#' @export
#' @family statistics
#' @return Numeric scalar with Moran's I.
#' @references
#' Moran's I. (2015, September 3). In Wikipedia, The Free Encyclopedia.
#' Retrieved 06:23, December 22, 2015, from
#' \url{
#' https://en.wikipedia.org/w/index.php?title=Moran%27s_I&oldid=679297766
#' }
#' @examples
#' \dontrun{
#' # Generating a small random graph
#' set.seed(123)
#' graph <- rgraph_ba(t = 4)
#' w <- igraph::distances(igraph::graph_from_adjacency_matrix(graph))
#' x <- rnorm(5)
#'
#' # Computing Moran's I
#' moran(x, w)
#'
#' # Comparing with the ape's package version
#' moran(x, w/rowSums(as.array(w)))
#' ape::Moran.I(x, w)
#' }
#' @author George G. Vega Yon
moran <- function(x, w, normalize.w=TRUE) {
  if (!inherits(w, "matrix") & !inherits(w, "dgCMatrix"))
    stop("-w- must be either a matrix or a dgCMatrix.")

  if (any(!is.finite(x)))
    stop("-x- must have only numbers.")

  if (inherits(w, "matrix"))
    w <- methods::as(w, "dgCMatrix")

  if (normalize.w)
    w <- w/(rowSums(w) + 1e-15)

  moran_cpp(x, w)
}

