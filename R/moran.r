#' Computes Moran's I correlation index
#'
#' Natively built for computing Moran's I on \code{dgCMatrix} objects, this
#' routine allows computing the I on large sparse matrices (graphs). Part of
#' its implementation was based on \code{\link[ape:Moran.I]{ape::Moran.I}},
#' which computes the I for dense matrices.
#'
#' @param x Numeric vector of size \eqn{n}.
#' @param w Numeric matrix of size \eqn{n\times n}{n * n}. Weights. It can be
#'  either a object of class \code{\link{matrix}} or \code{\link[Matrix:dgCMatrix-class]{dgCMatrix}}
#'  from the \code{\link[Matrix]{Matrix}} package.
#' @param normalize.w Logical scalar. When TRUE normalizes rowsums to one (or zero).
#' @param alternative	Character String. Specifies the alternative hypothesis that
#' is tested against the null of no autocorrelation; must be of one \code{"two.sided"},
#' \code{"less"}, or \code{"greater"}.
#' @export
#' @details
#' In the case that the vector \code{x} is close to constant (degenerate random
#' variable), the statistic becomes irrelevant, and furthermore, the standard error
#' tends to be undefined (\code{NaN}).
#'
#' @family statistics
#' @family Functions for inference
#' @return A list of class \code{diffnet_moran} with the following elements:
#' \item{observed}{Numeric scalar. Observed correlation index.}
#' \item{expected}{Numeric scalar. Expected correlation index equal to \eqn{-1/(N-1)}.}
#' \item{sd}{Numeric scalar. Standard error under the null.}
#' \item{p.value}{Numeric scalar. p-value of the specified \code{alternative}.}
#' @references
#'
#' Moran's I. (2015, September 3). In Wikipedia, The Free Encyclopedia.
#' Retrieved 06:23, December 22, 2015, from \url{https://en.wikipedia.org/w/index.php?title=Moran\%27s_I&oldid=679297766}
#'
#' @examples
#'
#' \dontrun{
#'
#'   # Generating a small random graph
#'   set.seed(123)
#'   graph <- rgraph_ba(t = 4)
#'   w <- approx_geodesic(graph)
#'   x <- rnorm(5)
#'
#'   # Computing Moran's I
#'   moran(x, w)
#'
#'   # Comparing with the ape's package version
#'   ape::Moran.I(x, as.matrix(w))
#'
#' }
#'
#' @author George G. Vega Yon
moran <- function(x, w, normalize.w=TRUE, alternative = "two.sided") {
  if (!inherits(w, "matrix") & !inherits(w, "dgCMatrix"))
    stop("-w- must be either a matrix or a dgCMatrix.")

  if (any(!is.finite(x)))
    stop("-x- must have only numbers.")

  if (inherits(w, "matrix"))
    w <- methods::as(w, "dgCMatrix")

  if (normalize.w)
    w <- w/(rowSums(w) + 1e-20)

  res <- moran_cpp(x, w)

  # Computing pval
  pv  <- with(res, pnorm(observed, mean = expected, sd = sd))

  # Checking alternatives
  alternative <- match.arg(
    alternative, c("two.sided", "less", "greater")
    )

  if (alternative == "two.sided")
    pv <- if (res$observed <= res$expected)
      2 * pv
  else 2 * (1 - pv)
  if (alternative == "greater")
    pv <- 1 - pv

  # Returning
  structure(
    c(res, p.value = pv),
    class = "diffnet_moran"
  )

}

