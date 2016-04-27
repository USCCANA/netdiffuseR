#' \pkg{netdiffuseR} default options
#' @details Set of default options used by the package. These can be retrieved
#' via \code{\link{getOption}} using the prefix \code{diffnet} (see examples)
#' @examples
#' getOption("diffnet.undirected")
#' getOption("diffnet.multiple")
#' getOption("diffnet.self")
#' @return The full list of options follows:
#' \item{undirected}{FALSE}
#' \item{self}{FALSE}
#' \item{multiple}{FALSE}
#' \item{tol}{1e-8 (used for package testing)}
#' \item{valued}{FALSE}
#' \item{outgoing}{TRUE}
#' \item{keep.isolates}{TRUE}
#' @name netdiffuseR-options
#' @author George G. Vega Yon
NULL
.onLoad <- function(libname, pkgname) {
  options(
    diffnet.undirected=FALSE,
    diffnet.self=FALSE,
    diffnet.multiple=FALSE,
    diffnet.tol=1e-8,
    diffnet.valued=FALSE,
    diffnet.outgoing=TRUE,
    diffnet.keep.isolates=TRUE
  )
}
