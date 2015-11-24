#' \pkg{netdiffuseR} default options
#' @details Set of default options used by the package. These can be retrieved
#' via \code{\link{getOption}} using the prefix \code{diffnet} (see examples)
#' @examples
#' getOption("diffnet.undirected")
#' getOption("diffnet.multiple")
#' getOption("diffnet.self")
#' @return The full list of options follows:
#' \item{undirected}{TRUE}
#' \item{self}{FALSE}
#' \item{multiple}{FALSE}
#' @name netdiffuseR-options
NULL
.onLoad <- function(libname, pkgname) {
  options(
    diffnet.undirected=TRUE,
    diffnet.self=FALSE,
    diffnet.multiple=FALSE
  )
}
