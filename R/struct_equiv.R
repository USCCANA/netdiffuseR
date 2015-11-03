#' Conmputes Structural Equivalence
#' @param adjmat Square matrix. Adjacency matrix
#' @param v Double. Constant for SE
#' @param ... Further arguments to be passed to \code{sna::geodist}
#' @return A square matrix indicating structural equivalence.
#' @export
struct_equiv <- function(adjmat, v=1, ...) {
  geod <- sna::geodist(adjmat, inf.replace = 0, ...)
  geod[["gdist"]] <- geod[["gdist"]]/max(geod[["gdist"]])
  struct_equiv_cpp(geod[["gdist"]], v)
}
