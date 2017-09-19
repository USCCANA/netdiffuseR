#' @param vertex.size Either a numeric scalar or vector of size \eqn{n}, or any
#' of the following values: "indegree", "degree", or "outdegree" (see details).
#' @param minmax.relative.size Passed to \code{\link{rescale_vertex_igraph}}.
#' @details
#' Plotting is done via the function \code{\link[igraph:plot.igraph]{plot.igraph}}.
#'
#' When \code{vertex.size} is either of \code{"degree"}, \code{"indegree"}, or
#' \code{"outdegree"}, \code{vertex.size} will be replace with \code{dgr(.,cmode = )}
#' so that the vertex size reflects the desired degree.
#'
#' The argument \code{minmax.relative.size} is passed to \code{\link{rescale_vertex_igraph}}
#' which adjusts \code{vertex.size} so that the largest and smallest vertices
#' have a relative size of \code{minmax.relative.size[2]} and
#' \code{minmax.relative.size[1]} respectively with respect to the x-axis.
