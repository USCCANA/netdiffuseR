#' Calculates degree of each node
#' @param adjmat Square matrix. Adjacency matrix
#' @param cmode Integer. 0 for indegree, 1 for outdegree and 2 for degree.
#' @param undirected Logical. TRUE when the graph is undirected.
#' @param self Logical. TRUE when self edges should not be considered.
#' @return A numeric vector with the degree of each node.
#' @export
degree <- function(adjmat, cmode=2, undirected=TRUE, self=FALSE) {
  degree_cpp(adjmat, cmode, undirected, self)
}

#' Ego exposure
#' @param graph An array of adjacency matrices
#' @param cumadopt nxT matrix. Cumulative adoption matrix obtained from
#' \code{\link{toa_mat}}
#' @param wtype Integer. Weighting type (see details).
#' @param v Double. Constant for Structural Equivalence.
#' @param undirected Logical. TRUE if the graph is undirected.
#' @return A matrix of size nxT with exposure for each node.
#' @export
exposure <- function(graph, cumadopt, wtype = 0, v = 1.0, undirected=TRUE) {
  exposure_cpp(graph, cumadopt, wtype, v, undirected)
}

#' Cummulative count of adopters
#' @param cumadopt nxT matrix. Cumulative adoption matrix obtained from
#' \code{\link{toa_mat}}
#' @export
cumulative_adopt_count <- function(cumadopt) {
  cumulative_adopt_count_cpp(cumadopt)
}

#' Calculates hazard rate
#' @param cumadopt nxT matrix. Cumulative adoption matrix obtained from
#' \code{\link{toa_mat}}
#' @return A vector indicating the hazard rate of each node
#' @export
hazard_rate <- function(cumadopt) {
  hazard_rate_cpp(cumadopt)
}

#' Calculates threshold
#' @param exposure nxT matrix. Exposure to the innovation obtained from
#' \code{\link{exposure}}.
#' @param toe Integer vector. Indicating the time of adoption of the innovation.
#' @return A vector of size n indicating the threshold of each node.
#' @export
threshold <- function(exposure, toe) {
  threshold_cpp(exposure, toe)
}

