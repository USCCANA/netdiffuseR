#' Indegree, outdegree and degree of the vertices
#'
#' Computes the requested degree measure for each node in the graph.
#'
#' @param graph Any class of accepted graph format (see \code{\link{netdiffuseR-graphs}}).
#' @param cmode Character. Either "indegree", "outdegree" or "degree".
#' @param undirected Logical. TRUE when the graph is undirected.
#' @param self Logical. TRUE when self edges should not be considered.
#' @return Either a numeric vector of size \eqn{n}{n} with the degree of each node (if graph is
#' a matrix), or a matrix of size \eqn{n\times T}{n * T}.
#' @export
#' @family statistics
#' @keywords univar
#' @aliases degree indegree outrdegree
#' @examples
#' # Creating a directed graph
#' graph <- rand_graph(undirected=FALSE)
#' graph
#'
#' # Comparing degree measurements
#'  data.frame(
#'    In=dgr(graph, "indegree", undirected = FALSE),
#'    Out=dgr(graph, "outdegree", undirected = FALSE),
#'    Degree=dgr(graph, "degree", undirected = FALSE)
#'  )
dgr <- function(graph, cmode="degree", undirected=getOption("diffnet.undirected"), self=getOption("diffnet.self")) {
  switch (class(graph),
    matrix = dgr.matrix(graph, cmode, undirected, self),
    array = dgr.array(graph, cmode, undirected, self),
    dgCMatrix = dgr.dgCMatrix(graph, cmode, undirected, self),
    list = dgr.list(graph, cmode, undirected, self)
  )
  # UseMethod("dgr")
}

# @rdname dgr
# @export
dgr.matrix <- function(graph, cmode="degree", undirected=getOption("diffnet.undirected"), self=getOption("diffnet.self")) {

  # Retrieving the number
  if      (cmode=="indegree")  cmode <- 0
  else if (cmode=="outdegree") cmode <- 1
  else if (cmode=="degree")    cmode <- 2
  else stop('Invalid -cmode- ',cmode,'. Should be either ',
            '"indegree", "outdegree" or "degree".')

  # Computing degree
  output <- degree_cpp(Matrix::Matrix(graph), cmode, undirected, self)
  if (length(dimnames(graph)[[1]]))
    rownames(output) <- dimnames(graph)[[1]]

  output
}

# @rdname dgr
# @export
dgr.dgCMatrix <- function(graph, cmode="degree", undirected=getOption("diffnet.undirected"), self=getOption("diffnet.self")) {

  # Retrieving the number
  if      (cmode=="indegree")  cmode <- 0
  else if (cmode=="outdegree") cmode <- 1
  else if (cmode=="degree")    cmode <- 2
  else stop('Invalid -cmode- ',cmode,'. Should be either ',
            '"indegree", "outdegree" or "degree".')

  # Computing degree
  output <- degree_cpp(graph, cmode, undirected, self)
  rownames(output) <- rownames(graph)

  output
}

# @rdname dgr
# @export
dgr.list <- function(graph, cmode="degree", undirected=getOption("diffnet.undirected"), self=getOption("diffnet.self")) {
  n <- ncol(graph[[1]])
  t <- length(graph)
  output <- matrix(ncol=t, nrow=n)

  for(i in 1:t)
    output[,i] <- dgr(graph[[i]], cmode, undirected, self)

  # Adding names
  if (length(names(graph)))
    colnames(output) <- names(graph)

  if (length(rownames(graph[[1]])))
    rownames(output) <- rownames(graph[[1]])

  output
}

# @rdname dgr
# @export
dgr.array <- function(graph, cmode="degree", undirected=getOption("diffnet.undirected"), self=getOption("diffnet.self")) {
  n <- dim(graph)[1]
  t <- dim(graph)[3]
  output <- matrix(ncol=t, nrow=n)

  for(i in 1:t)
    output[,i] <- dgr(Matrix::Matrix(graph[,,i]), cmode, undirected, self)

  # Adding names
  if (length(dimnames(graph)[[3]]))
    colnames(output) <- dimnames(graph)[[3]]

  if (length(dimnames(graph)[[1]]))
    rownames(output) <- dimnames(graph)[[1]]

  output
}

#' Ego exposure
#'
#' Calculates the level of exposure to the innovation that each node has in the
#' graph.
#'
#' @param graph A dynamic graph (see \code{\link{netdiffuseR-graphs}}).
#' @param cumadopt nxT matrix. Cumulative adoption matrix obtained from
#' \code{\link{toa_mat}}
#' @param wtype Integer. Weighting type (see details).
#' @param v Double. Constant for Structural Equivalence.
#' @param undirected Logical. TRUE if the graph is undirected.
#' @param normalized Logical. When true, the exposure will be between zero
#' and one (see details).
#' @details
#'
#' When \code{wtype=0} (default), exposure is defined as follows
#'
#' \deqn{E_{n(t)}=\frac{S_{nn}\times A_{n(t)}}{S_{n+}}}{%
#'       E(n,t)  =[     S(n,n) x     A(n,t)] / S(n+)}
#'
#' where \eqn{E_{n(t)}}{E(n,t)} is the exposure of the network at time t,
#' \eqn{S_{nn}}{S(n,n)} is the social network, \eqn{A_{n(t)}}{A(n,t)} is the
#' cumulative adoption matrix at time t, and \eqn{S_{n+}}{S(n+)} is the row-wide
#' sum of the graph.
#'
#' If \code{wtype=1} (Structural Equivalence), exposure is now computed as:
#'
#' \deqn{E_{n(t)}=\frac{SE_{nn}^{-1}\times A_{n(t)}}{SE_{n+}^{-1}}}{%
#'       E(n,t)  =[     SE(n,n)^[-1] x     A(n,t)] / SE(n+)^[-1]}
#'
#' where SE stands for Structural Equivalence (see \code{\link{struct_equiv}}).
#' Otherwise, for any value above of \code{wtype}--2, 3 or 4, which stands for
#' indegree, outdegree and degree respectively-- is calculated accordingly to
#'
#' \deqn{E_{n(t)}=\frac{S_{nn}\times (A_{n(t)} / D_n)}{S_{n+}}}{%
#'       E(n,t)  =[     S(n,n) x     [A(n,t)/D(n)] / S(n+)}
#'
#' where \eqn{D_n}{D(n)} is a column vector of size n containing the degree of
#' each node.
#'
#' Finally, note that whenever \code{normalized=TRUE}, the resulting output is only the
#' numerator of the above formulas.
#' @references
#' Burt, R. S. (1987). "Social Contagion and Innovation: Cohesion versus Structural
#' Equivalence". American Journal of Sociology, 92(6), 1287.
#' \url{http://doi.org/10.1086/228667}
#'
#' Valente, T. W. (1995). "Network models of the diffusion of innovations"
#'  (2nd ed.). Cresskill N.J.: Hampton Press.
#'
#' @family statistics
#' @keywords univar
#' @return A matrix of size nxT with exposure for each node.
#' @export
exposure <- function(graph, cumadopt, wtype = 0, v = 1.0, undirected=getOption("diffnet.undirected"), normalized=TRUE) {
  switch (class(graph),
    array = exposure.array(graph, cumadopt, wtype, v, undirected, normalized),
    list = exposure.list(graph, cumadopt, wtype, v, undirected, normalized)
  )
}

# @rdname exposure
# @export
exposure.array <- function(graph, cumadopt, wtype = 0, v = 1.0, undirected=getOption("diffnet.undirected"), normalized=TRUE) {

  # Preparing the data
  n <- nrow(graph)
  t <- dim(graph)[3]
  graphl <- vector("list", t)
  for (i in 1:t)
    graphl[[i]] <- graph[,,i]

  # Calculating the exposure, and asigning names
  output <- exposure_cpp(graphl, cumadopt, wtype, v, undirected, normalized, n, t)
  dimnames(output) <- list(rownames(graph), dimnames(graph)[[3]])
  output
}

# @rdname exposure
# @export
exposure.list <- function(graph, cumadopt, wtype = 0, v = 1.0, undirected=getOption("diffnet.undirected"), normalized=TRUE) {
  n <- nrow(graph[[1]])
  t <- length(graph)
  output <- exposure_cpp(graph, cumadopt, wtype, v, undirected, normalized, n, t)
  dimnames(output) <- list(rownames(graph[[1]]), names(graph))
  output
}

#' Cummulative count of adopters
#'
#' For each period, calculates the number of adopters, the proportion of adopters,
#' and the adoption rate.
#'
#' @param cumadopt A \eqn{n\times T}{n * T} matrix. Cumulative adoption matrix obtained from
#' \code{\link{toa_mat}}
#' @details
#'
#' The rate of adoption--returned in the 3rd row out the resulting
#' matrix--is calculated as
#'
#' \deqn{\frac{q_t - q_{t-1}}{q_{t-1}}}{[q(t) - q(t-1)]/q(t-1)}
#'
#' where \eqn{q_i}{q(i)} is the number of adopters in time \eqn{t}. Note that
#' it is only calculated fot \eqn{t>1}.
#' @return A 3xT matrix, where its rows contain the number of adoptes, the proportion of
#' adopters and the rate of adoption respectively, for earch period of time.
#' @family statistics
#' @keywords univar
#' @export
cumulative_adopt_count <- function(cumadopt) {
  x <- cumulative_adopt_count_cpp(cumadopt)
  dimnames(x) <- list(c("num", "prop", "rate"), colnames(cumadopt))
  return(x)
}

#' Graph Hazard Rate
#'
#' Computes the hazard rate of the network at each period of time.
#'
#' @param cumadopt A \eqn{n\times T}{n * T}. Cumulative adoption matrix obtained from
#' \code{\link{toa_mat}}
#' @details
#' For \eqn{t>1}, hazard rate is calculated as
#'
#' \deqn{\frac{q_t - q_{t-1}}{n - q_{t-1}}}{[q(t) - q(t-1)]/[n - q(t-1)]}
#'
#' where \eqn{q_i}{q(i)} is the number of adopters in time \eqn{t}, and \eqn{n} is
#' the number of vertices in the graph.
#' @return A row vector of size \eqn{T} with hazard rates for \eqn{t>1}.
#' @family statistics
#' @keywords univar
#' @export
hazard_rate <- function(cumadopt) {
  x <- hazard_rate_cpp(cumadopt)
  dimnames(x) <- list("hazard", colnames(cumadopt))
  x
}

#' Retrive threshold levels from the exposure matrix
#'
#' Threshold as the exposure of vertex by the time of the adoption (see \code{\link{exposure}}).
#'
#' @param exposure A \eqn{n\times T}{n * T} matrix. Exposure to the innovation obtained from
#' \code{\link{exposure}}.
#' @param times Integer vector. Indicating the time of adoption of the innovation.
#' @param times.recode Logical. TRUE when time recoding must be done.
#' @return A vector of size \eqn{n} indicating the threshold for each node.
#' @family statistics
#' @seealso Threshold can be visualized using \code{\link{plot_threshold}}
#' @keywords univar
#' @examples
#' # Generating a random graph with random Times of Adoption
#' set.seed(783)
#' graph <- rand_graph(n=5, t=4)
#' toa <- sample.int(4, 5, TRUE)
#'
#' # Computing exposure using Structural Equivalnece (wtype=1)
#' adopt <- toa_mat(toa)
#' expo <- exposure(graph, adopt$cumadopt, wtype = 1)
#'
#' # Retrieving threshold
#' threshold(expo, toa)
#' @export
threshold <- function(exposure, times, times.recode=TRUE) {
  if (times.recode) times <- times - min(times) + 1L
  output <- threshold_cpp(exposure, times)
  dimnames(output) <- list(rownames(exposure), "threshold")
  output
}
