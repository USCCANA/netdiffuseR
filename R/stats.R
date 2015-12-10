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
#' # Creating an undirected graph
#' graph <- rgraph_ba()
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
  output <- degree_cpp(methods::as(graph, "dgCMatrix"), cmode, undirected, self)

  # Naming
  rn <- rownames(graph)
  if (!length(rn)) rn <- 1:nrow(graph)
  rownames(output) <- rn

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

  # Naming
  rn <- rownames(graph)
  if (!length(rn)) rn <- 1:nrow(graph)
  rownames(output) <- rn

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
  cn <- names(graph)
  if (!length(cn)) cn <- 1:length(graph)
  colnames(output) <- cn

  # Naming
  rn <- rownames(graph[[1]])
  if (!length(rn)) rn <- 1:nrow(graph[[1]])
  rownames(output) <- rn

  output
}

# @rdname dgr
# @export
dgr.array <- function(graph, cmode="degree", undirected=getOption("diffnet.undirected"), self=getOption("diffnet.self")) {
  n <- dim(graph)[1]
  t <- dim(graph)[3]
  output <- matrix(ncol=t, nrow=n)

  for(i in 1:t)
    output[,i] <- dgr(methods::as(graph[,,i], "dgCMatrix"), cmode, undirected, self)

  # Adding names
  cn <- dimnames(graph)[[3]]
  if (!length(cn)) cn <- 1:dim(graph)[3]
  colnames(output) <- cn

  # Naming
  rn <- dimnames(graph)[[1]]
  if (!length(rn)) rn <- 1:nrow(graph)
  rownames(output) <- rn

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

#' Network Hazard Rate
#'
#' Calculate and plot the hazard rate of the network.
#' @aliases plot_hazarrate
#' @param cumadopt A \eqn{n\times T}{n * T}. Cumulative adoption matrix obtained from
#' \code{\link{toa_mat}}
#' @param x An object of class \code{diffnet_hr}.
#' @param y ignored.
#' @param main Title of the plot
#' @param xlab x-axis lab
#' @param ylab y-axis lab
#' @param include.grid Logical scalar. When TRUE includes a grid on the plot.
#' @param bg Character scalar. Color of the points.
#' @param no.plot Logical scalar. When TRUE, suppress plotting (only returns hazard rates).
#' @param ... further arguments to be passed to \code{\link{plot}}
#' @details
#'
#' This function computes hazard rate, plots it and returns the hazard rate vector
#' invisible (so is not printed on the console). For \eqn{t>1}, hazard rate is calculated as
#'
#' \deqn{\frac{q_t - q_{t-1}}{n - q_{t-1}}}{[q(t) - q(t-1)]/[n - q(t-1)]}
#'
#' where \eqn{q_i}{q(i)} is the number of adopters in time \eqn{t}, and \eqn{n} is
#' the number of vertices in the graph.
#' @return A row vector of size \eqn{T} with hazard rates for \eqn{t>1} of class \code{diffnet_hr}.
#' The class of the object is only used by the S3 plot method.
#' @family statistics
#' @family visualizations
#' @keywords univar
#' @examples
#' # Creating a random vector of times of adoption
#' toa <- sample(2000:2005, 20, TRUE)
#'
#' # Computing cumulative adoption matrix
#' cumadopt <- toa_mat(toa)$cumadopt
#'
#' # Visualizing the hazard rate
#' hazard_rate(cumadopt)
#' @export
hazard_rate <- function(cumadopt, no.plot=FALSE, include.grid=TRUE, ...) {
  x <- hazard_rate_cpp(cumadopt)
  dimnames(x) <- list("hazard", colnames(cumadopt))
  class(x) <- c("diffnet_hr", class(x))
  if (!no.plot) plot.diffnet_hr(x, include.grid=include.grid, ...)
  invisible(x)
}

#' @rdname hazard_rate
#' @export
plot.diffnet_hr <- function(x,y, main="Hazard Rate", xlab="Time", ylab="Hazard Rate",
                            include.grid=TRUE, bg="lightblue",
                            ...) {
  plot(y=t(x), x=colnames(x), type="l", main=main, xlab=xlab, ylab=ylab, ...)
  if (include.grid) grid()
  points(y=t(x), x=colnames(x), pch=21, bg=bg)
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
#' graph <- rgraph_er(n=5, t=4)
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
