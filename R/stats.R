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
#' @aliases degree indegree outdegree
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
#' @author Vega Yon
dgr <- function(graph, cmode="degree", undirected=getOption("diffnet.undirected"), self=getOption("diffnet.self")) {

  switch (class(graph),
    matrix = dgr.matrix(graph, cmode, undirected, self),
    array = dgr.array(graph, cmode, undirected, self),
    dgCMatrix = dgr.dgCMatrix(graph, cmode, undirected, self),
    list = dgr.list(graph, cmode, undirected, self),
    diffnet = dgr.list(graph$graph, cmode, undirected = graph$meta$undirected),
    stopifnot_graph(graph)
  )

}

# @rdname dgr
# @export
dgr.matrix <- function(graph, cmode="degree", undirected=getOption("diffnet.undirected"), self=getOption("diffnet.self")) {

  # Checking dimensions
  dm <- dim(graph)
  if (dm[1] != dm[2]) stop("-graph- must be a square matrix.")

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

  # Checking dimensions
  dm <- dim(graph)
  if (dm[1] != dm[2]) stop("-graph- must be a square matrix.")

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
#' Calculates exposure to adoption over time via multiple different types of weight
#' matrices.  The basic  model is exposure to adoption by immediate neighbors
#' (outdegree) at the time period prior to ego’s adoption. This exposure can be
#' based on (1) incoming ties, (2) structural equivalence, (3) indirect ties, (4)
#' attribute weighted (5) network-metric weighted (e.g., central nodes have more
#' influence), and attribute-weighted (e.g., based on homophily or tie strength).
#'
#' @param graph A dynamic graph (see \code{\link{netdiffuseR-graphs}}).
#' @param cumadopt nxT matrix. Cumulative adoption matrix obtained from
#' \code{\link{toa_mat}}
#' @param attrs A numeric matrix of size \eqn{n\times T}{n * T}. Weighting for each time, period (see details).
#' @param alt.graph A dynamic graph that should be used instead of \code{graph} (see details).
#' @param outgoing Logical scalar. When TRUE, computed using outgoing ties.
#' @param valued Logical scalar. When FALSE, values of \code{graph} are set to one.
#' @param normalized Logical scalar. When true, the exposure will be between zero
#' and one (see details).
#' @details
#' DETAILS NOT UP TO DATE.
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
#' @return A matrix of size \eqn{n\times T}{n * T} with exposure for each node.
#' @export
#' @author Vega Yon, Dyal, Hayes & Valente
exposure <- function(graph, cumadopt, attrs = NULL, alt.graph=NULL, outgoing=TRUE,
                     valued=getOption("diffnet.valued"), normalized=TRUE) {

  # Checking cumadopt mat
  if (missing(cumadopt))
    if (!inherits(graph, "diffnet")) {
      stop("-cumadopt- should be provided when -graph- is not of class 'diffnet'")
    } else {
      cumadopt <- graph$cumadopt
    }

  # Checking attrs
  if (!length(attrs)) {
    attrs <- matrix(1, ncol=ncol(cumadopt), nrow=nrow(cumadopt))
  }

  # Checking alt graph
  if (length(alt.graph)) graph <- alt.graph

  switch (class(graph),
    array = exposure.array(graph, cumadopt, attrs, outgoing, valued, normalized),
    list = exposure.list(graph, cumadopt, attrs, outgoing, valued, normalized),
    diffnet = exposure.list(graph$graph, cumadopt, attrs, outgoing, valued, normalized),
    stopifnot_graph(graph)
  )
}

# @rdname exposure
# @export
exposure.array <- function(
  graph, cumadopt, attrs,
  outgoing, valued, normalized) {

  # Preparing the data
  n <- nrow(graph)
  t <- dim(graph)[3]
  graphl <- vector("list", t)
  for (i in 1:t)
    graphl[[i]] <- methods::as(graph[,,i], "dgCMatrix")

  # attrs can be either
  #  degree, indegree, outdegree, or a user defined vector.
  #  by default is user equal to 1
  da <- dim(attrs)
  if (!length(da)) stop("-attrs- must be a matrix of size n by T.")
  if (any(da != dim(cumadopt))) stop("Incorrect size for -attrs-. ",
                                     "It must be of size that -cumadopt-.")

  # Dimnames
  rn <- rownames(cumadopt)
  if (!length(rn)) rn <- 1:nrow(cumadopt)

  tn <- colnames(cumadopt)
  if (!length(tn)) tn <- 1:ncol(cumadopt)

  # Calculating the exposure, and asigning names
  output <- exposure_cpp(graphl, cumadopt, attrs, outgoing, valued, normalized)
  dimnames(output) <- list(rn, tn)
  output
}

# @rdname exposure
# @export
exposure.list <- function(
  graph, cumadopt, attrs,
  outgoing, valued, normalized) {

  # attrs can be either
  #  degree, indegree, outdegree, or a user defined vector.
  #  by default is user equal to 1
  da <- dim(attrs)
  if (!length(da)) stop("-attrs- must be a matrix of size n by T.")
  if (any(da != dim(cumadopt))) stop("Incorrect size for -attrs-. ",
                                     "It must be of size that -cumadopt-.")

  n <- nrow(graph[[1]])
  t <- length(graph)
  output <- exposure_cpp(graph, cumadopt, attrs, outgoing, valued, normalized)

  rn <- rownames(cumadopt)
  if (!length(rn)) rn <- 1:nrow(cumadopt)

  tn <- colnames(cumadopt)
  if (!length(tn)) tn <- 1:ncol(cumadopt)

  dimnames(output) <- list(rn, tn)
  output
}

#' Cummulative count of adopters
#'
#' For each time period, calculates the number of adopters, the proportion of adopters,
#' and the adoption rate.
#'
#' @param obj A \eqn{n\times T}{n * T} matrix (Cumulative adoption matrix obtained from
#' \code{\link{toa_mat}}) or a \code{\link{diffnet}} object.
#' @details
#'
#' The rate of adoption--returned in the 3rd row out the resulting
#' matrix--is calculated as
#'
#' \deqn{\frac{q_t - q_{t-1}}{q_{t-1}}}{[q(t) - q(t-1)]/q(t-1)}
#'
#' where \eqn{q_i}{q(i)} is the number of adopters in time \eqn{t}. Note that
#' it is only calculated fot \eqn{t>1}.
#' @return A \eqn{3\times T}{3 * T} matrix, where its rows contain the number of adoptes, the proportion of
#' adopters and the rate of adoption respectively, for earch period of time.
#' @family statistics
#' @keywords univar
#' @export
#' @author Vega Yon, Dyal, Hayes & Valente
cumulative_adopt_count <- function(obj) {

  if (inherits(obj, "diffnet")) obj <- obj$cumadopt

  x <- cumulative_adopt_count_cpp(obj)
  dimnames(x) <- list(c("num", "prop", "rate"), colnames(obj))
  return(x)
}

#' Network Hazard Rate
#'
#' The hazard rate is the instantaneous probability of adoption at each time
#' representing the likelihood members will adopt at that time (Allison 1984).
#' The shape of the hazard rate indicates the pattern of new adopters over time.
#' Rapid diffusion with convex cumulative adoption curves will have hazard functions
#' that peak early and decay over time whereas slow concave cumulative adoption
#' curves will have hazard functions that are low early and rise over time.
#' Smooth hazard curves indicate constant adoption whereas those that oscillate
#' indicate variability in adoption behavior over time.
#' @aliases plot_hazarrate
#' @param obj A \eqn{n\times T}{n * T} matrix (Cumulative adoption matrix obtained from
#' \code{\link{toa_mat}}) or a \code{\link{diffnet}} object.
#' @param x An object of class \code{diffnet_hr}.
#' @param y ignored.
#' @param main Character scalar. Title of the plot
#' @param xlab Character scalar. x-axis label.
#' @param ylab Character scalar. y-axis label.
#' @param include.grid Logical scalar. When TRUE includes a grid on the plot.
#' @param bg Character scalar. Color of the points.
#' @param pch Integer scalar. See \code{\link{par}}.
#' @param type Character scalar. See \code{\link{par}}.
#' @param no.plot Logical scalar. When TRUE, suppress plotting (only returns hazard rates).
#' @param add Logical scalar. When TRUE it adds the hazard rate to the current plot.
#' @param ylim Numeric vector. See \code{\link{plot}}.
#' @param ... further arguments to be passed to \code{\link{par}}
#' @details
#'
#' This function computes hazard rate, plots it and returns the hazard rate vector
#' invisible (so is not printed on the console). For \eqn{t>1}, hazard rate is calculated as
#'
#' \deqn{\frac{q_t - q_{t-1}}{n - q_{t-1}}}{[q(t) - q(t-1)]/[n - q(t-1)]}
#'
#' where \eqn{q_i}{q(i)} is the number of adopters in time \eqn{t}, and \eqn{n} is
#' the number of vertices in the graph.
#'
#' In survival analysis, hazard rate is defined formally as
#'
#' \deqn{%
#' \lambda(t)=\lim_{h\to +0}\frac{F(t+h)-F(t)}{h}\frac{1}{1-F(t)} %
#' }{%
#' \lambda(t-1)= lim (t -> +0) [F(t+h)-F(t)]/h * 1/[1-F(t)] %
#' }
#'
#' Then, by approximating \eqn{h=1}, we can rewrite the equation as
#'
#' \deqn{%
#' \lambda(t)=\frac{F(t+1)-F(t)}{1-F(t)} %
#' }{%
#' \lambda(t-1)= [F(t+1)-F(t)]/[1-F(t)] %
#' }
#'
#' Furthermore, we can estimate \eqn{F(t)}, the probability of not having adopted
#' the innovation in time \eqn{t}, as the proportion of adopters in that time, this
#' is \eqn{F(t) \sim q_t/n}{F(t) ~ q(t)/n}, so now we have
#'
#' \deqn{%
#' \lambda(t)=\frac{q_{t+1}/n-q_t/n}{1-q_t/n} = \frac{q_{t+1} - q_t}{n - q_t} %
#' }{%
#' \lambda(t-1)= [q(t+1)/n-q(t)/n]/[1-q(t)/n] = [q(t+1) - q(t)]/[n - q(t)] %
#' }
#'
#' As showed above.
#'
#' The \code{plot_hazard} function is an alias for the \code{plot.diffnet_hr} method.
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
#' @references
#' Allison, P. (1984). Event history analysis regression for longitudinal event
#' data. Beverly Hills: Sage Publications.
#'
#' Wooldridge, J. M. (2010). Econometric Analysis of Cross Section and Panel Data
#' (2nd ed.). Cambridge: MIT Press.
#' @export
#' @author Vega Yon, Dyal, Hayes & Valente
hazard_rate <- function(obj, no.plot=FALSE, include.grid=TRUE, ...) {

  if (inherits(obj, "diffnet")) obj <- obj$cumadopt

  x <- hazard_rate_cpp(obj)
  dimnames(x) <- list("hazard", colnames(obj))
  class(x) <- c("diffnet_hr", class(x))
  if (!no.plot) plot.diffnet_hr(x, include.grid=include.grid, ...)
  invisible(x)
}

#' @rdname hazard_rate
#' @export
plot_hazard <- function(x,main="Hazard Rate", xlab="Time", ylab="Hazard Rate", type="b",
                        include.grid=TRUE, bg="lightblue", add=FALSE, ylim=c(0,1), pch=21,
                        ...) {
  hr <- hazard_rate(x, no.plot = TRUE)
  plot.diffnet_hr(x=hr, main=main, xlab=xlab, ylab=ylab, type=type, include.grid=include.grid, bg=bg,
                  add=add, ylim=ylim, pch=pch, ...)
}

#' @rdname hazard_rate
#' @export
plot.diffnet_hr <- function(x,y=NULL, main="Hazard Rate", xlab="Time",
                            ylab="Hazard Rate", type="b",
                            include.grid=TRUE, bg="lightblue", pch=21, add=FALSE, ylim=c(0,1),
                            ...) {

  if (add) {
    lines(y=t(x), x=colnames(x), type=type, bg=bg, pch=pch, ...)
  } else {
    plot(y=t(x), x=colnames(x), type=type, main=main, xlab=xlab, ylab=ylab,
         ylim=ylim, bg=bg, pch=pch,...)
    if (include.grid) grid()
  }

  invisible(x)
}

#' Retrive threshold levels from the exposure matrix
#'
#' Threshold as the exposure of vertex by the time of the adoption (see \code{\link{exposure}}).
#'
#' @param obj Either a \eqn{n\times T}{n * T} matrix (eposure to the innovation obtained from
#' \code{\link{exposure}}) or a \code{diffnet} object.
#' @param times Integer vector. Indicating the time of adoption of the innovation.
#' @param t0 Integer scalar. See \code{\link{toa_mat}}.
#' @param ... Further arguments to be passed to \code{\link{exposure}}.
#' @return A vector of size \eqn{n} indicating the threshold for each node.
#' @family statistics
#' @seealso Threshold can be visualized using \code{\link{plot_threshold}}
#' @keywords univar
#' @examples
#' # Generating a random graph with random Times of Adoption
#' set.seed(783)
#' toa <- sample.int(4, 5, TRUE)
#' graph <- rgraph_er(n=5, t=max(toa) - min(toa) + 1)
#'
#' # Computing exposure using Structural Equivalnece
#' adopt <- toa_mat(toa)
#' se <- struct_equiv(graph)
#' se <- lapply(se, function(x) methods::as((x$SE)^(-1), "dgCMatrix"))
#' expo <- exposure(graph, adopt$cumadopt, alt.graph=se)
#'
#' # Retrieving threshold
#' threshold(expo, toa)
#'
#' # We can do the same by creating a diffnet object
#' diffnet <- as_diffnet(graph, toa)
#' threshold(diffnet, alt.graph=se)
#' @export
#' @author Vega Yon, Dyal, Hayes & Valente
threshold <- function(obj, times, t0=min(times, na.rm = TRUE), ...) {

  if (inherits(obj, "diffnet")) {
    t0 <- min(obj$meta$pers)
    times <- obj$toa
    obj <- exposure(obj, ...)
  } else {
    if (missing(times))
      stop("-times- should be provided when -obj- is not of class 'diffnet'")
  }

  times <- times - t0 + 1L
  output <- threshold_cpp(obj, times)
  dimnames(output) <- list(rownames(obj), "threshold")
  output
}
