#' Structural Equivalence
#'
#' Computes structural equivalence between ego and alter in a network
#'
#' @param graph Any class of accepted graph format (see \code{\link{netdiffuseR-graphs}}).
#' @param v Numeric scalar. Cohesion constant (see details).
#' @param inf.replace Logical scalar. Passed to \code{\link[sna:geodist]{sna::geodist}}.
#' @param ... Further arguments to be passed to \code{\link[sna:geodist]{sna::geodist}}.
#' @family statistics
#' @keywords univar
#' @details
#'
#' Structure equivalence is computed as presented in Valente (1995), and Burt (1987),
#' in particular
#'
#' \deqn{%
#' SE_{ij} = \frac{(dmax_i - d_{ji})^v}{\sum_{k\neq i}^n(dmax_i-d_{ki})^v}
#' }{%
#' SE(ij) = [dmax(i) - d(ji)]^v/[\sum_k (dmax(i) - d(ki))^v]
#' }
#'
#' with the summation over \eqn{k\neq i}{k!=i}, and \eqn{d_{ji}}{d(ji)}, Eucledian distance in terms of geodesics, is defined as
#'
#' \deqn{%
#' d_{ji} = \left[(z_{ji} - z_{ij})^2 + \sum_k^n (z_{jk} - z_{ik})^2 +  \sum_k^n (z_{ki} - z_{kj})^2\right]^\frac{1}{2}
#' }{%
#' d(ji) = [(z(ji) - z(ij))^2 + \sum_k (z(jk) - z(ik))^2 +  \sum_k (z_(ki) - z_(kj))^2]^(1/2)
#' }
#'
#' with \eqn{z_{ij}}{z(ij)} as the geodesic (shortest path) from \eqn{i} to \eqn{j}, and
#' \eqn{dmax_i}{dmax(i)} equal to largest Euclidean distance between \eqn{i} and any other
#' vertex in the network. All summations are made over \eqn{k\not\in \{i,j\}}{k!={i,j}}
#'
#' Here, the value of \eqn{v} is interpreted as cohesion level. The higher its value,
#' the higher will be the influence that the closests alters will have over ego (see
#' Burt's paper in the reference).
#'
#' @return If \code{graph} is a static graph, a list with the following elements:
#' \item{\code{SE}}{Numeric Matrix of size \eqn{n\times n}{n * n} with Structural equivalence}
#' \item{\code{d}}{Numeric Matrix of size \eqn{n\times n}{n * n} Euclidean distances}
#' \item{\code{gdis}}{Numeric Matrix of size \eqn{n\times n}{n * n} Normalized geodesic distance}
#' In the case of dynamic graph, is a list of size \code{t} in which each element
#' contains a list as described before.
#'
#' @references Burt, R. S. (1987). "Social Contagion and Innovation: Cohesion versus
#' Structural Equivalence". American Journal of Sociology, 92(6), 1287–1335.
#' \url{http://doi.org/10.1086/228667}
#'
#' Valente, T. W. (1995). "Network models of the diffusion of innovations" (2nd ed.).
#' Cresskill N.J.: Hampton Press.
#' @export
#' @author Vega Yon, Dyal, Hayes & Valente
struct_equiv <- function(graph, v=1, inf.replace = 0,...) {
  switch (class(graph),
    matrix = struct_equiv.matrix(graph, v, inf.replace,...),
    dgCMatrix = struct_equiv.dgCMatrix(graph, v, inf.replace, ...),
    array = struct_equiv.array(graph, v, inf.replace, ...),
    list = struct_equiv.list(graph, v, inf.replace, ...),
    diffnet = struct_equiv.list(graph$graph, v, inf.replace, ...),
    stopifnot_graph(graph)
  )
}


# @rdname struct_equiv
# @export
struct_equiv.matrix <- function(graph, v, inf.replace,...) {
  geod <- sna::geodist(graph, inf.replace = inf.replace, ...)
  geod[["gdist"]] <- geod[["gdist"]]/max(geod[["gdist"]], na.rm = TRUE)
  output <- struct_equiv_cpp(methods::as(geod[["gdist"]], "dgCMatrix"), v)

  # Names
  rn <- rownames(graph)
  if (!length(rn)) rn <- 1:nrow(graph)
  output <- lapply(output, "dimnames<-", value=list(rn, rn))

  return(output)
}

# @rdname struct_equiv
# @export
struct_equiv.dgCMatrix <- function(graph, v, inf.replace,...) {
  # In order to use the SNA package functions, we need to coerce the graph
  # Into a -matrix.csc- object,
  geod <- sna::geodist(methods::as(graph, "matrix.csc"), inf.replace = inf.replace, ...)
  geod[["gdist"]] <- geod[["gdist"]]/max(geod[["gdist"]], na.rm = TRUE)
  output <- struct_equiv_cpp(methods::as(geod[["gdist"]], "dgCMatrix"), v)

  # Names
  rn <- rownames(graph)
  if (!length(rn)) rn <- 1:nrow(graph)
  output <- lapply(output, "dimnames<-", value=list(rn, rn))

  return(output)
}

# @rdname struct_equiv
# @export
struct_equiv.array <- function(graph, v, inf.replace,...) {
  t <- dim(graph)[3]
  output <- vector("list", t)
  for(i in 1:t) {
    output[[i]] <- struct_equiv.matrix(graph[,,i], v, inf.replace, ...)
  }

  # Naming
  tn <- dimnames(graph)[[3]]
  if (!length(tn)) tn <- 1:t

  names(output) <- tn

  output
}


# @rdname struct_equiv
# @export
struct_equiv.list <- function(graph, v, inf.replace, ...) {
  t <- length(graph)
  n <- nrow(graph[[1]])
  output <- vector("list", t)
  for(i in 1:t)
    output[[i]] <- struct_equiv.dgCMatrix(methods::as(graph[[i]], "dgCMatrix"), v, inf.replace, ...)

  # Naming
  tn <- dimnames(graph)[[3]]
  if (!length(tn)) tn <- 1:t

  names(output) <- tn

  output
}
