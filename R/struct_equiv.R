#' Structural Equivalence
#'
#' Computes structural equivalence between ego and alter in a network
#'
#' @param graph Any class of accepted graph format (see \code{\link{netdiffuseR-graphs}}).
#' @param v Numeric scalar. Cohesion constant (see details).
#' @param inf.replace Numeric scalar scalar. Replacing inf values obtained from \code{\link[igraph:distances]{igraph::distances}}.
#' @param groupvar Either a character scalar (if \code{graph} is diffnet), or a vector of size \eqn{n}.
#' @param ... Further arguments to be passed to \code{\link[igraph:distances]{igraph::distances}} (not valid for the print method).
#' @param mode Character scalar pased to \code{\link[igraph:distances]{igraph::distances}}
#' @param x A \code{diffnet_se} class object.
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
#' Structural equivalence can be computed either for the entire graph or by groups
#' of vertices. When, for example, the user knows before hand that the vertices
#' are distributed accross separated communities, he can make this explicit to
#' the function and provide a \code{groupvar} variable that accounts for this.
#' Hence, when \code{groupvar} is not \code{NULL} the algorithm will compute
#' structural equivalence within communities as marked by \code{groupvar}.
#'
#'
#' @return If \code{graph} is a static graph, a list with the following elements:
#' \item{\code{SE}}{Matrix of size \eqn{n\times n}{n * n} with Structural equivalence}
#' \item{\code{d}}{Matrix of size \eqn{n\times n}{n * n} Euclidean distances}
#' \item{\code{gdist}}{Matrix of size \eqn{n\times n}{n * n} Normalized geodesic distance}
#' In the case of dynamic graph, is a list of size \code{t} in which each element
#' contains a list as described before. When \code{groupvar} is specified, the
#' resulting matrices will be of class \code{\link[Matrix:dgCMatrix-class]{dgCMatrix}},
#' otherwise will be of class \code{\link{matrix}}.
#'
#' @references Burt, R. S. (1987). "Social Contagion and Innovation: Cohesion versus
#' Structural Equivalence". American Journal of Sociology, 92(6), 1287–1335.
#' \url{http://doi.org/10.1086/228667}
#'
#' Valente, T. W. (1995). "Network models of the diffusion of innovations" (2nd ed.).
#' Cresskill N.J.: Hampton Press.
#' @export
#'
#' @examples
#' # Computing structural equivalence for the fakedata -------------------------
#' data(fakesurvey)
#'
#' # Coercing it into a diffnet object
#' fakediffnet <- survey_to_diffnet(
#'    fakesurvey, "id", c("net1", "net2", "net3"), "toa", "group"
#' )
#'
#' # Computing structural equivalence without specifying group
#' se_all <- struct_equiv(fakediffnet)
#'
#' # Notice that pairs of individuals from different communities have
#' # non-zero values
#' se_all
#' se_all[[1]]$SE
#'
#' # ... Now specifying a groupvar
#' se_group <- struct_equiv(fakediffnet, groupvar="group")
#'
#' # Notice that pairs of individuals from different communities have
#' # only zero values.
#' se_group
#' se_group[[1]]$SE
#'
#'
#'
#' @author George G. Vega Yon, Stephanie R. Dyal, Timothy B, Hayes, Thomas W. Valente
struct_equiv <- function(graph, v=1, inf.replace = 0, groupvar=NULL, mode="out", ...) {

  # Checking groupvar
  if ((length(groupvar)==1) && inherits(graph, "diffnet"))
    groupvar <- graph[[groupvar]]

  output <- switch (class(graph),
    matrix = struct_equiv.matrix(graph, v, inf.replace, groupvar, mode, ...),
    dgCMatrix = struct_equiv.dgCMatrix(graph, v, inf.replace, groupvar, mode, ...),
    array = struct_equiv.array(graph, v, inf.replace, groupvar, mode, ...),
    list = struct_equiv.list(graph, v, inf.replace, groupvar, mode, ...),
    diffnet = struct_equiv.list(graph$graph, v, inf.replace, groupvar, mode, ...),
    stopifnot_graph(graph)
  )

  structure(output, class="diffnet_se", n = nnodes(graph), nper=nslices(graph),
            dyn = ifelse(class(graph) %in% c("diffnet", "list", "array"), TRUE, FALSE),
            inf.replace=inf.replace)
}

#' @export
#' @rdname struct_equiv
print.diffnet_se <- function(x, ...) {
  dyn <- attr(x, "dyn")
  cat(paste("Structural equivalence for a", ifelse(dyn, "dynamic", "static"),"graph"),
      paste(" # nodes :", attr(x, "n")),
      paste(" # of slices:",attr(x, "nper")),
      if (dyn) " Access elements via [[nslice]]$ (as a nested list). Available elements are:"
      else " Access elements via $ (as a list). Available elements are:",
      "  - SE    : Structural equivalence matrix (n x n)",
      "  - d     : Euclidean distances matrix (n x n) ",
      "  - gdist : Structural equivalence matrix (n x n)",
      sep="\n"
      )

  return(invisible(x))
}


# Apply per grouping variable
struct_equiv_by <- function(graph, v, inf.replace, groupvar, mode, ...) {
  # Checking length of the grouping variable
  if (length(groupvar) != nvertices(graph))
    stop("The length of -groupvar-, ",length(groupvar),
         " elements, must be equal to the number of vertices in the graph, ",
         nvertices(graph)," vertices.")

  # Checking that is complete
  test <- which(!complete.cases(groupvar))
  if (length(test))
    stop("-groupvar- must not have incomplete cases. Check the following elements:\n\t",
         paste0(test, collapse=", "), ".")

  # Does the graph has ids?
  if (!length(rownames(graph))) {
    rn <- 1:nvertices(graph)
    dimnames(graph) <- list(rn, rn)
  }

  # Splitting the data and computing structural equivalence per case
  G    <- unique(groupvar)
  n    <- nvertices(graph)
  SE   <- methods::new("dgCMatrix", Dim=c(n,n), p=rep(0L,n+1L))
  d    <- SE
  gdist <- SE

  rn <- NULL
  N  <- 0L
  for (g in G) {
    # Subsetting the graph
    test <- which(groupvar == g)
    subg <- graph[test, test, drop=FALSE]
    rn   <- c(rn, rownames(subg))
    ng   <- nvertices(subg)

    # Computing SE
    out <- struct_equiv(subg, v, inf.replace, groupvar=NULL, mode,...)

    # Appending values
    index <- (N + 1L):(N + ng)
    SE[index,index]   <- out$SE
    d[index,index]    <- out$d
    gdist[index,index] <- out$SE
    N <- N + ng
  }

  # Giving names
  dimnames(SE) <- list(rn, rn)
  dimnames(d) <- list(rn, rn)
  dimnames(gdist) <- list(rn, rn)

  # Ordening output
  index <- rownames(graph)
  SE    <- SE[index,index]
  d     <- d[index, index]
  gdist  <- gdist[index, index]

  return(list(SE=SE, d=d, gdist=gdist))

}

# @rdname struct_equiv
# @export
struct_equiv.matrix <- function(graph, v, inf.replace, groupvar, mode, ...) {

  # Running the algorithm
  if (length(groupvar)) {
    output <- struct_equiv_by(graph, v, inf.replace, groupvar, mode, ...)
  } else {

    geod <- igraph::distances(graph_from_adjacency_matrix(graph), mode=mode, ...)
    geod[!is.finite(geod)] <- inf.replace
    geod <- geod/max(geod, na.rm = TRUE)

    output <- struct_equiv_cpp(methods::as(geod, "dgCMatrix"), v)

    # Names
    rn <- rownames(graph)
    if (!length(rn)) rn <- 1:nrow(graph)
    output <- lapply(output, "dimnames<-", value=list(rn, rn))
  }

  return(output)
}

# @rdname struct_equiv
# @export
struct_equiv.dgCMatrix <- function(graph, v, inf.replace, groupvar, mode, ...) {

  if (length(groupvar)) {
    output <- struct_equiv_by(graph, v, inf.replace, groupvar, mode, ...)
  } else {
    # In order to use the SNA package functions, we need to coerce the graph
    # Into a -matrix.csc- object,
    geod <- igraph::distances(graph_from_adjacency_matrix(graph), mode=mode, ...)
    geod[!is.finite(geod)] <- inf.replace
    geod <- geod/max(geod, na.rm = TRUE)

    output <- struct_equiv_cpp(methods::as(geod, "dgCMatrix"), v)

    # Names
    rn <- rownames(graph)
    if (!length(rn)) rn <- 1:nrow(graph)
    output <- lapply(output, "dimnames<-", value=list(rn, rn))
  }
  return(output)
}

# @rdname struct_equiv
# @export
struct_equiv.array <- function(graph, v, inf.replace, groupvar, mode, ...) {
  t <- dim(graph)[3]
  output <- vector("list", t)
  for(i in 1:t) {
    output[[i]] <- struct_equiv.matrix(graph[,,i], v, inf.replace, groupvar, mode, ...)
  }

  # Naming
  tn <- dimnames(graph)[[3]]
  if (!length(tn)) tn <- 1:t

  names(output) <- tn

  output
}


# @rdname struct_equiv
# @export
struct_equiv.list <- function(graph, v, inf.replace, groupvar, mode, ...) {
  t <- length(graph)
  n <- nrow(graph[[1]])
  output <- vector("list", t)
  for(i in 1:t)
    output[[i]] <- struct_equiv.dgCMatrix(methods::as(graph[[i]], "dgCMatrix"),
                                          v, inf.replace, groupvar, mode,...)

  # Naming
  tn <- dimnames(graph)[[3]]
  if (!length(tn)) tn <- 1:t

  names(output) <- tn

  output
}
