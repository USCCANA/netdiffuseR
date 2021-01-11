euclidean_distance <- function(d) {

  ans <- (d - t(d)) ^ 2

  # Computing sum(z[ij] - z[ik])^2
  ids  <- t(utils::combn(1:nrow(d), 2)) # as.matrix(expand.grid(1:nrow(d), 1:nrow(d)))
  z_ik <- (d[ids[, 1], ] - d[ids[, 2], ])^2

  # Removing i,j
  z_ik[cbind(1:nrow(z_ik), ids[, 1])] <- 0
  z_ik[cbind(1:nrow(z_ik), ids[, 2])] <- 0
  ids <- cbind(ids, rowSums(z_ik))
  z_ik <- ans * 0
  z_ik[ids[,1:2]] <- ids[, 3]
  z_ik[ids[,2:1]] <- ids[, 3]

  # Computing sum(z[kj] - z[ki])^2
  z_ki <- (d[, ids[, 1]] - d[, ids[, 2]])^2

  # Removing i,j
  z_ki[cbind(ids[, 1], 1:ncol(z_ki))] <- 0
  z_ki[cbind(ids[, 2], 1:ncol(z_ki))] <- 0
  ids[, 3] <- colSums(z_ki)

  z_ki <- ans * 0
  z_ki[ids[,1:2]] <- ids[, 3]
  z_ki[ids[,2:1]] <- ids[, 3]

  # Final calculation
  sqrt(ans + z_ik + z_ki)

}

struct_equiv_new <- function(gdist, v = 1) {

  # Euclidean
  euclidean <- euclidean_distance(gdist)

  # Dmax
  n <- nrow(euclidean)
  d <- euclidean[cbind(1:n, max.col(euclidean))] #), nrow = n, ncol = n)

  # dmax(i) - d(ij)
  dmaxi_dij <- (d - euclidean) ^ v
  diag(dmaxi_dij) <- 0

  list(
    SE    = dmaxi_dij/(rowSums(dmaxi_dij) + 1e-15),
    d     = euclidean,
    gdist = gdist
  )

}

#' Structural Equivalence
#'
#' Computes structural equivalence between ego and alter in a network
#'
#' @template graph_template
#' @param v Numeric scalar. Cohesion constant (see details).
#' @param inf.replace Deprecated.
#' @param groupvar Either a character scalar (if \code{graph} is diffnet), or a
#' vector of size \eqn{n}.
#' @param ... Further arguments to be passed to \code{\link{approx_geodesic}}
#' (not valid for the print method).
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
#' \doi{10.1086/228667}
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
#' @author George G. Vega Yon & Thomas W. Valente
struct_equiv <- function(graph, v=1, inf.replace = 0, groupvar=NULL,  ...) {

  # Checking groupvar
  if ((length(groupvar)==1) && inherits(graph, "diffnet"))
    groupvar <- graph[[groupvar]]

  cls <- class(graph)

  output <- if ("matrix" %in% cls) {
      struct_equiv.dgCMatrix(methods::as(graph, "dgCMatrix"), v, inf.replace, groupvar,  ...)
    } else if ("dgCMatrix" %in% cls) {
      struct_equiv.dgCMatrix(graph, v, inf.replace, groupvar, ...)
    } else if ("array" %in% cls) {
      struct_equiv.list(apply(graph, 3, methods::as, Class="dgCMatrix"), v, inf.replace, groupvar, ...)
    } else if ("list" %in% cls) {
      struct_equiv.list(graph, v, inf.replace, groupvar, ...)
    } else if ("diffnet" %in% cls) {
      struct_equiv.list(graph$graph, v, inf.replace, groupvar, ...)
    } else
      stopifnot_graph(graph)


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
      "  - gdist : Geodesic distances matrix (n x n)",
      sep="\n"
      )

  return(invisible(x))
}


# Apply per grouping variable
struct_equiv_by <- function(graph, v, inf.replace, groupvar, ...) {
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
    out <- struct_equiv(subg, v, inf.replace, groupvar=NULL, ...)

    # Appending values
    index <- (N + 1L):(N + ng)
    SE[index,index]   <- out$SE
    d[index,index]    <- out$d
    gdist[index,index] <- out$gdist
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

#' Apply a function to a graph considering non-diagonal structural zeros
#'
#' When there are structural zeros given by groups, this function applies
#' a particular transformation function of a graph by groups returning a
#' square matrix of the same size of the original one with structural zeros
#' and the function applied by \code{INDICES}.
#'
#' @param graph A graph
#' @param INDICES A vector of length \eqn{n}.
#' @param fun A function. This function must return a matrix of class
#' \code{\link[Matrix:dgCMatrix-class]{dgCMatrix}} with the same dimension as
#' \code{dim(g)}.
#' @param ... Further arguments passed to \code{fun}
#'
#' @details The transformation function \code{fun} must return a square matrix
#' of size \eqn{m\times m}{m*m}, where \eqn{m} is the size of the subgroup
#' given by \code{INDICES}. See examples below
#'
#' @examples
#' # Rewiring a graph by community --------------------------------------------
#'
#' # Two Random graphs of different size
#' set.seed(123)
#' g0 <- rgraph_ba(m=2, self=FALSE)
#' g1 <- rgraph_ba(m=3, t=19, self=FALSE)
#'
#' # Need a place to store both networks together!
#' G <- methods::new(
#'   Class = "dgCMatrix",
#'   Dim   = c(1L,1L)*(nnodes(g0) + nnodes(g1)),
#'   p     = rep(0L, (nnodes(g0) + nnodes(g1)) + 1L)
#'   )
#'
#' # Filling the matrix
#' G[1:nnodes(g0),1:nnodes(g0)]                              <- g0
#' G[(nnodes(g0) + 1):nnodes(G), (nnodes(g0) + 1):nnodes(G)] <- g1
#'
#' # Creating an index (community)
#' indx <- c(rep(1, nnodes(g0)), rep(2, nnodes(g1)))
#'
#' # Apply the rewiring algorithm per group
#' ans <- transformGraphBy(G, indx, function(g, ...) {
#'   rewire_graph(g, 100, "swap")
#'   })
#'
#' ans
#'
#'
#' @export
transformGraphBy <- function(graph, INDICES, fun=function(g,...) g, ...)
  UseMethod("transformGraphBy")

#' @export
#' @rdname transformGraphBy
transformGraphBy.diffnet <- function(graph, INDICES, fun=function(g,...) g, ...) {
  for (per in as.character(graph$meta$pers))
    graph$graph[[per]] <- transformGraphBy(graph$graph[[per]], INDICES, fun, ...)

  return(graph)
}

#' @export
#' @rdname transformGraphBy
transformGraphBy.dgCMatrix <- function(graph, INDICES, fun=function(g,...) g, ...) {
  # Checking length of the grouping variable
  if (length(INDICES) != nvertices(graph))
    stop("The length of -INDICES-, ",length(INDICES),
         " elements, must be equal to the number of vertices in the graph, ",
         nvertices(graph)," vertices.")

  # Checking that is complete
  test <- which(!complete.cases(INDICES))
  if (length(test))
    stop("-INDICES- must not have incomplete cases. Check the following elements:\n\t",
         paste0(test, collapse=", "), ".")

  # Does the graph has ids?
  if (!length(rownames(graph))) {
    rn <- 1:nvertices(graph)
    dimnames(graph) <- list(rn, rn)
  }

  # Splitting the data and computing structural equivalence per case
  G    <- unique(INDICES)
  n    <- nvertices(graph)
  ans  <- methods::new("dgCMatrix", Dim=c(n,n), p=rep(0L,n+1L))

  rn <- NULL
  N  <- 0L
  for (g in G) {
    # Subsetting the graph
    test <- which(INDICES == g)
    subg <- graph[test, test, drop=FALSE]
    rn   <- c(rn, rownames(subg))
    ng   <- nvertices(subg)

    # Calculating indices
    i <- (N + 1L):(N + ng)

    # Computing
    ans[i,i] <- fun(subg,...)
    N <- N + ng
  }

  # Giving names
  dimnames(ans) <- list(rn, rn)

  # Ordening output
  i <- rownames(graph)
  return(ans[i,i])
}

# @rdname struct_equiv
# @export
struct_equiv.dgCMatrix <- function(graph, v, inf.replace, groupvar, ...) {

  if (length(groupvar)) {
    output <- struct_equiv_by(graph, v, inf.replace, groupvar, ...)
  } else {
    # In order to use the SNA package functions, we need to coerce the graph
    # Into a -matrix.csc- object,
    geod   <- do.call(approx_geodist, c(list(graph=graph), ...))
    geod@x <- geod@x/max(geod@x, na.rm = TRUE)

    output <- struct_equiv_new(geod, v)

    # Names
    rn <- rownames(graph)
    if (!length(rn)) rn <- 1:nrow(graph)
    output <- lapply(output, "dimnames<-", value=list(rn, rn))
  }
  return(output)
}

compare_matrix_and_vector <- function(m0, m1, v0, v1) {

  all((as.matrix(m0) - as.matrix(m1)) == 0) &
    all(v0 == v1)

}

# @rdname struct_equiv
# @export
struct_equiv.list <- function(graph, v, inf.replace, groupvar, ...) {
  t <- length(graph)
  n <- nrow(graph[[1]])
  output <- vector("list", t)

  # If groupvar is dynamic as well
  if (!is.list(groupvar)) {

    for(i in 1:t)
      output[[i]] <- struct_equiv.dgCMatrix(methods::as(graph[[i]], "dgCMatrix"),
                                            v, inf.replace, groupvar, ...)
  } else {

    for(i in 1:t)
      output[[i]] <- struct_equiv.dgCMatrix(methods::as(graph[[i]], "dgCMatrix"),
                                            v, inf.replace, groupvar[[i]], ...)
  }

  # Naming
  tn <- dimnames(graph)[[3]]
  if (!length(tn)) tn <- 1:t

  names(output) <- tn

  output
}
