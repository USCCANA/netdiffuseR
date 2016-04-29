#' Erdos-Renyi model
#'
#' Generates a bernoulli random graph.
#'
#' @param n Integer. Number of vertices
#' @param t Integer. Number of time periods
#' @param p Double. Probability of a link between ego and alter.
#' @param undirected Logical scalar. Whether the graph is undirected or not.
#' @param weighted Logical. Whether the graph is weighted or not.
#' @param self Logical. Whether it includes self-edges.
#' @param as.edgelist Logical. When TRUE the graph is presented as an edgelist
#' instead of an adjacency matrix.
#' @details
#' For each pair of nodes \eqn{\{i,j\}}{{i,j}}, an edge is created
#' with probability \eqn{p}, this is, \eqn{Pr\{Link i-j\} = Pr\{x<p\}}{%
#' Pr{Link i-j}}, where \eqn{x} is drawn from a \eqn{Uniform(0,1)}.
#'
#' When \code{weighted=TRUE}, the strength of ties is given by
#' the random draw \eqn{x} used to compare against \eqn{p}, hence, if \eqn{x < p}
#' then the strength will be set to \eqn{x}.
#'
#' In the case of dynamic graphs, the algorithm is repeated \eqn{t} times, so the
#' networks are uncorrelated.
#' @references
#' Barabasi, Albert-Laszlo. "Network science book" Retrieved November 1 (2015)
#' \url{http://barabasi.com/book/network-science}.
#' @return A graph represented by an adjacency matrix (if \code{t=1}), or an array of
#' adjacency matrices (if \code{t>1}).
#' @export
#' @aliases bernoulli
#' @concept Bernoulli Random graph
#' @note The resulting adjacency matrix is store as a dense matrix, not as a
#' sparse matrix, hence the user should be careful when choosing the size of
#' the network.
#' @examples
#' \dontrun{
#' # Setting the seed
#' set.seed(123)
#'
#' # Generating an directed graph
#' rgraph_er(undirected=FALSE)
#'
#' # Comparing P(tie)
#' x <- rgraph_er(1000, p=.1)
#' sum(x)/length(x)
#'
#' # Several period random gram
#' rgraph_er(t=5)
#' }
#' @keywords distribution
#' @concept Erdos-Renyi random graph
#' @family simulation functions
#' @include graph_data.R
#' @author George G. Vega Yon
rgraph_er <- function(n=10, t=1, p=0.3, undirected=getOption("diffnet.undirected"), weighted=FALSE,
                       self=getOption("diffnet.self"), as.edgelist=FALSE) {

  # Generating the random graph
  if (t==1) graph <- rgraph_er_cpp(n, p, undirected, weighted, self)
  else graph <- rgraph_er_dyn_cpp(n, t, p, undirected, weighted, self)

  if (as.edgelist) {
    graph <- adjmat_to_edgelist(graph, undirected)
    attr(graph, "undirected") <- undirected
    return(graph)
  }

  # Naming dimensions
  if (t==1) dimnames(graph) <- list(1:n, 1:n)
  else {
    names(graph) <- 1:t
    for (i in 1:t)
      dimnames(graph[[i]]) <- list(1:n, 1:n)
  }

  attr(graph, "undirected") <- undirected
  return(graph)

  return(graph)
}

#' Barabasi-Albert model
#'
#' Generates a scale-free random graph based on Bollabas et al. (2001), also know as
#' \emph{Linearized Chord Diagram} (LCD) which has nice mathematical propoerties.
#'
#' @param m0 Integer scalar. Number of initial vertices in the graph.
#' @param m Integer scalar. Number of new edges per vertex added.
#' @param t Integer scalar. Number of time periods (steps).
#' @param graph Any class of accepted graph format (see \code{\link{netdiffuseR-graphs}}).
#' @return If \code{graph} is not provided, a static graph, otherwise an expanded
#' graph (\code{t} aditional vertices) of the same class as \code{graph}.
#'
#' The resulting graph will have \code{graph$meta$undirected = FALSE} if it is of
#' class \code{diffnet} and \code{attr(graph, "undirected")=FALSE} otherwise.
#' @family simulation functions
#' @aliases scale-free
#' @concept Scale-free random graph
#' @concept Barabasi-Albert model
#' @concept Random graph
#' @keywords distribution
#' @details
#' Based on Ballobás et al. (2001) creates a directed random graph of size
#' \code{t + m0}. A big difference with B-A model
#' is that this allows for loops (self/auto edges) and further multiple links,
#' nevertheless, as \eqn{t} increases, the number of such cases reduces.
#'
#' By default, the degree of the first \code{m0} vertices is set to be 2 (loops).
#' When \code{m>1}, as described in the paper, each new link from the new vertex
#' is added one at a time
#' \dQuote{counting \sQuote{outward half} of the edge being added as already contributing to the degrees}.
#'
#' @examples
#' # Using another graph as a base graph
#' graph <- rgraph_ba()
#' graph
#'
#' graph <- rgraph_ba(graph=graph)
#' @export
#' @references
#' Bollobás, B´., Riordan, O., Spencer, J., & Tusnády, G. (2001). The degree
#' sequence of a scale-free random graph process. Random Structures & Algorithms,
#' 18(3), 279–290. \url{http://doi.org/10.1002/rsa.1009}
#'
#' Albert-László Barabási, & Albert, R. (1999). Emergence of Scaling in Random
#' Networks. Science, 286(5439), 509–512. \url{http://doi.org/10.1126/science.286.5439.509}
#'
#' Albert-László Barabási. (2016). Network Science: (1st ed.). Cambridge University Press.
#' Retrieved from \url{http://barabasi.com/book/network-science}
#' @author George G. Vega Yon
rgraph_ba <- function(m0=1L, m=1L, t=10L, graph=NULL) {
  # When the graph is not null, then use it as a seed (starting point)
  if (length(graph)) {
    d <- dgr(graph)
    d[d==0] <- 1

    # Checking undirected (if exists)
    checkingUndirected(graph)

    # Parsing the class
    out <- switch(class(graph),
           matrix = rgraph_ba_cpp(methods::as(graph, "dgCMatrix"), d, m, t),
           dgCMatrix =rgraph_ba_cpp(graph, d, m, t),
           list = lapply(graph, function(x) rgraph_ba_cpp(x, dgr(x), m, t)),
           diffnet = lapply(graph$graph, function(x) rgraph_ba_cpp(x, dgr(x), m, t)),
           array = {
             out <- vector("list", dim(graph)[3])
             for (i in 1:dim(graph)[3])
               out[[i]] <- rgraph_ba_cpp(methods::as(graph[,,i], "dgCMatrix"),
                                         d[,i], m, t)
             out
           },
           stopifnot_graph(graph)
           )

    # BA is undirected by definition
    attr(out, "undirected") <- FALSE

    # If it is a static graph
    if (inherits(graph, c("dgCMatrix", "matrix"))) {

      n    <- nrow(graph)
      nnew <- nrow(out)

      rn <- rownames(graph)
      if (!length(rn)) rn <- 1:n

      ids <- c(rn, (n+1):nnew)

      dimnames(out) <- list(ids, ids)
      return(out)
    } else if (inherits(graph, "diffnet")) { # If it is a dynamic graph
      # BA is undirected by definition
      attr(out, "undirected") <- NULL
      graph$meta$undirected <- FALSE

      names(out) <- graph$meta$pers

      graph$graph <- out

      # TOA MAT
      n <- graph$meta$n
      nnew <- nrow(out[[1]])

      # Meta attributes
      graph$meta$n <- nnew
      graph$meta$ids <- c(graph$meta$ids, (n+1):nnew)

      graph$adopt <- rbind(
        graph$adopt,
        matrix(0, nrow=nnew-n, ncol=graph$meta$nper)
      )
      graph$cumadopt <- rbind(
        graph$cumadopt,
        matrix(0, nrow=nnew-n, ncol=graph$meta$nper)
      )

      graph$toa <- c(graph$toa, rep(NA, nnew-n))
      names(graph$toa) <- graph$meta$ids

      # Names
      dimnames(graph$adopt) <- list(graph$meta$ids, graph$meta$pers)
      dimnames(graph$cumadopt) <- list(graph$meta$ids, graph$meta$pers)

      for (i in 1:length(out))
        dimnames(graph$graph[[i]]) <- list(graph$meta$ids, graph$meta$ids)

      return(graph)

    } else if (inherits(graph, "list")) {

      tn <- names(graph)
      if (!length(tn)) tn <- 1:length(graph)
      names(out) <- tn

      n <- nrow(graph[[1]])
      nnew <- nrow(out[[1]])

      rn <- rownames(graph[[1]])
      if (!length(rn)) rn <- 1:n

      ids <- c(rn, (n+1):nnew)
      for (i in 1:length(out))
        dimnames(out[[i]]) <- list(ids, ids)

      return(out)
    } else { # In the case of an array
      n <- nrow(graph)
      nnew <- nrow(out[[1]])
      rn <- rownames(graph)

      if (!length(rn)) rn <- 1:n
      ids <- c(rn, (n+1):nnew)

      tn <- dimnames(graph)[[3]]
      if (!length(tn)) tn <- 1:dim(graph)[3]
      names(out) <- tn

      for (i in 1:length(out))
        dimnames(out[[i]]) <- list(ids, ids)

      return(out)
    }
  } else {
    out <- rgraph_ba_new_cpp(m0, m, t)
    ids <- 1:(m0+t)

    attr(out, "undirected") <- FALSE
    dimnames(out) <- list(ids, ids)
    return(out)
  }
}



#' Watts-Strogatz model
#'
#' Generates a small-world random graph.
#'
#' @param n Integer scalar. Set the size of the graph.
#' @param k Integer scalar. Set the initial degree of the ring (must be less than \eqn{n}).
#' @param p Numeric scalar/vector of length \eqn{T}. Set the probability of changing an edge.
#' @param both.ends Logical scalar. When \code{TRUE} rewires both ends.
#' @param self Logical scalar. When \code{TRUE}, allows loops (self edges).
#' @param multiple Logical scalar. When \code{TRUE} allows multiple edges.
#' @param undirected Logical scalar. Passed to \code{\link{ring_lattice}}
#' @return A random graph of size \eqn{n\times n}{n*n} following the small-world
#' model. The resulting graph will have \code{attr(graph, "undirected")=FALSE}.
#' @family simulation functions
#' @aliases small-world
#' @export
#' @details Implemented as in Watts and Strogatz (1998). Starts from an
#' undirected ring with \eqn{n} vertices all with degree \eqn{k} (so it must
#' be an even number), and then rewire each edge by setting the endpoint (so
#' now you treat it as a digraph) randomly any vertex in \eqn{N \setminus {i}}{N \\ {i}}
#' avoiding multiple links (by default) using the rewiring algorithm described on
#' the paper.
#'
#' @references
#' Watts, D. J., & Strogatz, S. H. (1998). Collective dynamics of "small-world"
#' networks. Nature, 393(6684), 440–2. \url{http://dx.doi.org/10.1038/30918}
#'
#' Newman, M. E. J. (2003). The Structure and Function of Complex Networks.
#' SIAM Review, 45(2), 167–256. \url{http://doi.org/10.1137/S003614450342480}
#' @author George G. Vega Yon
#' @examples
#'
#' library(igraph)
#' set.seed(7123)
#' x0 <- graph_from_adjacency_matrix(rgraph_ws(10,2, 0))
#' x1 <- graph_from_adjacency_matrix(rgraph_ws(10,2, .3))
#' x2 <- graph_from_adjacency_matrix(rgraph_ws(10,2, 1))
#'
#' oldpar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(x0, layout=layout_in_circle, edge.curved=TRUE, main="Regular")
#' plot(x1, layout=layout_in_circle, edge.curved=TRUE, main="Small-world")
#' plot(x2, layout=layout_in_circle, edge.curved=TRUE, main="Random")
#' par(oldpar)
#'
#' @include rewire.R
rgraph_ws <- function(n,k,p, both.ends=FALSE, self=FALSE, multiple=FALSE,
                      undirected=FALSE) {
  # out <- rewire_endpoints(ring_lattice(n, k, TRUE), p, both.ends,
  #                  self, multiple, undirected)
  out <- rewire_ws(ring_lattice(n, k, TRUE), k, p, self, multiple)

  # WS model is directed
  attr(out, "undirected") <- FALSE
  out
}


