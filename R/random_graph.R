#' Erdos-Reyi model
#'
#' Several random graphs algorithms.
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
#' \url{http://barabasi.com/networksciencebook/}.
#' @return A graph represented by an adjacency matrix (if t=1), or an array of
#' adjacency matrices (if t>1).
#' @export
#' @concept Random graph
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
rgraph_er <- function(n=10, t=1, p=0.3, undirected=getOption("diffnet.undirected"), weighted=FALSE,
                       self=getOption("diffnet.self"), as.edgelist=FALSE) {

  # Generating the random graph
  if (t==1) graph <- rgraph_er_cpp(n, p, undirected, weighted, self)
  else graph <- rgraph_er_dyn_cpp(n, t, p, undirected, weighted, self)

  if (as.edgelist) return(adjmat_to_edgelist(graph, undirected))

  # Naming dimensions
  if (t==1) dimnames(graph) <- list(1:n, 1:n)
  else {
    names(graph) <- 1:t
    for (i in 1:t)
      dimnames(graph[[i]]) <- list(1:n, 1:n)
  }

  return(graph)
}

#' Barabasi-Albert model
#' @param m0 Integer scalar. Number of initial vertices in the graph.
#' @param m Integer scalar. Number of new edges per vertex added.
#' @param t Integer scalar. Number of time periods (steps).
#' @param graph Any class of accepted graph format (see \code{\link{netdiffuseR-graphs}}).
#' @return A graph.
#' @family simulation functions
#' @concept Scale-free random graph
#' @concept Barabasi-Albert model
#' @concept Random graph
#' @keywords distribution
#' @details
#' Creates an undirected random graph of size \code{t + m0}.
#' @examples
#' # Using another graph as a base graph
#' graph <- rgraph_ba()
#' graph
#'
#' graph <- rgraph_ba(graph=graph)
#' @export
rgraph_ba <- function(m0=1L, m=1L, t=10L, graph=NULL) {
  # When the graph is not null, then use it as a seed (starting point)
  if (length(graph)) {
    d <- dgr(graph)
    d[d==0] <- 1
    rgraph_ba_cpp(graph, d, m, t)
  }
  else rgraph_ba_new_cpp(m0, m, t)
}

#' Watts-Strogatz model
#' @param n Integer scalar. Set the size of the graph.
#' @param k Integer scalar. Set the initial degree of the ring (must be less than \eqn{n}).
#' @param p Numeric scalar. Set the probability of changing an edge.
#' @param both.ends Logical scalar. When \code{TRUE} rewires both ends.
#' @param self Logical scalar. When \code{TRUE}, allows loops (self edges).
#' @param multiple Logical scalar. When \code{TRUE} allows multiple edges.
#' @return A random graph of size \eqn{n\times n}{n*n} following the small-world
#' model.
#' @export
rgraph_ws <- function(n,k,p, both.ends=FALSE, self=FALSE, multiple=FALSE) {
  rewire_graph_cpp(ring_lattice(n, k), p, both.ends,
                   self, multiple, TRUE)
}

#' Rewires a graph
#' @inheritParams rgraph_ws
#' @param graph Any class of accepted graph format (see \code{\link{netdiffuseR-graphs}})
#' @export
rewire_graph <- function(graph, p, both.ends=FALSE, self=FALSE, multiple=FALSE,
                         undirected=getOption("diffnet.undirected")) {
  out <- switch(class(graph),
    dgCMatrix = rewire_graph.dgCMatrix(graph, p, both.ends, self, multiple, undirected),
    list = rewire_graph.list(graph, p, both.ends, self, multiple, undirected),
    matrix = rewire_graph.dgCMatrix(
      methods::as(graph, "dgCMatrix"), p, both.ends, self, multiple, undirected),
    diffnet = rewire_graph.list(graph$graph, p, both.ends, self, multiple,
                                graph$meta$undirected),
    array = rewire_graph.array(graph, p, both.ends, self, multiple, undirected),
    stopifnot_graph(graph)
  )

  # If diffnet, then it must return the same object but rewired
  if (inherits(graph, "diffnet")) {
    graph$graph <- out
    return(graph)
  }

  return(out)
}

# @rdname rewire_graph
rewire_graph.list <- function(graph, p, both.ends=FALSE, self=FALSE, multiple=FALSE,
                              undirected=getOption("diffnet.undirected")) {
  t   <- length(graph)
  out <- vector("list", t)

  # Names
  names(ngraph) <- names(graph)

  for (i in 1:t) {
    out[[i]] <- rewire_graph_cpp(graph[[i]], p, both.ends, self, multiple,
                                    undirected)
    # Names
    dimnames(out[[i]]) <- dimnames(graph[[i]])
  }

  out
}

# @rdname rewire_graph
rewire_graph.dgCMatrix <- function(graph, p, both.ends=FALSE, self=FALSE, multiple=FALSE,
                         undirected=getOption("diffnet.undirected")) {
  out <- rewire_graph_cpp(graph, p, both.ends, self, multiple, undirected)
  dimnames(out) <- dimnames(graph)
}

# @rdname rewire_graph
rewire_graph.array <-function(graph, p, both.ends=FALSE, self=FALSE, multiple=FALSE,
                              undirected=getOption("diffnet.undirected")) {
  n   <- dim(graph)[1]
  t   <- dim(graph)[3]
  out <- vector("list", t)

  # Checking time names
  tn <- dimnames(graph)[[3]]
  if (!length(tn)) tn <- 1:t
  names(out) <- tn

  # Rewiring
  for(i in 1:t) {
    out[[i]] <- rewire_graph_cpp(
      methods::as(graph[,,i], "dgCMatrix"), p, both.ends, self, multiple, undirected)
    dimnames(out[[i]]) <- dimnames(graph[,,i])
  }

  out
}
