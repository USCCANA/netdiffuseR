#' Random graphs
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
#' \code{rgraph_er} Creates an Erdos-Renyi graph
#'
#' \code{rgraph_ba} Creates a Barabasi-Albert graph
#' @section Erdos-Renyi (bernoulli):
#' Basically, for each pair of nodes \eqn{\{i,j\}}{{i,j}}, an edge is created
#' with probability \eqn{p}, this is, \eqn{Pr\{Link i-j\} = Pr\{x<p\}}{%
#' Pr{Link i-j}}, where \eqn{x} is drawn from a \eqn{Uniform(0,1)}.
#'
#' When \code{weighted=TRUE}, the strength of ties is given by
#' the random draw \eqn{x} used to compare against \eqn{p}, hence, if \eqn{x < p}
#' then the strength will be set to \eqn{x}.
#'
#' In the case of dynamic graphs, the algorithm is repeated \eqn{t} times, so the
#' networks are uncorrelated.
#' @section Barabasi-Albert:
#' Creates an undirected random graph of size \code{t + m0}
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
#' @concept Scale-free random graph
#' @concept Barabasi-Albert model
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
