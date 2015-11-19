#' Erdos-Renyi (bernoulli) random graph
#'
#' Using the \eqn{G(N,p)} model, generates a random graph.
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
#'
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
#'
#' @references
#' Barabasi, Albert-Laszlo. "Network science book" Retrieved November 1 (2015)
#' \url{http://barabasi.com/networksciencebook/}.
#' @return A graph represented by an adjacency matrix (if t=1), or an array of
#' adjacency matrices (if t>1).
#' @export
#' @note The resulting adjacency matrix is store as a dense matrix, not as a
#' sparse matrix, hence the user should be careful when choosing the size of
#' the network.
#' @examples
#' \dontrun{
#' # Setting the seed
#' set.seed(123)
#'
#' # Generating an directed graph
#' rand_graph(undirected=FALSE)
#'
#' # Comparing P(tie)
#' x <- rand_graph(1000, p=.1)
#' sum(x)/length(x)
#'
#' # Several period random gram
#' rand_graph(t=5)
#' }
#' @keywords distribution
rand_graph <- function(n=10, t=1, p=0.3, undirected=TRUE, weighted=FALSE,
                       self=FALSE, as.edgelist=FALSE) {

  # Generating the random graph
  if (t==1) graph <- rand_graph_cpp(n, p, undirected, weighted, self)
  else graph <- rand_dyn_graph_cpp(n, t, p, undirected, weighted, self)

  if (as.edgelist) return(adjmat_to_edgelist(graph, undirected))
  else return(graph)
}
