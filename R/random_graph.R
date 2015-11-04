#' Erdos-Renyi (bernoulli) random graph
#'
#' Follows the G(N,p) model
#'
#' @param n Integer. Number of vertices
#' @param t Integer. Number of time periods
#' @param p Double. Probability of connection between ego and alter.
#' @param undirected Logical scalar. Whether the graph is undirected or not.
#' @param weighted LoWhether the graph is weighted or not.
#' @param self Wheter it includes self-edges.
#' @param as.edgelist Logical. When TRUE the graph is presented as an edgelist
#' instead.
#' @details
#'
#' Basically, for each pair of nodes \{i,j\}, an edge is created with probability
#' p, this is, P\{Link i-j\} = P\{x<p\}, where x is drawn from a Uniform(0,1).
#'
#' When \code{weighted=TRUE}, the strength of ties is given by
#' the random draw x used to compare against p, hence, if x < p then the strength
#' will be set to x.
#'
#' In the case of dynamic graphs, the algorithm is repeated t times.
#'
#' @references
#' Barabási, Albert-László. "Network science book" Retrieved November 1 (2015)
#' \url{http://barabasi.com/networksciencebook/}.
#' @return A graph represented by an adjacency matrix (if t=1), or an array of
#' adjacency matrices (if t>1).
#' @export
#' @note The resulting adjacency matrix is dense (hence, be careful with the size)
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
rand_graph <- function(n=10, t=1, p=0.3, undirected=TRUE, weighted=FALSE,
                       self=FALSE, as.edgelist=FALSE) {

  # Generating the random graph
  if (t==1) graph <- rand_graph_cpp(n, p, undirected, weighted, self)
  else graph <- rand_dyn_graph_cpp(n, t, p, undirected, weighted, self)

  if (as.edgelist) return(adjmat_to_edgelist(graph, undirected))
  else return(graph)
}
