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
rand_graph <- function(n=10, t=0, p=0.3, undirected=TRUE, weighted=FALSE, self=FALSE) {
  if (t==1) return(rand_graph_cpp(n, p, undirected, weighted, self))
  else return(rand_dyn_graph_cpp(n, t, p, undirected, weighted, self))
}
