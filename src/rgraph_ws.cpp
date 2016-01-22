// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

arma::umat sparse_indexes(const arma::sp_mat & mat) {

  int n = mat.n_nonzero;
  arma::umat indices(n,2);
  int curcol = (int) mat.col_ptrs[0]/mat.n_rows;

  // If the matrix is empty (which makes no sense)
  if (n == curcol) return indices;

  int j = 0;
  int cumrow = 0;
  for (int i=0;i<n;i++) {
    // Figuring out what column
    if (cumrow >= mat.col_ptrs[j+1]) curcol = mat.col_ptrs[++j];

    // Asigning indexes
    indices.at(i,0) = mat.row_indices[i];
    indices.at(i,1) = j;
    ++cumrow;
  }
  // return indices;
  return indices;
}


//' Ring lattice graph
//'
//' Creates a ring lattice with \eqn{n} vertices, each one of degree (at most) \eqn{k}
//' as an undirected graph. This is the basis of \code{\link{rgraph_ws}}.
//' @param n Integer scalar. Size of the graph.
//' @param k Integer scalar. Degree of each vertex.
//' @details Since the created graph is undirected, the degree of each node always
//' even. So if \code{k=3}, then the degree will be \code{2}.
//' @return A sparse matrix of class \code{\link[Matrix:dgCMatrix-class]{dgCMatrix}} of size
//' \eqn{n\times n}{n * n}.
//' @references Watts, D. J., & Strogatz, S. H. (1998). Collective dynamics of
//' “small-world” networks. Nature, 393(6684), 440–2. \url{http://doi.org/10.1038/30918}
//' @export
// [[Rcpp::export]]
arma::sp_mat ring_lattice(int n, int k) {

  if ((n-1) < k)
    stop("k can be at most n - 1");

  arma::sp_mat graph(n,n);

  // Connecting to k next & previous neighbour
  for (int i=0;i<n;i++) {
    for (int j=1;j<=k/2;j++) {
      // Next neighbor
      int l = i+j;
      if (l >= n) l = l - n;
      graph.at(i,l) += 1.0;

      // Previous neighbor
      l = i-j;
      if (l < 0) l = l + n;
      graph.at(i,l) += 1.0;
      // graph.aAbstract Submission is now closed. Thank you!t(l,i) += 1.0;
    }
  }
  return graph;
}

/***R
library(Matrix)
library(sna)
x <- ring_lattice(6,2)
x
gplot(as.matrix(x), displaylabels = TRUE, mode="circle", jitter = FALSE)

# x <- ring_cpp(10,6)
# x
# gplot(as.matrix(x), displaylabels = TRUE, mode="circle", jitter = FALSE)
#
# x <- ring_cpp(7,7)
# x
# gplot(as.matrix(x), displaylabels = TRUE, mode="circle", jitter = FALSE)
*/

// [[Rcpp::export]]
arma::sp_mat rewire_graph_cpp(
    const arma::sp_mat & graph, double p,
    bool both_ends=false,
    bool self=false, bool multiple=false,
    bool undirected=false) {

  // Clonning graph
  arma::sp_mat newgraph(graph);
  int n = graph.n_cols;

  // Getting the indexes
  arma::umat indexes = sparse_indexes(graph);

  for (int i= 0;i<indexes.n_rows; i++) {

    // Checking user interrupt
    if (i % 1000 == 0)
      Rcpp::checkUserInterrupt();

    if (unif_rand() > p) continue;

    // Indexes
    int j = indexes.at(i, 0);
    int k = indexes.at(i, 1);

    // Since it is undirected...
    if (undirected && (j > k) ) continue;

    // New end(s)
    int newk, newj;

    // If rewiring is in both sides, then we start from the left
    if (both_ends) newj = floor( (n-1)*unif_rand() );
    else newj = j;

    // Checking for conditions
    int wcount = 0;
    bool conditions = true;
    arma::uvec picked(n, arma::fill::zeros);

    while (conditions) {
      newk = floor( (n-1)*unif_rand() );

      // In the case that the individual actually is connected to everyone and
      // multiple is not allowed, this is needed to break out the loop
      if (picked.at(newk)) {
        if (++wcount >= n*2) break;
        continue;
      } else picked.at(newk) = 1;

      if (undirected && newj > newk) continue;

      if (!self && newj == newk) continue;
      if (!multiple && (graph.at(newj, newk) != 0)) continue;

      conditions = false;
    }

    // Setting zeros
    double w = newgraph.at(j,k);
    newgraph.at(j,k) = 0;
    if (undirected) newgraph.at(k,j) = 0;

    // Adding up
    newgraph.at(newj,newk) += w;
    if (undirected) newgraph.at(newk,newj) += w;

  }
  return newgraph;
}


/***R
rgraph_ws <- function(n,k,p, both_ends=FALSE, self=FALSE, multiple=FALSE) {
  rewire_graph_cpp(ring_lattice(n, k), p, both_ends,
                   self, multiple, true)
}

x <- ring_lattice(10, 2)
gplot(as.matrix(x), mode="circle", jitter=FALSE, usecurve = TRUE, gmode = "graph")
gplot(as.matrix(rewire_graph_cpp(x, .1)), mode="circle", jitter=FALSE, usecurve = TRUE, gmode = "graph")
gplot(as.matrix(rewire_graph_cpp(x, 1)), mode="circle", jitter=FALSE, usecurve = TRUE, gmode = "graph")
*/
