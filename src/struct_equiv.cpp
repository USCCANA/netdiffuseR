// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
List struct_equiv_cpp(
    const arma::sp_mat & graph, // Must be a geodesic distances graph
    double v = 1.0,
    bool unscaled = false,
    bool inv = false, double invrep = 0.0) {

  int n = graph.n_cols;
  NumericMatrix d(n,n);

  // Calculating Z vector as Z_i - Z_j = {z_ik - z_jk}
  NumericVector dmax(n, -1e5);
  for(int i=0;i<n;i++) {
    for(int j=0;j<i;j++) {

      // Computing sum(z_ik - z_jk)
      double sumik = 0.0;
      double sumki = 0.0;
      for(int k=0;k<n;k++) {
        // Summation accross all but i and j
        if (k == i || k == j) continue;
        sumik += pow(graph.at(i,k)-graph.at(j,k), 2.0);
        sumki += pow(graph.at(k,i)-graph.at(k,j), 2.0);
      }

      // Adding up the results
      d.at(i,j) = pow(pow(graph.at(i,j) - graph.at(j,i), 2.0) + sumik + sumki, 0.5 );

      // If only inverse required
      if (inv && unscaled) d.at(i,j) = 1/(d.at(i,j) + 1e-10);

      d.at(j,i) = d.at(i,j);
    }
  }

  // If only distance must be computed
  if (unscaled) return List::create(_["SE"]=d, _["d"]=d, _["gdist"]=graph);

  // Computing distances
  NumericMatrix SE(n,n);

  for(int i=0;i<n;i++) {

    // Getting the max of the line
    for(int j=0;j<n;j++) {
      if (i==j) continue;
      if (dmax[i] < d.at(i,j)) dmax[i] = d.at(i,j);
    }

    // Computing sum(dmax - dkj)
    double sumdmaxd = 0.0;
    for(int k=0;k<n;k++) {
      if (k==i) continue;
      sumdmaxd += pow(dmax[i] - d.at(k,i), v);
    }

    // Computing (dmax - d)/sum(dmax - d)
    for(int j=0;j<n;j++) {
      if (i==j) continue;
      SE.at(i,j) = pow(dmax[i] - d.at(j,i), v)/(sumdmaxd + 1e-10);
    }

    // If inverse required
    if (inv) {
      for(int j=0;j<n;j++) {
        SE.at(i,j) = 1/(SE.at(i,j) + 1e-10);
      }
    }
  }

  return List::create(_["SE"]=SE, _["d"]=d, _["gdist"]=graph);
}

/** *R
set.seed(1234)
adjmat <- rand_graph_cpp()

struct_equiv <- function(adjmat, v=1, ...) {
  geod <- sna::geodist(adjmat, inf.replace = 0, ...)
  geod[["gdist"]] <- geod[["gdist"]]/max(geod[["gdist"]])
  struct_equiv_cpp(geod[["gdist"]], v)
}

microbenchmark::microbenchmark(sna::sedist, struct_equiv, times=10000)

# Use this implementation
new <- struct_equiv(adjmat)
old <- sna::sedist(adjmat, method = 'euclidean')
rowSums(new$SE)
new
old
*/
