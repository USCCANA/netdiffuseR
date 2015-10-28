#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

#include "adjmat.hpp"

using namespace Rcpp;

// [[Rcpp::export]]
List struct_equiv_cpp(
    const arma::mat & gdist,
    double v = 1.0,
    bool donly = false,
    bool inv = false, double invrep = 0.0) {

  int n = gdist.n_cols;
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
        sumik += pow(gdist(i,k)-gdist(j,k), 2.0);
        sumki += pow(gdist(k,i)-gdist(k,j), 2.0);
      }

      // Adding up the results
      d(i,j) = pow(pow(gdist(i,j) - gdist(j,i), 2.0) + sumik + sumki, 0.5 );

      // If only inverse required
      if (inv) {
        if (d(i,j)==0.0) d(i,j) = invrep;
        else d(i,j) = 1/d(i,j);
      }
      d(j,i) = d(i,j);
    }
  }

  // Computing distances
  NumericMatrix SE(n,n);

  // If only distance must be computed
  if (donly) return List::create(_["SE"]=SE, _["d"]=d, _["gdist"]=gdist);

  for(int i=0;i<n;i++) {

    // Getting the max of the line
    for(int j=0;j<n;j++) {
      if (i==j) continue;
      if (dmax[i] < d(i,j)) dmax[i] = d(i,j);
    }

    // Computing sum(dmax - dkj)
    double sumdmaxd = 0.0;
    for(int k=0;k<n;k++) {
      if (k==i) continue;
      sumdmaxd += pow(dmax[i] - d(k,i), v);
    }

    // Computing (dmax - d)/sum(dmax - d)
    for(int j=0;j<n;j++) {
      if (i==j) continue;
      SE(i,j) = pow(dmax[i] - d(j,i), v)/sumdmaxd;
    }
  }

  return List::create(_["SE"]=SE, _["d"]=d, _["gdist"]=gdist);
}

/***R
set.seed(1234)
n <- 10
edgelist <- matrix(sample(1:n, size = n, replace = TRUE), ncol=2)
adjmat <- edgelist_to_adjmat_cpp(edgelist)

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
