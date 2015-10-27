#include <Rcpp.h>
#include "adjmat.hpp"

using namespace Rcpp;

// [[Rcpp::export]]
List struct_equiv_cpp(
    const NumericMatrix & z,
    double v = 1.0) {

  int n = z.ncol();
  NumericMatrix d(n,n);

  // Calculating Z vector as Z_i - Z_j = {z_ik - z_jk}
  NumericVector dmax(n, -1e5);
  for(int i=0;i<n;i++) {
    for(int j=0;j<i;j++) {
      // Cant be computed in the selfes (right?)
      if (i == j) continue;

      // Computing sum(z_ik - z_jk)
      double sumik = 0.0;
      double sumki = 0.0;
      for(int k=0;k<n;k++) {
        // Summation accross all but i and j
        if (k == i || k == j) continue;
        sumik += pow(z(i,k)-z(j,k), 2.0);
        sumki += pow(z(k,i)-z(k,j), 2.0);
      }

      // Adding up the results
      d(i,j) = pow(pow(z(i,j) - z(j,i), 2.0) + sumik + sumki, 0.5 );
      d(j,i) = d(i,j);
      if (dmax[i] < d(i,j)) dmax[i] = d(i,j);
    }
  }

  // Computing distances
  NumericMatrix SE(n,n);
  for(int i=0;i<n;i++) {

    for(int j=0;j<i;j++) {
      // Cant be computed in the selfes (right?)
      if (i == j) continue;

      // Computing sum(dmax - dkj)
      double sumdmaxd = 0.0;
      for(int k=0;k<n;k++) {
        if (k==i) continue;
        sumdmaxd += pow(dmax[i] - d(i,k), v);
      }

      SE(i,j) = pow(dmax[i] - d(i,j), v)/sumdmaxd ;
      SE(j,i) = SE(i,j);
    }
  }

  return List::create(_["SE"]=SE, _["d"]=d);
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

microbenchmark::microbenchmark(sna::sedist, struct_equiv)

# Use this implementation
z <- sna::geodist(adjmat, inf.replace = 0)
y <- struct_equiv(adjmat, 1)
y$d/1.118034
sna::sedist(adjmat, method = 'euclidean')
*/
