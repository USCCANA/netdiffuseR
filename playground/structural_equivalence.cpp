#include <Rcpp.h>
#include "adjmat.hpp"

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix struct_equiv(
    const NumericMatrix & z,
    const IntegerMatrix & use) {
  int n = z.ncol();
  NumericMatrix d(n,n);

  double dmax = -1e5;
  for(int i=0;i<n;i++)
    for(int j=0;j<n;j++) {
      // Computing sum(z_ik - z_jk)
      double sumiz = 0.0;
      double sumzi = 0.0;
      for(int k=0;k<n;k++) {
        // Only use distance when is accountable
        if (use(i,k) && use(j,k)) sumiz += pow(z(i,k)-z(j,k), 2.0);
        if (use(k,i) && use(k,j)) sumzi += pow(z(k,i)-z(k,j), 2.0);
      }

      // Adding up the results
      d(i,j) = pow( sumiz + sumzi, 0.5 );
      d(j,i) = d(i,j);
      if (dmax < d(i,j)) dmax = d(i,j);
    }

  NumericMatrix SE(n,n);
  for(int i=0;i<n;i++)
    for(int j=0;j<i;j++) {

      // Computing sum(dmax - dkj)
      double sumdmaxd = 0.0;
      for(int k=0;k<n;k++)
        sumdmaxd += dmax - d(k,j);

      SE(i,j) = d(i,j) /*(dmax - d(i,j))/sumdmaxd*/;
      SE(j,i) = SE(i,j);
    }

  return SE;
}

/***R
set.seed(1234)
n <- 10
edgelist <- matrix(sample(1:n, size = n, replace = TRUE), ncol=2)
adjmat <- edgelist_to_adjmat_cpp(edgelist)

# Use this implementation
z <- sna::geodist(adjmat)
sna::sedist(adjmat, method = 'euclidean')
round(struct_equiv(z$gdist, z$counts), 6)
*/
