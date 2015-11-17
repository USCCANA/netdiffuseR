// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec seq_cpp(double from, double to, int lengthout) {
  arma::vec out(lengthout);
  double step = (to - from)/(lengthout - 1);
  for(int i=0;i<lengthout;i++)
    out(i) = from + i*step;

  return out;

}

// [[Rcpp::export]]
List grid_distribution(const arma::vec & x, const arma::vec & y, int n=100) {

  // Checking sizes of the vectors
  int m = x.size();
  int s = y.size();

  if (m!=s) stop("x and y don't have the same length.");

  // Crating empty matrix jointly with the sequences
  arma::mat distmat(n,n, arma::fill::zeros);
  double xlim[2], ylim[2];
  xlim[0] = x.min()-1e-10;
  xlim[1] = x.max()+1e-10;
  ylim[0] = y.min()-1e-10;
  ylim[1] = y.max()+1e-10;

  arma::vec xseq = seq_cpp(xlim[0], xlim[1], n + 1);
  arma::vec yseq = seq_cpp(ylim[0], ylim[1], n + 1);

  for(int k=0;k<m;k++)
    for(int i=0;i<n;i++) {
      bool cnt = false;
      for(int j=0;j<n;j++)
        // Testing if x and y are in the range
        if ( ((x(k) <= xseq(i+1)) & (x(k) > xseq(i))) & ((y(k) <= yseq(j+1)) & (y(k) > yseq(j)))) {
          distmat(i,j) += 1;
          cnt=true;
          break;
        }
      if (cnt) break;
    }

  // Output class mark
  NumericVector xmark(n);
  NumericVector ymark(n);

  for(int i=0;i<n;i++)
    xmark[i] = (xseq(i) + yseq(i+1))/2,
      ymark[i] = (yseq(i) + yseq(i+1))/2;

  return List::create(_["x"]=xmark, _["y"]=ymark, _["z"]=distmat);
}


// arma::mat grid_dist(const arma::vec & xran, const arma::vec & yran, const & arma::vec)

/***R
library(microbenchmark)
microbenchmark(
  seq_cpp(0,1,100),
  seq(0,1,length.out = 100),
  times=1000
)

# Testing
set.seed(123)
x <- rnorm(200000)
y <- rnorm(200000)
z <- grid_distribution(x,y,40)

sum(z$z)
with(z, contour(x,y,z/sum(z), nlevels = 40))
with(z, persp3d(as.vector(x),as.vector(y),z/sum(z), col="lightblue"))
*/

