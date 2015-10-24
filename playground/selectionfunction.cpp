#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
List select_egoalter_cpp(
  const NumericMatrix & adjmat_t0,
  const NumericMatrix & adjmat_t1,
  const NumericVector & adopt_t0,
  const NumericVector & adopt_t1
) {

  int n = adjmat_t0.ncol();
  IntegerMatrix change_mat(n, n);

  // Analazing the change from t-1 to t
  //  0: Stable
  //  1: Added
  // -1: Dropped
  for(int i=0;i<n;i++)
    for(int j=0;j<n;j++) {
      double chg = adjmat_t1(i,j) - adjmat_t0(i,j);
      if      (chg > 0) change_mat(i,j) = 1;
      else if (chg < 0) change_mat(i,j) = -1;
    }

  // Classifies dinamics between 1 and 16 depending on whether alter and ego
  // changed behavior between t and t-1. Follows classification on Valente
  //    n     y
  //    n  y  n  y
  //n n 1  2  9 10
  //  y 3  4 11 12
  //y n 5  6 13 14
  //  y 7  8 15 16
  IntegerMatrix select_mat_a(n,16); // Added
  IntegerMatrix select_mat_d(n,16); // Dropped
  IntegerMatrix select_mat_s(n,16); // Stable

  int cat = 1;

  for(int i=0;i<n;i++)
    for(int j=0;j<n;j++) {
      if (i==j) continue;
      // Fitting the category
      cat =
        (adopt_t0[i] > 0)*4 + (adopt_t1[i] > 0)*2 +
        (adopt_t0[j] > 0)*8 + (adopt_t1[j] > 0)*1 + 1;

      if      (change_mat(i,j) > 0) select_mat_a(i, cat - 1) += 1;
      else if (change_mat(i,j) < 0) select_mat_d(i, cat - 1) += 1;
      else if (adjmat_t1(i,j) != 0) select_mat_s(i, cat - 1) += 1;
    }

  return List::create(
    _["select_a"] = select_mat_a,
    _["select_d"] = select_mat_d,
    _["select_s"] = select_mat_s);
}





/***R
Rcpp::sourceCpp("/home/george/Documents/usc/software/diffusiontest/playground/adjmat.cpp")

load("~/Documents/usc/software/diffusiontest/data/trade.rda")
load("~/Documents/usc/software/diffusiontest/data/ratify.rda")

adjmat <- with(trade, edgelist_to_adjmat(
  cbind(id,alter),
  times=year
))

ratByCon <- split(ratify, ratify$treaty)
adoptslist<-lapply(ratByCon, function(x) adoptMat(x$year))

out <- selectionFunctionEgoAlter(
  adjmat, adoptslist[[1]]$Adopt_mat, 2
)

select_egoalter <- function(adjmat_t0, adjmat_t1, adopt_t0, adopt_t1) {
  out <- select_egoalter_cpp(adjmat_t0, adjmat_t1, adopt_t0, adopt_t1)
  out <- do.call(cbind, out)
  colnames(out) <- c(
    sprintf('select_a_%02d', 1:16), sprintf('select_d_%02d', 1:16),
    sprintf('select_s_%02d', 1:16))
  cbind(time=1, id=1:nrow(out), out)
}

adopt <- adoptslist[[1]]$Adopt_mat
out2 <- select_egoalter(adjmat[,,1], adjmat[,,2], adopt[,1], adopt[,2])

# Comparing the 2 outputs
max(out[,-(1:2)] - out2[,-(1:2)])

library(microbenchmark)
microbenchmark(
  old = selectionFunctionEgoAlter(adjmat, adopt, 2),
  new = select_egoalter(adjmat[,,1], adjmat[,,2], adopt[,1], adopt[,2]), times=1
)

*/
