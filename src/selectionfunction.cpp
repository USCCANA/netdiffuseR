// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame select_egoalter_cpp(
  const arma::sp_mat & adjmat_t0,
  const arma::sp_mat & adjmat_t1,
  const NumericVector & adopt_t0,
  const NumericVector & adopt_t1
) {

  int n = adjmat_t0.n_cols;
  arma::sp_mat change_mat(n, n);

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

      // Observe that indexes are from 0, that's why we have to substract 1 from
      // the categories
      if      (change_mat(i,j) > 0) select_mat_a(i, cat - 1) += 1;
      else if (change_mat(i,j) < 0) select_mat_d(i, cat - 1) += 1;
      else if (adjmat_t1(i,j) != 0) select_mat_s(i, cat - 1) += 1;
    }

  return DataFrame::create(
    _["select_a"] = select_mat_a,
    _["select_d"] = select_mat_d,
    _["select_s"] = select_mat_s);
}


/***R
source("/home/george/Documents/usc/software/diffusiontest/playground/adjmat.R", echo = FALSE)

library(microbenchmark)
library(diffusiontest)

load("~/Dropbox/usc/software/diffusiontest/data/trade.rda")
load("~/Documents/usc/software/diffusiontest/data/ratify.rda")

adjmat <- with(trade, edgelist_to_adjmat(
  cbind(id,alter),
  times=year
))

ratByCon <- split(ratify, ratify$treaty)
adoptslist<-lapply(ratByCon, function(x) adoptMat(x$year))

out <- selectionFunctionEgoAlter(
  adjmat, adoptslist[[1]]$Adopt_mat, 3
)

select_egoalter <- function(...) UseMethod("select_egoalter")

select_egoalter.array <- function(adjmat, adopt, period=NULL) {

  # Computing selection mat and coersing into a single matrix
  nper <- dim(adjmat)[3]
  if (length(period) && !(period %in% 2:nper))
    stop('Invalid period selected. Should be between -', 2,'- and -', nper,'-')

  # Creating column names
  cn <- c('time', 'id',
    sprintf('select_a_%02d', 1:16), sprintf('select_d_%02d', 1:16),
    sprintf('select_s_%02d', 1:16))

  # Output parameters
  n   <- dim(adjmat)[1]
  ids <- 1:n

  if (length(period)) {
    out <- cbind(per=period, ids, do.call(cbind, select_egoalter_cpp(
      adjmat[,,period-1], adjmat[,,period],
      adopt[,period-1], adopt[,period])))

    # Assigning colnames
    colnames(out) <- cn

    return(out)
  }

  # Looping over periods
  out <- vector("list", nper-1)
  for (i in 2:nper) {
    out[[i-1]] <- cbind(time=i, ids, do.call(cbind, select_egoalter_cpp(
      adjmat[,,i-1], adjmat[,,i],adopt[,i-1], adopt[,i])))

    # Assigning colnames
    colnames(out[[i-1]]) <- cn
  }

  return(array(unlist(out), dim=c(n,length(cn),nper)))
}

adopt <- adoptslist[[1]]$Adopt_mat
out2 <- select_egoalter(adjmat, adopt, period = 3)
out3 <- select_egoalter(adjmat, adopt)

# Comparing the 2 outputs
max(out[,-(1:2)] - out2[,-(1:2)])

library(microbenchmark)
microbenchmark(
  old = selectionFunctionEgoAlter(adjmat, adopt, 2),
  new = select_egoalter(adjmat, adopt, 2), times=5
)

# Unit: milliseconds
#   expr         min          lq        mean      median         uq         max neval cld
#     old 2134.850649 2202.994444 2464.911493 2254.061699 2464.91809 4090.787622   100   b
#     new    1.800653    1.877095    2.106013    1.944777    2.02516    3.591877   100  a
# Around 1000 times faster!
*/
