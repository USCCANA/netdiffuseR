#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
List select_egoalter(
  const NumericMatrix & adjmat_t0,
  const NumericMatrix & adjmat_t1,
  const NumericVector & adopt_t0,
  const NumericVector & adopt_t1
) {

  int n = adjmat_t0.ncol();
  NumericMatrix change_mat(n,n);

  // Analazing the change from t-1 to t
  //  0: Stable
  //  1: Added
  // -1: Dropped
  for(int i=0;i<n;i++)
    for(int j=0;j<n;j++) {
      double chg = adjmat_t1(i,j) - adjmat_t0(i,j);
      if      (chg > 0) change_mat(i,j) = 1.0;
      else if (chg < 0) change_mat(i,j) = -1.0;
    }

  // Classifies dinamics between 1 and 16 depending on whether alter and ego
  // changed behavior between t and t-1. Follows classification on Valente
  //    n     y
  //    n  y  n  y
  //n n 1  2  9 10
  //  y 3  4 11 12
  //y n 5  6 13 14
  //  y 7  8 15 16
  NumericMatrix select_mat_a(n,16); // Added
  NumericMatrix select_mat_d(n,16); // Dropped
  NumericMatrix select_mat_s(n,16); // Stable

  int cat = 1;

  for(int i=0;i<n;i++)
    for(int j=0;j<n;j++) {
      if (i==j) continue;
      // Fitting the category
      cat =
        (adopt_t0[i] > 0)*4 + (adopt_t1[i] > 0)*2 +
        (adopt_t0[j] > 0)*8 + (adopt_t0[j] > 0)*1 + 1;

      if (!change_mat(i,j)) select_mat_s(i, cat) += 1;
      else if (change_mat(i,j)<0) select_mat_d(i, cat) += 1;
      else select_mat_a(i, cat) += 1;
    }

  return List::create(
    _["select_a"] = select_mat_a,
    _["select_d"] = select_mat_d,
    _["select_s"] = select_mat_s);
}

/* selectionFunctionEgoAlter <- function(all_nets, Adopt_mat, Time){
  n <- nrow(all_nets)
  if (Time>1) {
    change_mat   <- (all_nets[,,(Time-1)] + all_nets[,,(Time)])
    change_mat_s <- (all_nets[,,(Time-1)] - all_nets[,,(Time)])
    net_stabl  <- matrix(as.integer(change_mat ==2), n, n)
    net_added  <- matrix(as.integer(change_mat_s ==-1), n, n)
    net_dropd  <- matrix(as.integer(change_mat_s ==1), n, n)

    adopts <- Adopt_mat[,Time]
    adoptstm1 <- Adopt_mat[,(Time-1)]
    select_mat <- matrix(0, n, n)

    for (i in 1:n) {
      for (j in 1:n) {
        if (i==j) {next}

        if ((adoptstm1[i]==0) & (adopts[i]==0) & (adoptstm1[j]==0) & (adopts[j]==0)) {select_mat[i,j] <-1}
        if ((adoptstm1[i]==0) & (adopts[i]==0) & (adoptstm1[j]==0) & (adopts[j]==1)) {select_mat[i,j] <-2}
        if ((adoptstm1[i]==0) & (adopts[i]==1) & (adoptstm1[j]==0) & (adopts[j]==0)) {select_mat[i,j] <-3}
        if ((adoptstm1[i]==0) & (adopts[i]==1) & (adoptstm1[j]==0) & (adopts[j]==1)) {select_mat[i,j] <-4}
        if ((adoptstm1[i]==1) & (adopts[i]==0) & (adoptstm1[j]==0) & (adopts[j]==0)) {select_mat[i,j] <-5}
        if ((adoptstm1[i]==1) & (adopts[i]==0) & (adoptstm1[j]==0) & (adopts[j]==1)) {select_mat[i,j] <-6}
        if ((adoptstm1[i]==1) & (adopts[i]==1) & (adoptstm1[j]==0) & (adopts[j]==0)) {select_mat[i,j] <-7}
        if ((adoptstm1[i]==1) & (adopts[i]==1) & (adoptstm1[j]==0) & (adopts[j]==1)) {select_mat[i,j] <-8}
        if ((adoptstm1[i]==0) & (adopts[i]==0) & (adoptstm1[j]==1) & (adopts[j]==0)) {select_mat[i,j] <-9}
        if ((adoptstm1[i]==0) & (adopts[i]==0) & (adoptstm1[j]==1) & (adopts[j]==1)) {select_mat[i,j] <-10}
        if ((adoptstm1[i]==0) & (adopts[i]==1) & (adoptstm1[j]==1) & (adopts[j]==0)) {select_mat[i,j] <-11}
        if ((adoptstm1[i]==0) & (adopts[i]==1) & (adoptstm1[j]==1) & (adopts[j]==1)) {select_mat[i,j] <-12}
        if ((adoptstm1[i]==1) & (adopts[i]==0) & (adoptstm1[j]==1) & (adopts[j]==0)) {select_mat[i,j] <-13}
        if ((adoptstm1[i]==1) & (adopts[i]==0) & (adoptstm1[j]==1) & (adopts[j]==1)) {select_mat[i,j] <-14}
        if ((adoptstm1[i]==1) & (adopts[i]==1) & (adoptstm1[j]==1) & (adopts[j]==0)) {select_mat[i,j] <-15}
        if ((adoptstm1[i]==1) & (adopts[i]==1) & (adoptstm1[j]==1) & (adopts[j]==1)) {select_mat[i,j] <-16}
      }
    }
    select_mat_d  <- (select_mat * net_dropd)
      select_mat_s  <- (select_mat * net_stabl)
      select_mat <- (select_mat * net_added)

#apply(select_mat, 1, table (x, responseName=select))

      select1 <- rowSums((select_mat) == 1)
        select2 <- rowSums((select_mat) == 2)
        select3 <- rowSums((select_mat) == 3)
        select4 <- rowSums((select_mat) == 4)
        select5 <- rowSums((select_mat) == 5)
        select6 <- rowSums((select_mat) == 6)
        select7 <- rowSums((select_mat) == 7)
        select8 <- rowSums((select_mat) == 8)
        select9 <- rowSums((select_mat) == 9)
        select10 <- rowSums((select_mat) == 10)
        select11 <- rowSums((select_mat) == 11)
        select12 <- rowSums((select_mat) == 12)
        select13 <- rowSums((select_mat) == 13)
        select14 <- rowSums((select_mat) == 14)
        select15 <- rowSums((select_mat) == 15)
        select16 <- rowSums((select_mat) == 16)

        selectd1 <- rowSums((select_mat_d) == 1)
        selectd2 <- rowSums((select_mat_d) == 2)
        selectd3 <- rowSums((select_mat_d) == 3)
        selectd4 <- rowSums((select_mat_d) == 4)
        selectd5 <- rowSums((select_mat_d) == 5)
        selectd6 <- rowSums((select_mat_d) == 6)
        selectd7 <- rowSums((select_mat_d) == 7)
        selectd8 <- rowSums((select_mat_d) == 8)
        selectd9 <- rowSums((select_mat_d) == 9)
        selectd10 <- rowSums((select_mat_d) == 10)
        selectd11 <- rowSums((select_mat_d) == 11)
        selectd12 <- rowSums((select_mat_d) == 12)
        selectd13 <- rowSums((select_mat_d) == 13)
        selectd14 <- rowSums((select_mat_d) == 14)
        selectd15 <- rowSums((select_mat_d) == 15)
        selectd16 <- rowSums((select_mat_d) == 16)

        selects1 <- rowSums((select_mat_s) == 1)
        selects2 <- rowSums((select_mat_s) == 2)
        selects3 <- rowSums((select_mat_s) == 3)
        selects4 <- rowSums((select_mat_s) == 4)
        selects5 <- rowSums((select_mat_s) == 5)
        selects6 <- rowSums((select_mat_s) == 6)
        selects7 <- rowSums((select_mat_s) == 7)
        selects8 <- rowSums((select_mat_s) == 8)
        selects9 <- rowSums((select_mat_s) == 9)
        selects10 <- rowSums((select_mat_s) == 10)
        selects11 <- rowSums((select_mat_s) == 11)
        selects12 <- rowSums((select_mat_s) == 12)
        selects13 <- rowSums((select_mat_s) == 13)
        selects14 <- rowSums((select_mat_s) == 14)
        selects15 <- rowSums((select_mat_s) == 15)
        selects16 <- rowSums((select_mat_s) == 16)

        select <- cbind(time=Time, id=1:n, select1_=select1, select2_=select2, select3_=select3, select4_=select4, select5_=select5, select6_=select6, select7_=select7, select8_=select8, select9_=select9, select10_=select10,
                        select11_=select11, select12_=select12, select13_=select13, select14_=select14, select15_=select15, select16_=select16,
                        selectd1_=selectd1, selectd2_=selectd2, selectd3_=selectd3, selectd4_=selectd4, selectd5_=selectd5, selectd6_=selectd6, selectd7_=selectd7, selectd8_=selectd8, selectd9_=selectd9, selectd10_=selectd10,
                        selectd11_=selectd11, selectd12_=selectd12, selectd13_=selectd13, selectd14_=selectd14, selectd15_=selectd15, selectd16_=selectd16,
                        selects1_=selects1, selects2_=selects2, selects3_=selects3, selects4_=selects4, selects5_=selects5, selects6_=selects6, selects7_=selects7, selects8_=selects8, selects9_=selects9, selects10_=selects10,
                        selects11_=selects11, selects12_=selects12, selects13_=selects13, selects14_=selects14, selects15_=selects15, selects16_=selects16)
        return(select)



  }
}
*/



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

adopt <- adoptslist[[1]]$Adopt_mat
out2 <- select_egoalter(adjmat[,,1], adjmat[,,2], adopt[,1], adopt[,2])

library(microbenchmark)
microbenchmark(
  selectionFunctionEgoAlter(adjmat, adopt, 2),
  select_egoalter(adjmat[,,1], adjmat[,,2], adopt[,1], adopt[,2])
)

*/
