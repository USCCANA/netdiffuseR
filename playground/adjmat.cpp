#include <Rcpp.h>
using namespace Rcpp;

/*
// [[Rcpp::export]]
NumericMatrix adjmat(const NumericMatrix edgelist) {

  // Getting the dimensions and creating objects
  int m = edgelist.nrow();
  int k = edgelist.ncol();

  // Identifying the unique nodes
  NumericVector nodes;
  for(int i=0;i<m;m++)
    nodes.push_back()

  NumericMatrix mat(n,n);



  NumericVector nodes =
}*/

// [[Rcpp::export]]
NumericVector uniquecpp(NumericVector x) {return unique(x);}

/***R
x <- c(1,1,2,3,4,4,5,1)
uniquecpp(x)
*/
