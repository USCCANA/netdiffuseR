#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "adjmat.hpp"
#include "struct_equiv.hpp"

using namespace Rcpp;

/* which:
 * 0: Unweighted
 * 1: Structure equiv weighted
 * 2: Indegree
 * 3: Outdegree
 * 4: Degree
 */

// [[Rcpp::export]]
arma::cube inputoutputarray(NumericVector x) {
  IntegerVector dims = x.attr("dim");
  arma::cube y(x.begin(), dims[0], dims[1], dims[2], false);
  return y;
}

// [[Rcpp::export]]
arma::mat exposure_cpp(NumericVector dynmat, const arma::mat adopt,
                  int wtype = 0) {

  IntegerVector dims=dynmat.attr("dim");
  const arma::cube dynadjmat(dynmat.begin(), dims[0], dims[1], dims[2], false);

  // Variables initialization
  const int n = dynadjmat.n_cols;
  const int T = dynadjmat.n_slices;
  arma::mat exposure(n,T, arma::fill::zeros);

  List se(3);
  arma::colvec degree(n);
  arma::mat sedist(n,n);

  arma::colvec NUMERATOR(n);
  arma::colvec DENOMINATOR(n);
  arma::colvec RATIO(n);

  arma::mat subgraph(n,n);

  for(int i=0;i<T;i++) {
    // Computing weights
    if (wtype==0) { // Unweighted
      NUMERATOR = dynadjmat.slice(i) * adopt.col(i);
    }
    /*else if (wtype == 1) {// SE
      se = struct_equiv_cpp(subgraph, 1.0, true, true);
      sedist = as< arma::mat >(se["d"]);
      NUMERATOR = (subgraph % sedist) * adopt.col(i);
    }
    else if (wtype > 1 && wtype <= 4) { // Degree
      degree = degree_cpp(subgraph, wtype - 2);
      NUMERATOR = subgraph * (adopt.col(i) % degree);
    }
    else {
      stop("Invalid weight code.");
    }*/
    // if (i==1) return exposure;
    DENOMINATOR = sum(dynadjmat.slice(i),1);
    // arma::colvec suma=NUMERATOR / sum(dynadjmat.slice(t),1);

    Rprintf("Num ncols:%d nrows:%d\n", NUMERATOR.n_cols, NUMERATOR.n_rows);
    Rprintf("Den ncols:%d nrows:%d\n", DENOMINATOR.n_cols, DENOMINATOR.n_rows);
    for(int j=0;j<n;j++) Rprintf("%02.4f,",NUMERATOR[j]);
    Rprintf("\n");
    for(int j=0;j<n;j++) Rprintf("%02.4f,",DENOMINATOR[j]);
    Rprintf("\n");

    // return exposure;
    // exposure.col(i) = (NUMERATOR/DENOMINATOR);

  }

  return exposure;
}
/*
 ExposureCalc <- function(all_nets, Adopt_mat){
 n <- dim(all_nets)[1]
 maxTime <- ncol(Adopt_mat)
 Exposure   <- matrix(data = 0, nrow = n, ncol = (maxTime), byrow = TRUE)
 ExposureSE <- matrix(data = 0, nrow = n, ncol = (maxTime), byrow = TRUE)
 ExposureC  <- matrix(data = 0, nrow = n, ncol = (maxTime), byrow = TRUE)
 dynmatfor(int j=0;j<n;j++) Rprintf("%02.4f,",NUMERATOR[j]);
  for(Time in 1:maxTime){
#Anything that is multiplied by adoption is specific by convention
    if(Time < maxTime){
      adjmat_mat <- all_nets[,,Time]
      semat_mat <- (1 /(sedist(adjmat_mat, method="euclidean")))
      semat_mat[is.infinite(semat_mat[,])]<- 0
      diag(semat_mat) <- 0
      in_deg <- (degree(as.network(adjmat_mat), cmode="indegree"))
      ExposureC[,Time+1]    <-((adjmat_mat %*% (Adopt_mat[,Time] * in_deg)) / (rowSums(adjmat_mat)+.0001))
      Exposure[,Time+1]   <-((adjmat_mat %*% Adopt_mat[,Time]) / (rowSums(adjmat_mat)+.0001))
      ExposureSE[,Time+1] <-((semat_mat %*% Adopt_mat[,Time]) / (rowSums(semat_mat)+.0001))
    }
  }arma::colvec NUMERATOR(n);

  res <- list(Exposure = Exposure, ExposureSE = ExposureSE, ExposureC = ExposureC)
    return(res)
}*/

/** *R
library(sna)
library(network)
library(diffusiontest)
set.seed(123)
grapj <- array(unlist(lapply(1:3, function(...) rand_graph_cpp())), c(10,10,3))
adopt <- adopt_mat_cpp(sample(1:3, 10, TRUE))
exposure_cpp(grapj, adopt$adoptmat)

diffusiontest::ExposureCalc(graph, adopt$adoptmat)
*/
