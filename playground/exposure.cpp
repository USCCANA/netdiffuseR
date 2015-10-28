// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

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
arma::mat exposure_cpp(const arma::mat & adjmat, const arma::mat adopt,
                  int which = 0) {
  int n = adjmat.n_cols;
  int T = adopt.n_cols;

  List SE_OBJ(3);
  if (which == 1)
    SE_OBJ = struct_equiv_cpp(adjmat, 1.0, true, true);

  arma::colvec degree(n);
  if (which > 1)
    degree = degree_cpp(adjmat, which - 2);

  // Output
  arma::mat exposure(n,T);

  // Computing weights
  arma::colvec NUMERATOR(n);
  for(int t=0;t<T;t++) {
    switch ( which ) {
    case 0: {// Unweighted
      NUMERATOR = adjmat * adopt.col(t);
      break;
    }
    case 1: {// SE
      NUMERATOR = as<arma::mat>(SE_OBJ["d"]);
      NUMERATOR = (adjmat % NUMERATOR) * adopt.col(t);
      break;
    }
    default:
      if (which > 1 && which <= 4) { // Degree
        NUMERATOR = adjmat * (adopt.col(t) % degree);
        break;
      }
      stop("Invalid weight code.");
    }

    exposure.col(t) = NUMERATOR / sum(adjmat,1);
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
  }

  res <- list(Exposure = Exposure, ExposureSE = ExposureSE, ExposureC = ExposureC)
    return(res)
}*/

/***R
set.seed(123)
adjmat <- rand_graph_cpp()
adopt <- adopt_mat_cpp(sample(1:3, 10, TRUE))
exposure_cpp(adjmat, adopt$adoptmat)

adjmatarray <- array(adjmat, dim=c(10,10,3))
library(sna)
library(network)
diffusiontest::ExposureCalc(adjmatarray, adopt$adoptmat)
*/
