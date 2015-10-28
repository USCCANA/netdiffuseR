#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "adjmat.hpp"
#include "struct_equiv.hpp"

using namespace Rcpp;

// [[Rcpp::export]]
List exposure_cpp(const arma::mat & adjmat, const arma::colvec adopt) {
  int n = adjmat.n_cols;
  int T = adopt.n_cols;

  // Exposure weighted by SE
  List SE_OBJ = struct_equiv_cpp(adjmat, 1.0, true, true);
  arma::mat SE = SE_OBJ["d"];
  arma::mat EX = ((adjmat%SE)*adopt)/sum(adjmat,1);

  // Exposure weighted by Indegree
  arma::mat degree = degree_cpp(adjmat,);
  return List::create(_["exp_se"]=EX,_["exp"]);
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
adjmat <- sna::rgraph(10)
adopt <- adopt_mat_cpp(sample(1, 10, TRUE))
exposure_cpp(adjmat, adopt$adoptmat[,1])

adjmatarray <- array(adjmat, dim=c(10,10,1))
library(sna)
diffusiontest::ExposureCalc(adjmatarray, adopt$adoptmat)
*/
