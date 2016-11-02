/*******************************************************************************
* rgraph.h header function for rgraph.cpp: Functions for gerating random graphs
* and rewiring graphs.
*******************************************************************************/

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#ifndef NETDIFFUSER_RGRAPH_
#define NETDIFFUSER_RGRAPH_
using namespace Rcpp;

arma::sp_mat rgraph_er_cpp(
    int n=10,
    double p = 0.3,
    bool undirected=true,
    bool weighted=false,
    bool self=false
  );

arma::sp_mat ring_lattice(
    int n,
    int k,
    bool undirected=true
  );

arma::sp_mat rewire_swap(
    const arma::sp_mat & graph,
    int nsteps=100,
    bool self=false,
    bool multiple=false,
    bool undirected=false,
    double pr_rewire=0.5 //,bool althexagons=false
);

arma::sp_mat rewire_endpoints(
    const arma::sp_mat & graph,
    double p,
    bool both_ends=false,
    bool self=false,
    bool multiple=false,
    bool undirected=false
  );

arma::sp_mat permute_graph_cpp(
    const arma::sp_mat & x,
    bool self = false,
    bool multiple = false
  );

arma::sp_mat rgraph_ba_cpp(
    arma::sp_mat graph,
    arma::colvec dgr,
    int m = 1,
    int t = 10,
    bool self=true
  );

arma::sp_mat rgraph_ba_new_cpp(
    int m0 = 1,
    int m = 1,
    int t = 10,
    bool self=true
  );

arma::sp_mat rgraph_sf_homo(
    const arma::colvec & eta,
    const arma::sp_mat & graph,
    const arma::colvec & dgr,
    int m = 1,
    int t = 10,
    bool self=true
  );

arma::sp_mat rgraph_sf_homo_new(
    const arma::colvec & eta,
    int m0 = 1,
    int m = 1,
    int t = 10,
    bool self=true
  );

#endif
