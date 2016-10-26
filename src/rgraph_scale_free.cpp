// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "netdiffuser_extra.h"
using namespace Rcpp;

// [[Rcpp::export]]
arma::sp_mat rgraph_ba_cpp(
    const arma::sp_mat & graph,
    const arma::colvec & dgr, int m = 1, int t = 10, bool self=true) {

  // Creating the empty graph
  // m0: Size of the graph (changes throughout time).
  // n: Final size of the graph.
  int m0 = graph.n_cols;
  int n = m0 + t;

  // Creating new graph and vector of degree
  arma::sp_mat graph_new(n,n);
  graph_new.submat(0,0,m0-1, m0-1) = graph;

  arma::colvec dgr_new(n, arma::fill::zeros);
  dgr_new.subvec(0, m0-1) = dgr;

  // Start the process, K is sum(dgr)
  int K = sum(dgr);

  int m_trunc;
  double randdraw, cump;

  // If self=true, then the prob are computed over m0+1, otherwise only over m0
  int extra = self? 1 : 0;
  for(int i=0;i<t;i++) {
    // Checling user interrup
    if (i % 1000 == 0)
      Rcpp::checkUserInterrupt();

    // The number of conections is trucated by the number of vertices in the graph
    // so -m_trunc- is actually -m-, but if currently there are less vertices than m,
    // then its truncated.
    m_trunc = (m > m0)? m0 : m;
    // if (m > m0) m0_trunc = m0;

    for (int j=0;j<m_trunc;j++) {
      // Incrementing the degree of the one that is been added
      // one by one until having degree m0.
      // Notice that m0 is updated each time, hence is equiv
      // to m0 = i + graph.n_cols
      dgr_new.at(m0) += 1.0;

      // std::cout << j << " Iter, "<< m << " m\n";
      // Random selection
      randdraw = unif_rand();

      // Calculating probabilities of been drawn. -cump- is the cumsum of the
      // probabilities
      cump = 0.0;
      for (int k=0; k<m0+extra; k++) {

        // In the case that in iter i the total degree is zero (no links)
        // then all individuals are equally likely to receive a new link.
        if (K != 0) cump += dgr_new.at(k)/(double)(K + extra);
        else cump += 1.0/(double)(m0 + extra);

        // DEBUGGING LINES, everything looks OK.
        // Rprintf("Prob (i:%02i, m:%02i): %-4.2g, randdraw: %-4.2g, m0:%02i\n", i+1, j+1, cump, randdraw, m0);
        // dgr_new.t().print();

        // Links to the set of previous vertices
        if (randdraw <= cump) {
          graph_new.at(m0, k) += 1.0; // , graph_new.at(k, m0) += 1.0
          dgr_new.at(k) += 1.0;

          // Sumation of degrees
          K += 2; // outgoing + incoming (2)
          break;
        }
      }
    }
    ++m0;
  }

  return graph_new;
}

// [[Rcpp::export]]
arma::sp_mat rgraph_ba_new_cpp(int m0 = 1, int m = 1, int t = 10, bool self=true) {
  int n = m0;

  arma::sp_mat graph(n, n);
  arma::colvec dgr(n);

  if (!self) {
    graph.diag().fill(0.0);
    dgr.fill(0.0);
  }
  else {
    graph.diag().fill(1.0);
    dgr.fill(2.0);
  }

  return rgraph_ba_cpp(graph, dgr, m, t, self);
}

// Implements the model described in:
// De Almeida, M. L., Mendes, G. A., Madras Viswanathan, G., & Da Silva, L. R. (2013).
// Scale-free homophilic network. European Physical Journal B, 86(2).
// \url{}http://doi.org/10.1140/epjb/e2012-30802-x}
// [[Rcpp::export]]
arma::sp_mat rgraph_sf_homo_cpp(
    const arma::colvec & eta,
    const arma::sp_mat & graph,
    const arma::colvec & dgr,
    int m = 1, int t = 10, bool self=true
  ) {

  // Creating the empty graph
  // m0: Size of the graph (changes throughout time).
  // n: Final size of the graph.
  int m0 = graph.n_cols;
  int n = m0 + t;

  // Creating new graph and vector of degree
  arma::sp_mat graph_new(n,n);
  graph_new.submat(0,0,m0-1, m0-1) = graph;

  arma::colvec dgr_new(n, arma::fill::zeros);
  dgr_new.subvec(0, m0-1) = dgr;

  // Computing similarity
  arma::colvec Ai(n - (self? 0 : 1) ); // similitude
  arma::colvec etanorm = (eta - min(eta))/(max(eta) - min(eta) + 1e-15);
  arma::colvec K1Ai(Ai.size());

  int m_trunc;
  double randdraw, cump, sum_1AK;

  // If self=true, then the prob are computed over m0+1, otherwise only over m0
  int extra = self? 1 : 0;
  for(int i=0;i<t;i++) {
    // Checling user interrup
    if (i % 1000 == 0)
      Rcpp::checkUserInterrupt();

    // The number of conections is trucated by the number of vertices in the graph
    // so -m0_trunc- is actually -m-, but if currently there are less vertices than m,
    // then its truncated.
    m_trunc = (m > m0)? m0 : m;
    // if (m > m0) m0_trunc = m0;

    // Calculating similitude as A(i) = |eta - eta_i|
    for (int j=0;j<(m0 + extra);j++)
      Ai.at(j) = fabs(etanorm.at(j) - etanorm.at(i));

    // If there are no links in the graph then all individuals are likeli
    // prop to eta.
    double dgrsum=0.0;
    if (i==0) dgrsum = sum(dgr_new);
    if (i==0 && (dgrsum == 0.0))  {
      for (int j=0; j<m0+extra;j++)
        K1Ai.at(j) = (1-Ai.at(j));
    } else {
      for (int j=0; j<m0+extra;j++)
        K1Ai.at(j) = (1-Ai.at(j)) * dgr_new.at(j);
    }

    // Denominator sum_j A(ij)*K(j)
    if (self) K1Ai.at(m0) = 1.0;

    for (int j=0;j<m_trunc;j++) {

      // Incrementing the degree of the one that is been added
      // one by one until having degree m0_trunk.
      dgr_new.at(m0) += 1.0;

      // std::cout << j << " Iter, "<< m << " m\n";
      // Random selection
      randdraw = unif_rand();

      // Calculating probabilities of been drawn. -cump- is the cumsum of the
      // probabilities
      cump = 0.0;
      sum_1AK = sum(K1Ai.subvec(0,m0+extra-1));
      for (int k=0; k<m0+extra; k++) {

        // Accumulating probability
        cump += K1Ai.at(k)/(sum_1AK);

        if (cump > 1.0+1e-15) printf("over1 (%04d/%04d): %9.4g\n", k+1, m0+1,cump );

        // Links to the set of previous vertices
        if (randdraw <= cump) {
          graph_new.at(m0, k) += 1.0; // , graph_new.at(k, m0) += 1.0
          dgr_new.at(k) += 1.0;

          // Probabilities
          K1Ai.at(k) = (1 - Ai.at(k))*dgr_new.at(k);
          if (self) K1Ai.at(m0) = dgr_new.at(m0);

          break;
        }
      }
    }
    ++m0;
  }

  return graph_new;
}

// [[Rcpp::export]]
arma::sp_mat rgraph_sf_homo_new_cpp(const arma::colvec & eta, int m0 = 1, int m = 1, int t = 10, bool self=true) {
  int n = m0;

  arma::sp_mat graph(n, n);
  arma::colvec dgr(n);
  // Rcpp::stop("-eta- should be of length m0+t.");
  // eta should be of length m0 + t
  if (((int) eta.n_elem) != (m0 + t))
    Rcpp::stop("-eta- should be of length m0+t.");

  if (!self) {
    graph.diag().fill(0.0);
    dgr.fill(0.0);
  }
  else {
    graph.diag().fill(1.0);
    dgr.fill(2.0);
  }

  return rgraph_sf_homo_cpp(eta, graph, dgr, m, t, self);
}

/** *R
 library(Matrix)
 set.seed(123)
 n <- 1e3
 eta <- cbind(as.numeric(sample(c(0,1,1000,1001), n, TRUE)))
 ans <- rgraph_sf_homo_new_cpp(eta, 1,5,n-1)
# ans

 ig <- igraph::graph_from_adjacency_matrix(ans)
 plot(ig, vertex.size=5, vertex.color=factor(eta), layout=igraph::layout_with_fr, vertex.label=NA)
 */
