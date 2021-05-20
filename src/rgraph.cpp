// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "netdiffuser_extra.h"
using namespace Rcpp;

// -----------------------------------------------------------------------------
//
// Bernoulli graphs
//
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
arma::sp_mat rgraph_er_cpp(
    int n=10, double p = 0.3, bool undirected=true,
    bool weighted=false, bool self=false) {

  std::vector< unsigned int > source;
  std::vector< unsigned int > target;
  std::vector< double > vals;

  for(int i=0;i<n;++i) {
    // Checling user interrup
    if (i % 200 == 0)
      Rcpp::checkUserInterrupt();

    /* Setting the length of the subloop acordingly to type of graph */
    int m = n;
    if (undirected) m=i;
    for(int j=0;j<m;++j) {

      /* Assessing if include self */
      if (!self && (i==j)) continue;

      /* Setting the value of the tie */
      double w = unif_rand(); //R::runif(0, 1);
      if (w > (1-p)) {
        if (!weighted) w=1.0;

        source.push_back(i);
        target.push_back(j);
        vals.push_back(w);

        // Adding the mirror
        if (undirected) {
          source.push_back(j);
          target.push_back(i);
          vals.push_back(w);
        }
      }
    }
  }

  // Storing the data as sparse matrix
  arma::sp_mat graph(
      // Creating a umat
      arma::join_cols(
        arma::conv_to< arma::urowvec >::from(source),
        arma::conv_to< arma::urowvec >::from(target)
    ),
    // Values
    arma::conv_to< arma::colvec >::from(vals),
    // Size, sorting and checking for zeros
    n, n, true, false);

  return graph;
}

// -----------------------------------------------------------------------------
//
// Small world graphs
//
// -----------------------------------------------------------------------------


//' Ring lattice graph
//'
//' Creates a ring lattice with \eqn{n} vertices, each one of degree (at most) \eqn{k}
//' as an undirected graph. This is the basis of \code{\link{rgraph_ws}}.
//' @param n Integer scalar. Size of the graph.
//' @param k Integer scalar. Out-degree of each vertex.
//' @param undirected Logical scalar. Whether the graph is undirected or not.
//' @details when \code{undirected=TRUE}, the degree of each node always
//' even. So if \code{k=3}, then the degree will be \code{2}.
//' @return A sparse matrix of class \code{\link[Matrix:dgCMatrix-class]{dgCMatrix}} of size
//' \eqn{n\times n}{n * n}.
//' @references Watts, D. J., & Strogatz, S. H. (1998). Collective dynamics of
//' “small-world” networks. Nature, 393(6684), 440–2. \doi{10.1038/30918}
//' @export
//' @family simulation functions
// [[Rcpp::export]]
arma::sp_mat ring_lattice(int n, int k, bool undirected=false) {

  if ((n-1) < k)
    stop("k can be at most n - 1");

  arma::sp_mat graph(n,n);

  // Adjusting k
  if (undirected)
    if (k>1) k = (int) floor((double) k/2.0);

  // Connecting to k/2 next & previous neighbour
  for (int i=0;i<n;++i) {
    for (int j=1;j<=k;++j) {
      // Next neighbor
      int l = i+j;
      if (l >= n) l = l - n;

      graph.at(i,l) += 1.0;
      if (undirected) graph.at(l,i) += 1.0;
    }
  }

  return graph;
}

/** *R
library(Matrix)
library(sna)
x <- ring_lattice(6,2)
x
gplot(as.matrix(x), displaylabels = TRUE, mode="circle", jitter = FALSE)

# x <- ring_cpp(10,6)
# x
# gplot(as.matrix(x), displaylabels = TRUE, mode="circle", jitter = FALSE)
#
# x <- ring_cpp(7,7)
# x
# gplot(as.matrix(x), displaylabels = TRUE, mode="circle", jitter = FALSE)
*/

// [[Rcpp::export]]
arma::sp_mat rewire_endpoints(
    const arma::sp_mat & graph, double p,
    bool both_ends=false,
    bool self=false, bool multiple=false,
    bool undirected=false) {

  // Clonning graph
  int n = graph.n_cols;
  arma::sp_mat newgraph(graph);

  // Getting the indexes

// for (unsigned i= 0;i<indexes.n_rows; ++i) {
  unsigned int i = 0;
  for (arma::sp_mat::const_iterator it = graph.begin(); it != graph.end(); ++it) {

    // Checking user interrupt
    if (++i % 1000 == 0)
      Rcpp::checkUserInterrupt();

    // Checking whether to change it or not
    if (unif_rand() > p) continue;

    // Indexes
    int j = it.row();
    int k = it.col();

    // In the case of undirected graphs, we only modify the lower triangle
    // The upper triangle part will be rewritten during the rand.
    if (undirected && (j < k)) continue;

    // New end(s)
    int newk, newj;

    // If rewiring is in both sides, then we start from the left
    if (both_ends) newj = (int) floor( unif_rand()*n );
    else           newj = j;

    if (undirected) newk = (int) floor( unif_rand()*(newj + 1));
    else            newk = (int) floor( unif_rand()*n);

    // Self edges are not allowed, again, must check this on the new graph
    if (!self && (newj == newk)) continue;

    // Multiple edges are not allowed. Must check this on the new graph
    if (!multiple && (newgraph.at(newj, newk) != 0)) continue;

    // Adding up
    double w = graph.at(j,k);

    newgraph.at(j,k) = 0;
    if (undirected) newgraph.at(k,j) = 0;

    newgraph.at(newj,newk) += w;
    if (undirected) newgraph.at(newk,newj) += w;

  }
  return newgraph;
}

// [[Rcpp::export]]
arma::sp_mat rewire_swap(
    const arma::sp_mat & graph, int nsteps=100,
    bool self=false, bool multiple=false,
    bool undirected=false, double pr_rewire=0.5 //,bool althexagons=false
) {

  // Clonning graph
  arma::sp_mat newgraph(graph);

  // Getting the indexes
  arma::umat indexes(graph.n_nonzero, 2);
  unsigned int m = 0;
  for (arma::sp_mat::const_iterator it = newgraph.begin(); it != newgraph.end(); ++it) {

    // Checking cases
    if      (!self      && (it.row() == it.col())) continue;
    else if (undirected && (it.row() <  it.col())) continue;

    // Filling the matrix
    indexes.at(m,0) = it.row();
    indexes.at(m++,1) = it.col();

  }

  // Shedding rows
  if (m < indexes.n_rows)
    indexes.shed_rows(m, indexes.n_rows - 1u);

  // double dens = graph.n_nonzero/(graph.n_cols*graph.n_cols);
  int s = 0;
  while (s++ < nsteps) {

    // Does it need to be rewired?
    if (unif_rand() > pr_rewire) continue;

    // Checking user interrupt
    if (s % 1000 == 0)
      Rcpp::checkUserInterrupt();

    // Indexes
    int ij = floor(unif_rand()*indexes.n_rows);
    int i = indexes.at(ij, 0);
    int j = indexes.at(ij, 1);

    // If self, then remove from the indexes and go for the next
    if (!self && i==j) continue;

    // New end(s) (and values)
    int newij, newi, newj;
    // int i1,i2,i3;
    // double v2;
    double v0,v1;

    // Choosing the other rand
    newij = floor( unif_rand()*indexes.n_rows );

    // The newij shouldn't be the same as ij!
    if (ij == newij) continue;

    newi = indexes.at(newij, 0);
    newj = indexes.at(newij, 1);

    bool ismultiple = !multiple &&
      (newgraph.at(i, newj) != 0) | (newgraph.at(newi, j) != 0);

    /*// Alternating Hexagons
    // Ramachandra Rao, et al, The Indian Journal of Statistics
    if (ismultiple && althexagons && (unif_rand() < 0.5)) {
    // Case 1: will switch i1i2 and i2i3 (so j == newi)
    if ((j==newi) && (newgraph.at(newj, i))) {

    i1 = i, i2 = j, i3 = newj;

    // All alt are zero
    if ((newgraph.at(i1,i3) != 0) ||
    (newgraph.at(i2,i1) != 0) ||
    (newgraph.at(i3,i2) != 0)) continue;

    v0 = newgraph.at(i1, i2);
    v1 = newgraph.at(i2, i3);
    v2 = newgraph.at(i3, i1);

    newgraph.at(i1, i2) = 0;
    newgraph.at(i2, i3) = 0;
    newgraph.at(i3, i1) = 0;

    newgraph.at(i1, i3) = v0;
    newgraph.at(i2, i1) = v1;
    newgraph.at(i3, i2) = v2;

    } else if ((i==newj) && (newgraph.at(j, newi))) {

    i1 = i, i2 = j, i3 = newi;

    // All alt are zero
    if ((newgraph.at(i1,i3) != 0) ||
    (newgraph.at(i2,i1) != 0) ||
    (newgraph.at(i3,i2) != 0)) continue;

    v0 = newgraph.at(i1, i2);
    v1 = newgraph.at(i2, i3);
    v2 = newgraph.at(i3, i1);

    newgraph.at(i1, i2) = 0;
    newgraph.at(i2, i3) = 0;
    newgraph.at(i3, i1) = 0;

    newgraph.at(i1, i3) = v0;
    newgraph.at(i2, i1) = v1;
    newgraph.at(i3, i2) = v2;

    }

    continue;
    }*/

    // Checking multiple
    if (ismultiple)
      continue;

    // If any of the edges coincides and no self edge is allowed,
    // then we must check for self pointing
    if (!self && ((i == newj) | (newi == j))) continue;

    // Cleaning 2: Changing the graph patterns
    v0 = newgraph.at(i, j);
    v1 = newgraph.at(newi, newj);

    newgraph.at(i, j) = 0.0;
    if (undirected) newgraph.at(j, i) = 0.0;

    newgraph.at(newi, newj) = 0.0;
    if (undirected) newgraph.at(newj, newi) = 0.0;

    // Setting new values
    newgraph.at(i, newj) += v0;
    if (undirected) newgraph.at(newj, i) += v0;

    newgraph.at(newi, j) += v1;
    if (undirected) newgraph.at(j, newi) += v1;

    // And changing indexes
    indexes.at(ij,1) = newj;
    indexes.at(newij,1) = j;
  }

  return newgraph;
}



/** *R
x <- netdiffuseR::ring_lattice(10, 2)
set.seed(19237)

t(cbind(netdiffuseR::dgr(x, "indegree"), netdiffuseR::dgr(x, "outdegree")))
x2 <- rewire_swap(x, 1)
t(cbind(netdiffuseR::dgr(x2, "indegree"), netdiffuseR::dgr(x2, "outdegree")))

# Benchmark
library(igraph)
library(microbenchmark)
library(netdiffuseR)

set.seed(1133)
x <- barabasi.game(1e4)
y <- as_adj(x)


microbenchmark(
  ig      = rewire(x, keeping_degseq(niter = 100)),
  nd      = netdiffuseR:::rewire_swap(y, 100),
  unit="relative"
)

*/


// [[Rcpp::export]]
arma::sp_mat rewire_ws(arma::sp_mat G, int K, double p=0.0,
                       bool self=false, bool multiple=false) {

  arma::sp_mat out(G);
  int n = G.n_rows;

  // First half
  for(int k=1;k<=K;++k) {
    for(int i=0;i<n;++i) {
      // Clock wise choose
      int j;
      if (k <= K/2) j = ((i + k) < n)? i + k: k - (n - i);
      else {
        int tmpk = k - K/2;
        j = ((i - tmpk) < 0)? n - tmpk + i: i - tmpk;
      }

      // If not rewire then leave as is
      if (unif_rand() > p) continue;

      // Picking the new random obs excluding i
      int newj;
      if (!self) newj = unif_rand_w_exclusion(n, i);
      else       newj = floor(unif_rand()*n);

      // If multiple
      std::vector< bool > checked(n);
      int nchecked = 0;
      while (!multiple && out.at(i,newj) != 0) {
        // Picking a new one
        if (!self) newj = unif_rand_w_exclusion(n, i);
        else       newj = floor(unif_rand()*n);

        // Has it already been drawn?
        if (!checked.at(newj)) {
          checked.at(newj) = true;
          ++nchecked;

        } else {
          if      ( self && (nchecked >= n))       break;
          else if (!self && (nchecked >= (n - 1))) break;
        }
      }

      // If multiple, then continue
      if (multiple && out.at(i, newj) != 0) continue;


      // Changing values
      double v = G.at(i, j);
      out.at(i,j) = 0;
      out.at(i, newj) = v;


    }
  }
  // // Second half
  // for(int k=1;k<=K/2;++k) {
  //   for(int i=0;i<G.n_cols;++i) {
  //     // Clockwise choose
  //     int j = ((i - k) < 0)? G.n_cols - k + i: i - k;
  //     // // Rprintf("(%d, %d)\n", i, j);
  //     // out.at(i,j) = 1;
  //   }
  // }

  return out;
}


/** *R
rgraph_ws <- function(n,k,p, both_ends=FALSE, self=FALSE, multiple=FALSE) {
rewire_endpoints(ring_lattice(n, k), p, both_ends,
                 self, multiple, true)
}

x <- ring_lattice(14, 2)
gplot(as.matrix(x), mode="circle", jitter=FALSE, usecurve = TRUE, gmode = "graph")
gplot(as.matrix(netdiffuseR:::rewire_endpoints(x, .1)), mode="circle", jitter=FALSE, usecurve = TRUE, gmode = "graph")
gplot(as.matrix(netdiffuseR:::rewire_endpoints(x, 1)), mode="circle", jitter=FALSE, usecurve = TRUE, gmode = "graph")
*/

// [[Rcpp::export]]
arma::sp_mat permute_graph_cpp(const arma::sp_mat & x,
                               bool self = false,
                               bool multiple = false) {

  // Initialization
  typedef arma::sp_mat::const_iterator spiter;
  spiter beg = x.begin();
  spiter end = x.end();
  int n = x.n_rows;
  arma::sp_mat ans(x.n_rows, x.n_cols);
  bool keeplooking;
  int niter = 0;
  for(spiter iter=beg;iter!=end;++iter) {

    // Checking user interrupt
    if (++niter % 1000 == 0)
      Rcpp::checkUserInterrupt();

    keeplooking = true;
    while (keeplooking) {

      // Proposal
      int newi = unif_rand_w_exclusion(n,-1);
      int newj = unif_rand_w_exclusion(n, (self? -1 : newi));

      if (!multiple && ans.at(newi, newj) != 0) continue;

      ans.at(newi, newj) += *iter;
      keeplooking = false;

    }
  }

  return ans;
}
//
// /***R
// set.seed(1123)
// N <- 1e2
// g <- netdiffuseR:::rgraph_ba(t=9)
//
// sapply(1:N, function(x) sum(permute_graph(g)))
//
// */

// -----------------------------------------------------------------------------
//
// Scale-free graphs
//
// -----------------------------------------------------------------------------

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
  arma::colvec dgr_new(n, arma::fill::zeros);
  dgr_new.subvec(0, m0-1) = dgr;

  std::vector< std::vector<unsigned int> > source(n);

  // Setting the initial values
  int nlocations = graph.n_nonzero;
  arma::sp_mat::const_iterator iter;
  for (iter = graph.begin(); iter != graph.end(); ++iter) {
    source[iter.row()].push_back(iter.col());
  }


  // Start the process, K is sum(dgr)
  int K = sum(dgr);

  int m_trunc;
  double randdraw, cump;

  // If self=true, then the prob are computed over m0+1, otherwise only over m0
  int extra = self? 1 : 0;
  for(int i=0;i<t;++i) {
    // Checling user interrup
    if (i % 1000 == 0)
      Rcpp::checkUserInterrupt();

    // The number of conections is trucated by the number of vertices in the graph
    // so -m_trunc- is actually -m-, but if currently there are less vertices than m,
    // then its truncated.
    m_trunc = (m > m0)? m0 : m;
    // if (m > m0) m0_trunc = m0;

    for (int j=0;j<m_trunc;++j) {
      // Incrementing the degree of the one that is been added
      // one by one until having degree m0.
      // Notice that m0 is updated each time, hence is equiv
      // to m0 = i + graph.n_cols
      dgr_new.at(m0) += 1.0;

      // std::cout << j << " Iter, "<< m << " m\n";
      // Random selection
      randdraw = R::runif(0,1); //unif_rand();

      // Calculating probabilities of been drawn. -cump- is the cumsum of the
      // probabilities
      cump = 0.0;
      for (int k=0; k<m0+extra; ++k) {

        // In the case that in iter i the total degree is zero (no links)
        // then all individuals are equally likely to receive a new link.
        if (K != 0) cump += dgr_new.at(k)/(double)(K + extra);
        else cump += 1.0/(double)(m0 + extra);

        // DEBUGGING LINES, everything looks OK.
        // Rprintf("Prob (i:%02i, m:%02i): %-4.2g, randdraw: %-4.2g, m0:%02i\n", i+1, j+1, cump, randdraw, m0);
        // dgr_new.t().print();

        // Links to the set of previous vertices
        if (randdraw <= cump) {

          source[m0].push_back(k);
          ++nlocations;

          dgr_new.at(k) += 1.0;

          // Sumation of degrees
          K += 2; // outgoing + incoming (2)
          break;
        }
      }
    }
    ++m0;
  }

  // Now need to coerce into a -sp_mat-
  arma::umat locations(2, nlocations);
  arma::colvec values(nlocations, arma::fill::ones);

  int curloc = 0;
  for (int i = 0; i<n; ++i) {
    for (unsigned int j = 0; j < source[i].size(); ++j) {
      locations.at(0, curloc) = i;
      locations.at(1, curloc++) = source[i][j];
    }

  }

  // Creating the sparse matrix
  arma::sp_mat graph_new(true, locations, values, n, n, true, false);

  return graph_new;
}


// The bag algorithm should be straight forward.
// For each time that node i receives a link, add it to the bag, which we
// know what is the size in the beginning.
// Later, to draw a new vertex, it suffices to bag[runif()*length(bag)],
// so we get him from the bag.

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
arma::sp_mat rgraph_sf_homo(
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
  for(int i=0;i<t;++i) {
    // Checling user interrup
    if (i % 1000 == 0)
      Rcpp::checkUserInterrupt();

    // The number of conections is trucated by the number of vertices in the graph
    // so -m0_trunc- is actually -m-, but if currently there are less vertices than m,
    // then its truncated.
    m_trunc = (m > m0)? m0 : m;
    // if (m > m0) m0_trunc = m0;

    // Calculating similitude as A(i) = |eta - eta_i|
    for (int j=0;j<(m0 + extra);++j)
      Ai.at(j) = fabs(etanorm.at(j) - etanorm.at(i));

    // If there are no links in the graph then all individuals are likeli
    // prop to eta.
    double dgrsum=0.0;
    if (i==0) dgrsum = sum(dgr_new);
    if (i==0 && (dgrsum == 0.0))  {
      for (int j=0; j<m0+extra;++j)
        K1Ai.at(j) = (1-Ai.at(j));
    } else {
      for (int j=0; j<m0+extra;++j)
        K1Ai.at(j) = (1-Ai.at(j)) * dgr_new.at(j);
    }

    // Denominator sum_j A(ij)*K(j)
    if (self) K1Ai.at(m0) = 1.0;

    for (int j=0;j<m_trunc;++j) {

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
      for (int k=0; k<m0+extra; ++k) {

        // Accumulating probability
        cump += K1Ai.at(k)/(sum_1AK);

        // if (cump > 1.0+1e-15) printf("over1 (%04d/%04d): %9.4g\n", k+1, m0+1,cump );

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
arma::sp_mat rgraph_sf_homo_new(const arma::colvec & eta, int m0 = 1, int m = 1, int t = 10, bool self=true) {
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

  return rgraph_sf_homo(eta, graph, dgr, m, t, self);
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
