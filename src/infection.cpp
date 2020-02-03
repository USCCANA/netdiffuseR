// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector infection_cpp(
    List graph,
    const arma::colvec & times,
    bool normalize = true,
    int K = 1,
    double r = 0.5,
    bool expdiscount = false,
    int n=0,
    bool valued=false,
    bool outgoing=true) {

//   // Coersing a NumericVector into a cube for ease of use
//   IntegerVector dims=graph.attr("dim");
//   const arma::cube graph_cube(graph.begin(), dims[0], dims[1], dims[2], false);

  // Variables initialization
//   const int n = graph_cube.n_rows;
//   const int T = graph_cube.n_slices;
  NumericVector infect(n);

  // Variables to use within loop
  int ti, tj;

  // Creating discount variable, firts, must be truncated.
  if (K >= graph.size()) {
    warning("Too many periods selected, will be truncated to T-1.");
    K = graph.size() - 1;
  }

  // // Checking classes
  // for (int i = 0; i< graph.length();++i)
  //   if ()

  // The discount can be either exponential (1 + r)^(k-1), or
  // lineal in the form of k.
  // double * discount = new double[K]
  std::vector< double > discount(K);
  if (!expdiscount) {
    for(int k=1;k<=K;++k)
      discount[k-1] = k;
  }
  else {
    for(int k=1;k<=K;++k)
      discount[k-1] = pow((1.0 + r), k-1.0);
  }

  for(int i=0;i<n;++i) {
    // Checling user interrup
    if (i % 1000 == 0)
      Rcpp::checkUserInterrupt();

    // If NA (aka nan in Armadillo), then NA.
    if (!arma::is_finite(times(i))) {
      infect.at(i) = NA_REAL;
      continue;
    }

    // Capturing variables
    ti = times(i);
    if (ti == graph.size()) continue;

    double numerator = 0.0;
    double denominator = 0.0;

    // For the adjusted verion, see the mathematical supplement on Valente et al. (2015)
    double nadopt_t = 0;

    for(int k=1;k<=K;++k) {

      // Current time period from 1 to T
      int t = ti + k - 1;

      // If the required time period does not exists, then continue, recall that
      // vectors can be reach up to T - 1
      if (t >= graph.size()) continue;

      // Storing the graph
      arma::sp_mat graph_cube = graph.at(t);
      if (!valued)   graph_cube = arma::spones(graph_cube);
      if (!outgoing) graph_cube = graph_cube.t();

      for(int j=0;j<n;++j) {
        if (i==j) continue;

        tj = times(j);
        // Adding up for t+k iff a link between j and i exists
        // if (graph_cube(j,i,t) != 0) {
        if (graph_cube(j,i) != 0) {
          if (tj ==(ti+k)) {
            numerator += 1.0 / discount[k-1];
          }
          // Adding up for t+1 <= t <= T
          if (tj >=(ti+k)) denominator += 1.0 / discount[k-1];
        }

        // Has adopted so far? (discounted version)
        if (tj == (ti+k)) nadopt_t += 1.0 / discount[k-1];
      }
    }

    // Putting all together
    infect(i) = (numerator / (denominator + 1e-15));
    if (normalize) infect(i) = infect(i) / (nadopt_t + 1e-15);
  }

  return infect;
}

// [[Rcpp::export]]
NumericVector susceptibility_cpp(
    List graph,
    const arma::colvec & times,
    bool normalize = true,
    int K = 1,
    double r = 0.5,
    bool expdiscount = false,
    int n=0, bool valued=false,
    bool outgoing=true) {

  // Coersing a NumericVector into a cube for ease of use
//   IntegerVector dims=graph.attr("dim");
//   const arma::cube graph_cube(graph.begin(), dims[0], dims[1], dims[2], false);

  // Variables initialization
//   const int n = graph_cube.n_rows;
//   const int T = graph_cube.n_slices;
  NumericVector suscep(n);

  // Variables to use within loop
  int ti, tj;

  // Creating discount variable, firts, must be truncated.
  if (K >= graph.size()) {
    warning("Too many periods selected, will be truncated to T-1.");
    K = graph.size() - 1;
  }

  // The discount can be either exponential (1 + r)^(k-1), or
  // lineal in the form of k.
  // double * discount = new double[K];
  std::vector< double > discount(K);
  if (!expdiscount) {
    for(int k=1;k<=K;++k)
      discount[k-1] = k;
  }
  else {
    for(int k=1;k<=K;++k)
      discount[k-1] = pow((1.0 + r), k-1.0);
  }

  for(int i=0;i<n;++i) {
    // Checling user interrup
    if (i % 1000 == 0)
      Rcpp::checkUserInterrupt();

    // If NA (aka nan in Armadillo), then NA.
    if (!arma::is_finite(times(i))) {
      suscep.at(i) = NA_REAL;
      continue;
    }

    // Capturing variables
    ti = times(i);
    if (ti == 1) continue;

    double numerator = 0.0;
    double denominator = 0.0;

    // For the adjusted verion, see the mathematical supplement on Valente et al. (2015)
    double nadopt_t = 0;

    for(int k=1;k<=K;++k) {

      // Current time period from 1 to T
      int t = ti - k + 1;

      // If the required time period does not exists, then continue, recall that
      // vectors can be reached starting 0. If t=1, then in ++C it is equivalent
      // to 0, so we need to reach time of adoption at -1 (which does not
      // exists, or we don't know if it exists).
      if (t <= 1) continue;

      // Storing the graph
      arma::sp_mat graph_cube = graph.at(t - 1);
      if (!valued)   graph_cube = arma::spones(graph_cube);
      if (!outgoing) graph_cube = graph_cube.t();

      for(int j=0;j<n;++j) {
        if (i==j) continue;

        tj = times(j);
        // Adding up for t+k iff a link between j and i exists. Notice that the
        // t - 1 is because t is in [1;T], and we actually want t-1. If t=1, then
        // in ++C it is 0.
        // if (graph_cube(i,j,t - 1) != 0) {
        if (graph_cube(i,j) != 0) {
          if (tj == (ti - k)) {
            numerator += 1.0 / discount[k-1];
          }
          // Adding up for t+1 <= t <= T
          if (tj <=(ti - k)) denominator += 1.0 / discount[k-1];
        }

        // Has adopted so far? (discounted version)
        if (tj == (ti - k)) nadopt_t += 1.0 / discount[k-1];
      }
    }

    // Putting all together
    suscep(i) = (numerator / (denominator + 1e-15)) ;
    if (normalize) suscep(i) = suscep(i) / (nadopt_t + 1e-15);
  }

  return suscep;
}

/***R
library(netdiffuseR)
set.seed(123)
graph <- rand_graph(undirected = FALSE)
graph <- array(unlist(lapply(1:3, function(x,...) graph)), dim=c(10,10,3))
times <- sample(1:3, 10, TRUE)
infection(graph, times)
plot_diffnet(graph, toa_mat(times)$cumadopt)

# Paper data
graph <- edgelist_to_adjmat(matrix(c(5,1,5,2,5,3,5,4,6,5,7,5,8,5), byrow = TRUE, ncol=2), simplify = FALSE)
graph <- array(unlist(lapply(1:5, function(x,...) graph)), dim=c(8,8,5))
times <- c(1,2,5,5,3,4,5,5)
plot_diffnet(graph, toa_mat(times)$cumadopt)
infection(graph, times, TRUE)
susceptibility(graph, times, TRUE)
*/


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
  for(int i = 0; i < n; ++i)
    for(int j = 0; j < n; ++j) {
      double chg = adjmat_t1(i, j) - adjmat_t0(i, j);
      if      (chg > 0) change_mat(i, j) = 1;
      else if (chg < 0) change_mat(i, j) = -1;
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

  for(int i = 0; i < n; ++i) {

    for(int j = 0; j < n; ++j) {

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
