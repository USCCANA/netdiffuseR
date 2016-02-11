/*******************************************************************************
 adjmat.cpp: FUNCTIONS FOR NETWORK DATA MANAGEMENT

 created: nov 2, 2015
 copyright: MIT + LICENCE file

 This set of functions helps handling network data, in particular, diffussion
 network data. Importing edgelist as adjacency matrices and generating auxiliary
 matrices used for the model.

*******************************************************************************/

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

void is_square_sp_mat(const arma::sp_mat & mat) {
  if (mat.n_cols != mat.n_rows)
    stop("It must be a square matrix");
  return;
}

// Returns min from an integer vector accounting for NA objects
double min_int_na_cpp(const IntegerVector & x, const LogicalVector & isna) {
  int n =x.size();
  double min = DBL_MAX;

  for(int i=0;i<n;i++) {
    if (isna[i]) continue;
    if (min > x[i]) min = x[i];
  }
  return min;
}

// Returns max from an integer vector accounting for NA objects
double max_int_na_cpp(const IntegerVector & x, const LogicalVector & isna) {
  int n =x.size();
  double max = -DBL_MAX;

  for(int i=0;i<n;i++) {
    if (isna[i]) continue;
    if (max < x[i]) max = x[i];
  }
  return max;
}

// [[Rcpp::export]]
List toa_mat_cpp(const IntegerVector & year, int t0, int t1) {

  // Pin down NAs
  LogicalVector isna = is_na(year);

  // Measuring time
  int n = year.size();

  // Creating output
  List out(2);
  arma::mat adopt(n,t1 - t0 + 1, arma::fill::zeros);

  for(int i=0;i<n;i++) {
    if (isna[i]) continue;
    adopt(i,year[i]-t0) = 1.0;
  }

  arma::mat cumadopt = cumsum(adopt, 1);

  /* Adopt_mat -> cumadopt
    Adopt_mat1 -> adopt
  */
  return List::create(_["adopt"]=adopt, _["cumadopt"]=cumadopt);
}

// [[Rcpp::export]]
arma::sp_mat edgelist_to_adjmat_cpp(
    const arma::mat & edgelist,
    NumericVector weights = NumericVector::create(),
    int n = 0,
    bool undirected = false,
    bool self = false,
    bool multiple = false) {

  // Getting the dimensions and creating objects
  int m = edgelist.n_rows;

  // Checking out weights
  NumericVector w(m,1.0);
  if (weights.size() == m) w = clone(weights);

  // Identifying the unique nodes and creating the adjmat (base)
  // If n is not provided (0), then we need to calculate the number of
  // nodes in the network.
  if (n == 0) {
    n = (int) max(max(edgelist));
  }

  // Creating output mat and aux mat. The counts mat is used in logical expressions
  // since the arma::mat class can lead to error when compared to 0 (because of
  // the precision)
  arma::sp_mat adjmat(n,n);

  for(int i=0;i<m;i++) {
    // Ids of the vertices
    int ego   = ((int) edgelist(i,0)) - 1;
    int alter = ((int) edgelist(i,1)) - 1;

    // If self edges are not allowed
    if (!self && (ego == alter)) continue;

    // If undirected, the order does not matters
    if (undirected) {
      int tmp=ego;
      if (ego > alter) ego = alter, alter=tmp;
    }

    // If multiple edges are not allowed.
    if (!multiple && (adjmat(ego, alter) != 0)) continue;

    adjmat.at(ego, alter) += w[i];

    // If undirected, must include in the switch
    if (undirected) adjmat.at(alter, ego) += w[i];
  }

  return adjmat;
}


// [[Rcpp::export]]
arma::mat adjmat_to_edgelist_cpp(const arma::sp_mat & adjmat, bool undirected = true) {
  std::vector< double > ego;
  std::vector< double > alter;

  int n = adjmat.n_cols;

  // CAN BE REWRITTEN ACCESSING VALUES OF THE SPMAT DIRECTLY

  for(int i=0;i<n;i++) {
    /* Setting the length of the subloop acordingly to type of graph */
    int m = n;
    if (undirected) m=i;
    for(int j=0;j<m;j++) {
      if (adjmat(i,j))
        ego.push_back(i+1.0), alter.push_back(j+1.0);
    }
  }

  // Creating colvectors to be used with join_rows.
  arma::mat egom(ego);
  arma::mat alterm(alter);

  arma::mat edgelist = join_rows(egom, alterm);

  return edgelist;
}


// // [[Rcpp::export]]
// arma::mat adjmat_to_dyn_edgelist_cpp(NumericVector adjmat, bool undirected=true) {
//
//   // Coersing a NumericVector into a cube for ease of use
//   IntegerVector dims=adjmat.attr("dim");
//   const arma::cube adjmat_cube(adjmat.begin(), dims[0], dims[1], dims[2], false);
//
//   int T = adjmat_cube.n_slices;
//   int n = adjmat_cube.n_cols;
//
//   std::vector< double > ego;
//   std::vector< double > alter;
//   std::vector< double > time;
//
//   for(int t=0;t<T;t++)
//     for(int i=0;i<n;i++) {
//       /* Setting the length of the subloop acordingly to type of graph */
//       int m = n;
//       if (undirected) m=i;
//       for(int j=0;j<m;j++)
//         if (adjmat_cube(i,j,t))
//           ego.push_back(i+1.0), alter.push_back(j+1.0), time.push_back(t+1.0);
//     }
//
//   // Creating colvectors to be used with join_rows.
//   arma::mat egom(ego);
//   arma::mat alterm(alter);
//   arma::mat timem(time);
//
//   arma::mat edgelist = join_rows(join_rows(egom, alterm), timem);
//
//   return edgelist;
// }


// [[Rcpp::export]]
IntegerMatrix toa_diff_cpp(const IntegerVector & year) {
  int n = year.size();

  IntegerMatrix diff(n,n);
  LogicalVector isna = is_na(year);

  for(int i=0;i<n;i++) {
    // Checling user interrup
    if (i % 1000 == 0)
      Rcpp::checkUserInterrupt();

    // If na, then fill the diff with NA
    if (isna[i]) {
      for(int j=0;j<n;j++)
        diff(i,j) = NA_INTEGER, diff(j,i) = NA_INTEGER;
      continue;
    }
    for(int j=0;j<i;j++) {
      if (isna[j]) continue;
      diff(i,j) = year[j]-year[i], diff(j,i)=year[i]-year[j];
    }
  }

  return diff;
}

// [[Rcpp::export]]
arma::icolvec isolated_cpp(
    const arma::sp_mat & adjmat,
    bool undirected=true) {

  int n = adjmat.n_cols;

  // Checking size
  is_square_sp_mat(adjmat);

  arma::icolvec isolated(n, arma::fill::ones);

  // Looping through (all) the matrix. Setting the value to 0
  // whenever there's a link (for both individuals).

  for(int i=0;i<n;i++) {
    // First, need to get the size of the loop
    int m=n;
    if (undirected) m = i;
    for(int j=0;j<m;j++)
      if ((adjmat.at(i,j)!=0))
        isolated[i]=0, isolated[j]=0;
  }

  return isolated;
}

// [[Rcpp::export]]
arma::sp_mat drop_isolated_cpp(
    const arma::sp_mat & adjmat,
    arma::icolvec isolated, bool undirected=true) {

  int n = adjmat.n_cols;

  // Checking size
  is_square_sp_mat(adjmat);

  if (isolated.n_rows==0) isolated = isolated_cpp(adjmat, undirected);
  else if (isolated.n_rows != n) stop("-isolated- must have the same length as nrow(adjmat)");

  int m = sum(isolated);
  arma::sp_mat newadjmat(n-m,n-m);

  // Rcpp::warning()

  // Indexes of the new adjacency matrix
  int ii=0;
  int ji=0;
  for(int i=0;i<n;i++) {

    // If an isolated was found, continue next
    if (isolated[i]) continue;
    for(int j=0;j<n;j++) {
      if (isolated[j]) continue;
      newadjmat.at(ii,ji++)=adjmat.at(i,j);
    }
    // Continue next
    ii++,ji=0;
  }

  return newadjmat;
}

/** *R
set.seed(123)
lonenet <- rand_graph_cpp(8, p=.8,undirected=TRUE)
lonenet[c(1,4),] <- 0
lonenet[,c(1,4)] <- 0
isolated_cpp(lonenet)
*/
