#include <Rcpp.h>

using namespace Rcpp;

// Combine two vectors and (if) get the unique vect
NumericVector vec_comb(const NumericVector & a, const NumericVector & b, bool uni = false) {

  // Parameters
  int n_a = a.size(), n_b = b.size();

  NumericVector c(n_a+n_b);

  for(int i=0;i<n_a;i++)
    c[i] = a[i];

  for(int i=0;i<n_b;i++)
    c[i+n_a] = b[i];

  // If returning a unique vector
  if (uni) return(unique(c));
  return c;
}

// [[Rcpp::export]]
List adopt_mat_cpp(const IntegerVector & year) {

  // Measuring time
  int T0 = min(year);
  int T = max(year) - T0 + 1;
  int n = year.size();

  // Creating output
  List out(2);
  IntegerMatrix adoptmat1(n,T);

  for(int i=0;i<n;i++)
    adoptmat1(i,year[i]-T0) = 1;

  IntegerMatrix adoptmat = clone(adoptmat1);
  for(int i=0;i<n;i++)
    for(int j=1;j<T;j++)
      adoptmat(i,j) += adoptmat(i,j-1);

  return List::create(_["adoptmat"]=adoptmat,_["adoptmat1"]=adoptmat1);
}

// [[Rcpp::export]]
NumericMatrix edgelist_to_adjmat_cpp(
    const NumericMatrix & data,
    NumericVector weights = NumericVector::create(),
    int n = 0,
    bool undirected = false) {

  // Getting the dimensions and creating objects
  int m = data.nrow();

  // Checking out weights
  NumericVector w(m,1.0);
  if (weights.size() == m) w = clone(weights);

  // Identifying the unique nodes and creating the adjmat (base)
  if (n == 0) {
    NumericVector nodes = vec_comb(data(_,0), data(_,1), true);
    n = max(nodes);
  }
  NumericMatrix mat(n,n);

  for(int i=0;i<m;i++) {
    mat(data(i,0)-1,data(i,1)-1) += w[i];
    if (undirected) mat(data(i,1)-1,data(i,0)-1) += w[i];
  }

  return mat;
}

// [[Rcpp::export]]
IntegerMatrix toa_mat_cpp(const IntegerVector & year) {
  int n = year.size();

  IntegerMatrix toa(n,n);

  for(int i=0;i<n;i++)
    for(int j=0;j<i;j++)
      toa(i,j) = year[j]-year[i], toa(j,i)=year[i]-year[j];

  return toa;
}

// [[Rcpp::export]]
IntegerVector isolated_cpp(const NumericMatrix & adjmat) {

  int n = adjmat.ncol();
  IntegerVector isolated(n,1);

  // Looping through (all) the matrix. Setting the value to 0
  // whenever there's a link (for both individuals).
  for(int i=0;i<n;i++)
    for(int j=0;j<n;j++)
      if (adjmat(i,j))
        isolated[i] = 0, isolated[j] = 0;

  return isolated;
}

