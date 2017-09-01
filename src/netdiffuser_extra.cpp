#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::sp_mat sp_trimatl(const arma::sp_mat & x) {
  // Getting start-end
  arma::sp_mat::const_iterator start = x.begin();
  arma::sp_mat::const_iterator end   = x.end();

  // Empty mat
  int n = x.n_cols;
  arma::sp_mat newx(n,n);

  for (arma::sp_mat::const_iterator it = start; it != end; ++it) {
    int i = it.row();
    int j = it.col();
    if (i >= j)
      newx.at(i,j) = *it;
  }

  return newx;
}

// [[Rcpp::export]]
arma::sp_mat sp_diag(const arma::sp_mat & x, const arma::vec & v) {
  // Checking dimensions
  if (x.n_cols != x.n_rows) Rcpp::stop("-x- must be square.");
  if (x.n_cols != v.n_elem) Rcpp::stop("length(v) must be equal to ncol(x)");

  arma::sp_mat out(x);
  out.diag() = v;

  return out;
}


// [[Rcpp::export]]
int unif_rand_w_exclusion(int n, int e) {
  double r = unif_rand();
  int num = floor(r*(n-1));
  if (e != n) {
    if (num >= e) ++num;
  }
  return(num);
}


// [[Rcpp::export]]
arma::sp_mat sp_as_undirected(const arma::sp_mat & x) {
  // Getting start-end
  arma::sp_mat::const_iterator start = x.begin();
  arma::sp_mat::const_iterator end   = x.end();

  // Empty mat
  arma::sp_mat newx(x);

  for (arma::sp_mat::const_iterator it = start; it != end; ++it) {
    int i = it.row();
    int j = it.col();
    newx.at(j,i) = *it;
  }

  return newx;
}

typedef double (*funcPtr)(double y0, double y1);

double st_dist(double y0, double y1) {return fabs(y0-y1);}
double st_quaddist(double y0, double y1) {return pow(y0-y1, 2.0);}
double st_greater(double y0, double y1) {return (double) (y0 > y1);}
double st_greaterequal(double y0, double y1) {return (double) (y0 >= y1);}
double st_smaller(double y0, double y1) {return (double) (y0 < y1);}
double st_smallerequal(double y0, double y1) {return (double) (y0 <= y1);}
double st_equal(double y0, double y1) {return (double) (y0 == y1);}

double st_min(double y0, double y1) {return fmin(y0, y1);}
double st_max(double y0, double y1) {return fmax(y0, y1);}
double st_mean(double y0, double y1) {return (y0 + y1) /2.0;}


// XPtr<funcPtr> st_getfun(std::string funname) {
void st_getfun(std::string funname, funcPtr & fun) {
  if      (funname == "distance")                           fun = &st_dist;
  else if ((funname == "quaddist") | (funname == "^2"))       fun = &st_quaddist;
  else if ((funname == "greater") | (funname == ">"))       fun = &st_greater;
  else if ((funname == "greaterequal") | (funname == ">=")) fun =  &st_greaterequal;
  else if ((funname == "smaller") | (funname == "<"))       fun =  &st_smaller;
  else if ((funname == "smallerequal") | (funname == "<=")) fun =  &st_smallerequal;
  else if ((funname == "equal") | (funname == "=="))        fun =  &st_equal;
  else if ((funname == "min") | (funname == "minimum"))        fun =  &st_min;
  else if ((funname == "max") | (funname == "maximum"))        fun =  &st_max;
  else if ((funname == "mean") | (funname == "avg"))        fun =  &st_mean;
  else Rcpp::stop("Unkown function.");

  return ;
}

// [[Rcpp::export]]
arma::sp_mat bootnet_fillself(
    arma::sp_mat & graph,
    const IntegerVector & index,
    const NumericVector & E
  ) {

  // Parameters: Size and density of the graph
  int n = index.length();
  int m = E.length();
  double dens = ((double) m) / (double) (n*n);

  // Finding repeated values
  std::vector< std::vector<int> > reps(n);
  for (int i=0; i<n; i++) {
    reps.at(index.at(i)-1).push_back(i);
  }

  // Sampling E
  NumericVector rand(2);
  for (int i=0; i<n; i++) {
    // If has no elements
    if (reps.at(i).size() < 2) continue;

    std::vector<int> r(reps.at(i));

    for (unsigned int j=0; j< r.size(); j++)
      for (unsigned int k=j; k<r.size(); k++) {

        if (j==k) continue;

        // Add accordingly to density
        if (unif_rand() <= dens)
          graph.at(r.at(j),r.at(k)) = E.at(floor(unif_rand()*m));

        // Add accordingly to density
        if (unif_rand() <= dens)
          graph.at(r.at(k),r.at(j)) = E.at(floor(unif_rand()*m));
      }

  }

  return graph;
}
