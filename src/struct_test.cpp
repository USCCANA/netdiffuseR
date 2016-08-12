// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

typedef double (*funcPtr)(double y0, double y1);

double st_dist(double y0, double y1) {return fabs(y0-y1);}
double st_greater(double y0, double y1) {return (double) (y0 > y1);}
double st_greaterequal(double y0, double y1) {return (double) (y0 >= y1);}
double st_smaller(double y0, double y1) {return (double) (y0 < y1);}
double st_smallerequal(double y0, double y1) {return (double) (y0 <= y1);}
double st_equal(double y0, double y1) {return (double) (y0 == y1);}

// XPtr<funcPtr> st_getfun(std::string funname) {
void st_getfun(std::string funname, funcPtr & fun) {
  if      (funname == "distance")                           fun = &st_dist;
  else if ((funname == "greater") | (funname == ">"))       fun = &st_greater;
  else if ((funname == "greaterequal") | (funname == ">=")) fun =  &st_greaterequal;
  else if ((funname == "smaller") | (funname == "<"))       fun =  &st_smaller;
  else if ((funname == "smallerequal") | (funname == "<=")) fun =  &st_smallerequal;
  else if ((funname == "equal") | (funname == "=="))        fun =  &st_equal;
  else Rcpp::stop("Unkown function.");

  return ;
}

// [[Rcpp::export]]
double struct_test_mean(NumericVector & y,
                        std::string funname, bool self=false) {

  int    n =y.size();
  double m = (self? n*n/2:n*(n-1)/2);
  double ans = 0.0;

  // Fetching function
  funcPtr fun;
  st_getfun(funname, fun);

  for (int i=0;i<n;i++)
    for (int j=i;j<n;j++) {
      if (self) ans+= fun(y[i], y[j])/m;
      else if (i!=j) ans+= fun(y[i], y[j])/m;
    }

    return ans;
}


// [[Rcpp::export]]
double struct_test_var(NumericVector & y, std::string funname, bool self=false)  {

  int n      = y.size();
  double m   = (self ? n*n*n : n*(n-1)*(n-1));
  double ans = 0.0;

  // Fetching function
  funcPtr fun;
  st_getfun(funname, fun);

  for (int i=0;i<n;i++)
    for (int j=0;j<n;j++)
      for (int k=0;k<n;k++) {
        if (self) ans+=fun(y[i],y[j])*fun(y[i],y[k])/m;
        else if (i!=j && i!=k) ans+=fun(y[i],y[j])*fun(y[i],y[k])/m;
      }

      return ans - pow(struct_test_mean(y, funname, self),2.0);
}



/***R
set.seed(1231)
Y <- rnorm(500)

# Greater
struct_test_mean(Y, "greater")
struct_test_var(Y, "greater")

# Diff
struct_test_mean(Y, "distance")
struct_test_var(Y, "distance")

# Equal
struct_test_mean(round(Y,1), "equal")
struct_test_var(round(Y,1), "equal")
#
# L <- letters[1:3]
# for(i in 1:3)
#   for(j in i:3)
#     for (k in i:3) {
#
#     }

library(netdiffuseR)
x     <- kfamilyDiffNet
index <- which(x$toa != max(x$toa, na.rm = TRUE))
x     <- x[index,]
Mobs <- mean(threshold(x), na.rm = TRUE)
# Mobs <- sum(threshold(x)*dgr(x, cmode = "outdegree")[,1]/sum(dgr(x, cmode = "outdegree")[,1]), na.rm = TRUE)
M <- struct_test_mean(x$toa, "smaller")
S <- struct_test_var(x$toa, "smaller")

pval <- pnorm(sqrt(nnodes(x))*(M-Mobs)/sqrt(S))
ifelse(pval >.5,1-pval,pval)*2

d <- dgr(x)[,1]
S <- sapply(sample(1:nnodes(x), 1e4, TRUE), function(i) {
  mean(sample(x$toa, d[i]) < x$toa[i])
})

mean(S, na.rm = TRUE)
var(S, na.rm = TRUE)
pval <- pnorm(sqrt(nnodes(x))*(Mobs-M)/sd(S))
pval <- ifelse(pval > .5, 1-pval, pval)*2
pval
# struct_test(x, function(x) mean(threshold(x), na.rm = TRUE), R=500)

x <- matrix(rnorm(5e5*3), ncol=3)

mean((x[,1] > x[,2])*(x[,1] > x[,3]))
mean(x[,1]>x[,2])
*/
