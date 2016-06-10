// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
List egonet_attrs_cpp(
    const arma::sp_mat & graph, const arma::uvec V,
    bool outer=true, bool self=true, bool self_attrs=false, bool valued=true) {

  // General variables
  int N = V.n_elem;

  // Column names
  CharacterVector cnames(2);

  cnames[0] = "value";
  cnames[1] = "id";

  // Depending on inner or outer edges:
  //  - outer: Accounts for the rows
  //  - inner: Accounts for the columns
  arma::sp_mat tgraph;
  if (outer) tgraph = graph.t();
  else tgraph = graph;

  List data(N);

  for (int v=0;v<N;v++) {
    // Index
    int e = V.at(v);

    // Analyzing the case when is undirected
    int rm = 0;
    if (graph.at(e,e) && !self) rm = 1;

    // Individual specific variables
    arma::sp_mat g = tgraph.col(e);
    NumericMatrix out(g.n_nonzero-rm + (self_attrs ? 1 : 0),2);
    // We add two so:
    //  - we can include the value of the edge
    //  - we cna include the id number of the vertex (wich goes from 1 to n)

    // Retrieving the desired set of attributes
    if (self_attrs) {
      out(0,0) = tgraph.at(v,v);
      out(0,1) = v + 1;
    }
    int nloop = (self_attrs ? 1 : 0);
    for (unsigned i = 0u;i < g.n_nonzero;i++) {

      // Edge index
      int index = g.row_indices[i];
      if (!self && index == e) continue;
      out(nloop,1) = index + 1;

      // Edge value
      if (valued) out(nloop,0) = g.values[i];
      else out(nloop,0) = 1.0;

      // Increasing after success
      nloop++;
    }

    colnames(out) = cnames;

    data[v] = out;
  }


  return data;
}

/***R
library(netdiffuseR)
set.seed(123)
n <- 1e3
graph <- rgraph_ba(t=n-1)
mat <- cbind(id = 1:n, rand=runif(n), d=dgr(graph))

diffnet_attrs <- function(diffnet) {
  attrs
}



onestep_attrs <- egonet_attrs_cpp(graph,1:n-1, self = FALSE)
# twostep_attrs <- nsteps_attrs_cpp(graph %*% graph,1:n, mat, self = FALSE)
stop()
# Computing the degree weighted rand
w <- lapply(onestep_attrs, function(x) {
    (t(x[,4]) %*% x[,5])/sum(x[,5])
  })

w <- do.call(rbind, w)
*/

/***R
library(sna)

g <- matrix(0, ncol=6, nrow=6)
g[2:3,1] <- g[1,2:3] <- 1
g[4,2] <- g[2,4] <-1
g[6,5] <- g[5,6] <- 1

oldpar <- par(no.readonly = TRUE)
coord <- gplot.layout.fruchtermanreingold(g, list())
par(mfrow=c(2,2))
gplot(g, displaylabels = TRUE, coord = coord, edge.lwd = 1)
gplot(g %*% g, displaylabels = TRUE, coord = coord, edge.lwd = 1)
gplot(g %*% g %*% g, displaylabels = TRUE, coord = coord, edge.lwd = 2)
gplot(g %*% g %*% g %*% g, displaylabels = TRUE, coord = coord, edge.lwd = 1)
par(oldpar)

*/
