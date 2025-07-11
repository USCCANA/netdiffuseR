% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/random_graph.R
\name{rgraph_ba}
\alias{rgraph_ba}
\alias{scale-free}
\title{Scale-free and Homophilic Random Networks}
\usage{
rgraph_ba(m0 = 1L, m = 1L, t = 10L, graph = NULL, self = TRUE, eta = NULL)
}
\arguments{
\item{m0}{Integer scalar. Number of initial vertices in the graph.}

\item{m}{Integer scalar. Number of new edges per vertex added.}

\item{t}{Integer scalar. Number of time periods (steps).}

\item{graph}{Any class of accepted graph format (see \code{\link{netdiffuseR-graphs}}).}

\item{self}{Logical scalar. When \code{TRUE} autolinks (loops, self edges) are allowed (see details).}

\item{eta}{Numeric vector of length \code{t+m0}. When specified, it generates
a scale-free homophilic network (see details).}
}
\value{
If \code{graph} is not provided, a static graph, otherwise an expanded
graph (\code{t} aditional vertices) of the same class as \code{graph}.

The resulting graph will have \code{graph$meta$undirected = FALSE} if it is of
class \code{diffnet} and \code{attr(graph, "undirected")=FALSE} otherwise.
}
\description{
Generates a scale-free random graph based on Bollabas et al. (2001), also know as
\emph{Linearized Chord Diagram} (LCD) which has nice mathematical propoerties.
And also scale-free homophilic networks when an vertex attribute \code{eta} is
passed.
}
\details{
Based on Ballobás et al. (2001) creates a directed random graph of size
\code{t + m0}. A big difference with B-A model
is that this allows for loops (self/auto edges) and further multiple links,
nevertheless, as \eqn{t} increases, the number of such cases reduces.

By default, the degree of the first \code{m0} vertices is set to be 2 (loops).
When \code{m>1}, as described in the paper, each new link from the new vertex
is added one at a time
\dQuote{counting \sQuote{outward half} of the edge being added as already contributing to the degrees}.

When \code{self=FALSE}, the generated graph is created without autolinks. This
means that at the beginning, if the number of links equals zero, all vertices
have the same probability of receiving a new link.

When \code{eta} is passed, it implements the model specified in De Almeida et al.
(2013), a scale-free homophilic network. To do so \code{eta} is rescaled to
be between 0 and 1 and the probability that the node \eqn{i} links to node \eqn{j}
is as follows:

\deqn{
\frac{(1-A_{ij})k_j}{\sum_j (1-A_{ij})k_j}
}{
[(1-A(ij))k(j)]/\sum_j[(1-A(ij))k(j)]
}

Where \eqn{A_{ij} = |\eta_i - \eta_j|}{A(ij) = |eta(i) - eta(j)|} and
\eqn{k_j}{k(j)} is the degree of the \eqn{j}-th vertex.
}
\examples{
# Using another graph as a base graph ---------------------------------------
graph <- rgraph_ba()
graph

graph <- rgraph_ba(graph=graph)

# Generating a scale-free homophilic graph (no loops) -----------------------
set.seed(112)
eta <- rep(c(1,1,1,1,2,2,2,2), 20)
ans <- rgraph_ba(t=length(eta) - 1, m=3, self=FALSE, eta=eta)

# Converting it to igraph (so we can plot it)
ig  <- igraph::graph_from_adjacency_matrix(ans)

# Neat plot showing the output
oldpar <- par(no.readonly = TRUE)
par(mfrow=c(1,2))
plot(ig, vertex.color=c("red","blue")[factor(eta)], vertex.label=NA,
    vertex.size=5, main="Scale-free homophilic graph")
suppressWarnings(plot(dgr(ans), main="Degree distribution"))
par(oldpar)

}
\references{
Bollobás, B´., Riordan, O., Spencer, J., & Tusnády, G. (2001). The degree
sequence of a scale-free random graph process. Random Structures & Algorithms,
18(3), 279–290. \doi{10.1002/rsa.1009}

Albert-László Barabási, & Albert, R. (1999). Emergence of Scaling in Random
Networks. Science, 286(5439), 509–512. \doi{10.1126/science.286.5439.509}

Albert-László Barabási. (2016). Network Science: (1st ed.). Cambridge University Press.
Retrieved from \url{https://networksciencebook.com}

De Almeida, M. L., Mendes, G. A., Madras Viswanathan, G., & Da Silva, L. R. (2013).
Scale-free homophilic network. European Physical Journal B, 86(2).
\doi{10.1140/epjb/e2012-30802-x}
}
\seealso{
Other simulation functions: 
\code{\link{permute_graph}()},
\code{\link{rdiffnet}()},
\code{\link{rewire_graph}()},
\code{\link{rgraph_er}()},
\code{\link{rgraph_ws}()},
\code{\link{ring_lattice}()}
}
\author{
George G. Vega Yon
}
\concept{Barabasi-Albert model}
\concept{Bollabas}
\concept{Homophilic random graph}
\concept{Random graph}
\concept{Scale-free random graph}
\concept{simulation functions}
\keyword{distribution}
