% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mentor.r
\name{mentor_matching}
\alias{mentor_matching}
\alias{leader_matching}
\alias{plot.diffnet_mentor}
\title{Optimal Leader/Mentor Matching}
\usage{
mentor_matching(
  graph,
  n,
  cmode = "indegree",
  lead.ties.method = "average",
  geodist.args = list()
)

leader_matching(
  graph,
  n,
  cmode = "indegree",
  lead.ties.method = "average",
  geodist.args = list()
)

\method{plot}{diffnet_mentor}(
  x,
  y = NULL,
  vertex.size = "degree",
  minmax.relative.size = getOption("diffnet.minmax.relative.size", c(0.01, 0.04)),
  lead.cols = grDevices::topo.colors(attr(x, "nleaders")),
  vshapes = c(Leader = "square", Follower = "circle"),
  add.legend = TRUE,
  main = "Mentoring Network",
  ...
)
}
\arguments{
\item{graph}{Any class of accepted graph format (see \code{\link{netdiffuseR-graphs}}).}

\item{n}{Number of leaders}

\item{cmode}{Passed to \code{\link{dgr}}.}

\item{lead.ties.method}{Passed to \code{\link{rank}}}

\item{geodist.args}{Passed to \code{\link{approx_geodesic}}.}

\item{x}{An object of class \code{diffnet_mentor}.}

\item{y}{Ignored.}

\item{vertex.size}{Either a numeric scalar or vector of size \eqn{n}, or any
of the following values: "indegree", "degree", or "outdegree" (see details).}

\item{minmax.relative.size}{Passed to \code{\link{rescale_vertex_igraph}}.}

\item{lead.cols}{Character vector of length \code{attr(x,"nleaders")}. Colors
to be applied to each group. (see details)}

\item{vshapes}{Character scalar of length 2. Shapes to identify leaders (mentors)
and followers respectively.}

\item{add.legend}{Logical scalar. When \code{TRUE} generates a legend to distinguish
between leaders and followers.}

\item{main}{Character scalar. Passed to \code{\link[graphics:title]{title}}}

\item{...}{Further arguments passed to \code{\link[igraph:plot.igraph]{plot.igraph}}}
}
\value{
An object of class \code{diffnet_mentor} and \code{data.frame} with the following columns:
\item{name}{Character. Labels of the vertices}
\item{degree}{Numeric. Degree of each vertex in the graph}
\item{iselader}{Logical. \code{TRUE} when the vertex was picked as a leader.}
\item{match}{Character. The corresponding matched leader.}

The object also contains the following attributes:

\item{nleaders}{Integer scalar. The resulting number of leaders (could be greater than \code{n})}.
\item{graph}{The original graph used to run the algorithm.}
}
\description{
Implementes the algorithm described in Valente and Davis (1999)
}
\details{
The algorithm works as follows:
\enumerate{
  \item Find the top \code{n} individuals ranking them by \code{dgr(graph, cmode)}.
  The rank is computed by the function \code{\link{rank}}. Denote this set \code{M}.
  \item Compute the geodesic matrix.
  \item For each \code{v in V} do:

  \enumerate{
    \item Find the mentor \code{m in M} such that is closest to \code{v}
    \item Were there a tie, choose the mentor that minimizes the average
    path length from \code{v}'s direct neighbors to \code{m}.
    \item If there are no paths to any member of \code{M}, or all have the
    same average path length to \code{v}'s neighbors, then assign one
    randomly.
  }
}

Plotting is done via the function \code{\link[igraph:plot.igraph]{plot.igraph}}.

When \code{vertex.size} is either of \code{"degree"}, \code{"indegree"}, or
\code{"outdegree"}, \code{vertex.size} will be replace with \code{dgr(.,cmode = )}
so that the vertex size reflects the desired degree.

The argument \code{minmax.relative.size} is passed to \code{\link{rescale_vertex_igraph}}
which adjusts \code{vertex.size} so that the largest and smallest vertices
have a relative size of \code{minmax.relative.size[2]} and
\code{minmax.relative.size[1]} respectively with respect to the x-axis.
}
\examples{
# A simple example ----------------------------------------------------------
set.seed(1231)
graph <- rgraph_ws(n=50, k = 4, p = .5)

# Looking for 3 mentors
ans <- mentor_matching(graph, n = 3)

head(ans)
table(ans$match) # We actually got 9 b/c of ties

# Visualizing the mentor network
plot(ans)

}
\references{
Valente, T. W., & Davis, R. L. (1999). Accelerating the Diffusion of
Innovations Using Opinion Leaders. The ANNALS of the American Academy of
Political and Social Science, 566(1), 55–67.
\doi{10.1177/000271629956600105}
}
