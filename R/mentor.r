#' Optimal Leader/Mentor Matching
#'
#' Implementes the algorithm described in Valente and Davis (1999)
#'
#' @template graph_template
#' @param n Number of leaders
#' @param cmode Passed to \code{\link{dgr}}.
#' @param lead.ties.method Passed to \code{\link{rank}}
#' @param geodist.args Passed to \code{\link{approx_geodesic}}.
#'
#' @details The algorithm works as follows:
#' \enumerate{
#'   \item Find the top \code{n} individuals ranking them by \code{dgr(graph, cmode)}.
#'   The rank is computed by the function \code{\link{rank}}. Denote this set \code{M}.
#'   \item Compute the geodesic matrix.
#'   \item For each \code{v in V} do:
#'
#'   \enumerate{
#'     \item Find the mentor \code{m in M} such that is closest to \code{v}
#'     \item Were there a tie, choose the mentor that minimizes the average
#'     path length from \code{v}'s direct neighbors to \code{m}.
#'     \item If there are no paths to any member of \code{M}, or all have the
#'     same average path length to \code{v}'s neighbors, then assign one
#'     randomly.
#'   }
#' }
#'
#' @return An object of class \code{diffnet_mentor} and \code{data.frame} with the following columns:
#' \item{name}{Character. Labels of the vertices}
#' \item{degree}{Numeric. Degree of each vertex in the graph}
#' \item{iselader}{Logical. \code{TRUE} when the vertex was picked as a leader.}
#' \item{match}{Character. The corresponding matched leader.}
#'
#' The object also contains the following attributes:
#'
#' \item{nleaders}{Integer scalar. The resulting number of leaders (could be greater than \code{n})}.
#' \item{graph}{The original graph used to run the algorithm.}
#'
#' @references
#' Valente, T. W., & Davis, R. L. (1999). Accelerating the Diffusion of
#' Innovations Using Opinion Leaders. The ANNALS of the American Academy of
#' Political and Social Science, 566(1), 55–67.
#' \url{http://journals.sagepub.com/doi/abs/10.1177/000271629956600105}
#' @examples
#' # A simple example ----------------------------------------------------------
#' set.seed(1231)
#' graph <- rgraph_ws(n=50, k = 4, p = .5)
#'
#' # Looking for 3 mentors
#' ans <- mentor_matching(graph, n = 3)
#'
#' head(ans)
#' table(ans$match) # We actually got 9 b/c of ties
#'
#' # Visualizing the mentor network
#' plot(ans)
#'
#' @export
mentor_matching <- function(
  graph,
  n,
  cmode            = "indegree",
  lead.ties.method = "average",
  geodist.args     = list()
) {

  cls <- class(graph)

  if (any(c("dgCMatrix", "matrix", "array") %in% cls)) {
    # Adding labels in case there aren't any
    if (!length(rownames(graph))) {
      warning("-graph- as no labels (rownames). We'll add some from 1 to n.")
      dimnames(graph) <- list(1:nnodes(graph), 1:nnodes(graph))
    }
  }

  if (any(c("dgCMatrix", "matrix") %in% cls)) {

    # Matrix method
    .mentor_matching(graph, n, cmode, lead.ties.method, geodist.args)

  } else if ("list" %in% cls) {

    # List method
    lapply(graph, .mentor_matching, n = n,
           cmode = cmode, lead.ties.method = lead.ties.method,
           geodist.args = geodist.args)

  } else if ("array" %in% cls) {

    # Array method
    apply(graph, 3, .mentor_matching, n = n,
          cmode = cmode, lead.ties.method = lead.ties.method,
          geodist.args = geodist.args)

  } else if ("diffnet" %in% cls) {

    # diffnet method
    g <- graph$graph
    g <- lapply(g, `dimnames<-`, value = list(nodes(graph), nodes(graph)))
    lapply(g, .mentor_matching, n = n,
           cmode = cmode, lead.ties.method = lead.ties.method,
           geodist.args = geodist.args)

  } else stopifnot_graph(graph)

}


.mentor_matching <- function(
  graph,
  n,
  cmode            = "indegree",
  lead.ties.method = "average",
  geodist.args     = list()
) {

  # Step 1. Find the pcent with highest
  d   <- dgr(graph, cmode = cmode)
  r   <- -rank(d, ties.method = lead.ties.method)
  r   <- as.integer(as.factor(r))
  top <- which(r <= n)

  # Step 2: Match each individual with their closest one
  G   <- do.call(
    approx_geodist,
    c(list(graph = as.matrix(graph)), geodist.args)
  )

  ans <- sapply(1:nnodes(graph), function(i) {
    x <- which(G[i,top] == min(G[i,top]))

    # If there are any ties, then solve them by taking a look at i's
    # neighbors
    if (length(x) > 1) {

      # Picking neighbors
      j <- which(graph[i,-top] != 0)

      # If all of them are top
      if (!length(j))
        return(sample(top, size = 1))

      # Average pathlength per top
      j <- sapply(top, function(h) {
        mean(G[j,h])
      })

      # Choose the min
      top[which.min(j)]

    } else top[x]
  })

  # Mentors should be assigned to them selfs
  ans[top] <- top

  # Returning
  structure(
    data.frame(
      name     = nodes(graph),
      degree   = d,
      isleader = 1:nnodes(graph) %in% top,
      match    = nodes(graph)[ans],
      stringsAsFactors = FALSE
    ),
    class    = c("diffnet_mentor", "data.frame"),
    nleaders = length(top),
    graph    = graph
  )

}

#' @export
#' @rdname mentor_matching
leader_matching <- mentor_matching

#' @export
#' @param x An object of class \code{diffnet_mentor}.
#' @param vertex.size Numeric vector of length \code{nnodes(attr(x, "graph"))}.
#' Passed to \code{\link[igraph:plot.igraph]{plot.igraph}}.
#' @param lead.cols Character vector of length \code{attr(x,"nleaders")}. Colors
#' to be applied to each group. (see details)
#' @param vertex.frame.color Passed to \code{\link[igraph:plot.igraph]{plot.igraph}}.
#' @param vertex.label.color Passed to \code{\link[igraph:plot.igraph]{plot.igraph}}.
#' @param edge.arrow.size Passed to \code{\link[igraph:plot.igraph]{plot.igraph}}.
#' @param vshapes Character scalar of length 2. Shapes to identify leaders (mentors)
#' and followers respectively.
#' @param add.legend Logical scalar. When \code{TRUE} generates a legend to distinguish
#' between leaders and followers.
#' @param y Ignored.
#' @param ... Further arguments passed to \code{\link[igraph:plot.igraph]{plot.igraph}}
#' @param main Character scalar. Passed to \code{\link[graphics:title]{title}}
#' @rdname mentor_matching
plot.diffnet_mentor <- function(
  x,
  y            = NULL,
  vertex.size  = rescale_vertex_igraph(dgr(attr(x, "graph")),
                                       minmax.relative.size = c(.01, .03)),
  lead.cols    = grDevices::topo.colors(attr(x, "nleaders")),
  vertex.frame.color = "gray",
  vertex.label.color = "black",
  vshapes      = c(Leader="square", Follower="circle"),
  edge.arrow.size = .5,
  add.legend   = TRUE,
  main         = "Mentoring Network",
  ...) {


  # Creating igraph obj
  ig <- cbind(
    as.character(x[["name"]]),
    as.character(x[["match"]])
  )

  ig <- igraph::graph_from_edgelist(ig)

  # Igraph's indices are not the same as the data!
  indx <- match(igraph::V(ig)$name, x[["name"]])
  ig <- igraph::permute(ig, indx)

  # Creating plot
  graphics::plot.new()
  graphics::plot.window(xlim = c(-1,1), ylim = c(-1,1))

  igraph::V(ig)$shape <- vshapes[2-x[["isleader"]]]

  igraph::plot.igraph(
    x   = ig,
    add = TRUE,
    vertex.color = lead.cols[as.integer(factor(x[["match"]]))],
    vertex.size = vertex.size,
    vertex.frame.color = vertex.frame.color,
    vertex.label.color = vertex.label.color,
    edge.arrow.size = edge.arrow.size,
    ...
  )


  if (add.legend) {

    # Checking names
    if (!length(names(vshapes)))
      names(vshapes) <- c("Leader", "Follower")

    A<-B<-1
    ig <- igraph::make_graph(~A,B)
    igraph::V(ig)$name  <- names(vshapes)
    igraph::V(ig)$color <- "gray"
    igraph::V(ig)$label.color <- "black"
    igraph::V(ig)$shape <- vshapes
    igraph::V(ig)$size  <- 2
    plot(ig, layout = rbind(c(-.4, -1.3), c(.15, -1.3)), add=TRUE,
         vertex.size = 5, rescale=FALSE, vertex.label.dist = 1, vertex.label.degree=0)

  }

  title(main=main)

  # Returning
  invisible(x)

}

# rm(list =ls())
# library(netdiffuseR)
# library(igraph)
#
# set.seed(8821)
# g <- netdiffuseR::rgraph_ba(m = 3, t=19)
# ans <- opinion_leaders(g, n = 3);ans
#
# ig <- graph_from_adjacency_matrix(g)
#
#
# dgr(g, cmode = "indegree", self = TRUE)
# degree(ig, mode="in")
#
# plot(ig, vertex.color = ans$isleader)



