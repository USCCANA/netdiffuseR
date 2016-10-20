#' @param graph <%= ifelse(exists("dynamic") && dynamic, "A dynamic graph", "Any class of accepted graph format") %> (see \code{\link{netdiffuseR-graphs}}).
#' <%=ifelse(exists("self") && self, "@param self Logical scalar. When \\code{TRUE} allows loops (self edges).", "") %>
#' <%=ifelse(exists("multiple") && multiple, "@param multiple Logical scalar. When \\code{TRUE} allows multiple edges.", "") %>
#' <%=ifelse(exists("valued") && valued, "@param valued Logical scalar. When \\code{TRUE} weights will be considered. Otherwise non-zero values will be replaced by ones.", "") %>
#' <%=ifelse(exists("undirected") && undirected, "@param undirected Logical scalar. When \\code{TRUE} only the lower triangle will be processed.", "") %>
#' <%=ifelse(exists("toa") && toa, "@param toa Integer vector of length \\eqn{n} with the times of adoption.", "") %>
#' <%=ifelse(exists("slice") && slice, "@param slice Integer scalar. Number of slice to use as baseline for drawing the graph.", "") %>
