#' @importFrom Rcpp evalCpp
NULL

#' @importFrom sna as.sociomatrix.sna
#' @importFrom igraph graph_from_adjacency_matrix set_vertex_attr
#'  any_multiple graph_attr_names as_adj is.loop set_graph_attr V permute make_graph
#'  layout_nicely graph_attr list.edge.attributes graph_from_adj_list
#' @importFrom network as.edgelist is.multiplex is.directed has.loops as.network
#'  get.network.attribute list.vertex.attributes
#' @importFrom networkDynamic networkDynamic network.extract network.collapse
NULL

#' @useDynLib netdiffuseR, .registration = TRUE
NULL

# Importing from the Matrix pkg ------------------------------------------------

#' @import methods
#' @importClassesFrom Matrix dgCMatrix
#' @importFrom Matrix Matrix rowSums colSums
#' @importMethodsFrom Matrix t
NULL

#' @importClassesFrom SparseM matrix.csc
NULL

# Importing from R CORE packages -----------------------------------------------

#' @importFrom grDevices grey rgb colorRampPalette blues9
#' @importFrom graphics grid par plot points symbols text layout legend lines
#'  matplot plot.new plot.window hist mtext polygon image title .filled.contour
#'  box rect
#' @importFrom stats complete.cases runif reshape setNames ftable sd pnorm var
#' as.formula optim nls coef
#' @importFrom utils getFromNamespace head str
#' @importFrom boot boot
#' @importFrom MASS bandwidth.nrd kde2d
#' @importFrom MatchIt matchit
NULL

release_questions <- function() {
  c(
    "Have you updated the inst/NEWS file?",
    "Have you changed the version+dates in DESCRIPTION and NEWS.md?",
    "Have you added all new files to GIT?",
    "Have you clean the vignettes file (source)?"
    )
}
