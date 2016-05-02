#' @importFrom Rcpp evalCpp
NULL

#' @importFrom sna gplot as.sociomatrix.sna
NULL

#' @importFrom igraph graph_from_adjacency_matrix set_vertex_attr distances
NULL

#' @useDynLib netdiffuseR
NULL

# Importing from the Matrix pkg ------------------------------------------------

#' @import methods
#' @importClassesFrom Matrix dgCMatrix
#' @importFrom Matrix Matrix rowSums colSums
NULL

#' @importClassesFrom SparseM matrix.csc
NULL

# Importing from R CORE packages -----------------------------------------------

#' @importFrom grDevices grey rgb
#' @importFrom graphics grid par plot points symbols text layout legend lines matplot plot.new plot.window hist mtext polygon
#' @importFrom stats complete.cases runif reshape setNames
#' @importFrom utils getFromNamespace head
#' @importFrom boot boot
NULL

release_questions <- function() {
  c(
    "Have you updated the inst/NEWS file?",
    "Have you changed the version+dates in DESCRIPTION and NEWS.md?",
    "Have you added all new files to GIT?",
    "Have you clean the vignettes file (source)?"
    )
}
