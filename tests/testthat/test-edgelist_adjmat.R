################################################################################
# edgelist_to_adjmat <--> adjmat_to_edgelist
################################################################################
context("Read/Write graphs edgelist to adjmat and vice-versa")

# Directed graph with attributes
edgelist <- cbind(
  c(2:4,1),
  1:4
)
w <- rowMeans(edgelist)
set.seed(123)
tim <- sample(1:4, 4, TRUE)

EL_digraph <- list(`edgelist matrix` = edgelist, `edgelist data.frame` = as.data.frame(edgelist))

# Undirected graph, generating the data
edgelist <- cbind(
  c(2,1,4,3),
  c(1,2,3,4)
)
EL_graph <- list(`edgelist matrix` = edgelist, `edgelist data.frame` = as.data.frame(edgelist))

#-Loop over classes
for (g in names(EL_digraph)) {

  # Arguments length
  test_that(paste0("Length of inputs in edgelist_adjmat should match (detecting error) - ",g), {
    expect_error(edgelist_to_adjmat(EL_digraph[[g]], weights = w[-1], undirected=FALSE), "Error.+same length")
    expect_error(edgelist_to_adjmat(EL_digraph[[g]], times = tim[-1], undirected=FALSE), "Error.+same length")
  })

  # Back and forth
  edgelist_recovered <- adjmat_to_edgelist(
    edgelist_to_adjmat(EL_digraph[[g]], undirected = FALSE), undirected = FALSE)
  EL_digraph[[g]] <- EL_digraph[[g]][order(EL_digraph[[g]][,1]),]
  test_that(paste0("Undirected static edgelist-adjmat-edgelist (should hold) - ",g), {
    expect_equivalent(edgelist_recovered, as.matrix(EL_digraph[[g]]))
  })


  # ----------------------------------------------------------------------------
  # Dynamic graphs (explicitly): As lists
  # ----------------------------------------------------------------------------

  # Creating the output graph
  edgelist_recovered <- adjmat_to_edgelist(
    edgelist_to_adjmat(EL_digraph[[g]], times = tim, undirected = FALSE),
    undirected = FALSE)

  ord <- order(edgelist_recovered$edgelist[,1], edgelist_recovered$edgelist[,2])
  edgelist_recovered$edgelist <-edgelist_recovered$edgelist[ord,]
  edgelist_recovered$times <-edgelist_recovered$times[ord]

  ord <- order(EL_digraph[[g]][,1], EL_digraph[[g]][,2])
  EL_digraph[[g]] <- EL_digraph[[g]][ord,]
  tim <- tim[ord]
  # Running the test
  test_that(paste0("Undirected dynamic edgelist-adjmat-edgelist (should hold) - ",g), {
    expect_equivalent(edgelist_recovered$edgelist, as.matrix(EL_digraph[[g]]))
  })

  # ----------------------------------------------------------------------------
  # Dynamic graphs (explicitly): As arrays
  # ----------------------------------------------------------------------------
  array_recovered <-
    lapply(edgelist_to_adjmat(EL_digraph[[g]], times = tim, undirected = FALSE),
           as.matrix)
  array_recovered <- array(unlist(array_recovered),
                              dim=c(dim(array_recovered[[1]]), length(array_recovered))
                              )



}

################################################################################
# edgelist_to_adjmat, adjmat_to_edgelist
################################################################################
