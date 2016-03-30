context("Structural equivalence")

# Basis graph
rn <- LETTERS[1:4]
graph <- matrix(c(0,0,1,0,1,0,1,0,0,1,0,1,0,0,0,0), ncol=4,
                dimnames = list(rn, rn))

dyngraph <- lapply(1:3,function(x) methods::as(graph, "dgCMatrix"))
# Static graphs
x <- list(
  `static matrix`=struct_equiv(graph),
  `static dgCMatrix`= struct_equiv(methods::as(graph, "dgCMatrix")),
  `dynamic array` = struct_equiv(array(graph, dim = c(4,4,3), dimnames = list(rn, rn, 1:3))),
  `dynamic list` = struct_equiv(dyngraph),
  `dynamic diffnet` = struct_equiv(as_diffnet(dyngraph, c(1L,1L,3L,2L)))
)

# `Manual` calculations

Z <- sna::geodist(graph, inf.replace = 0)$gdist
Z <- Z/max(Z)

d <- matrix(0, ncol=4, nrow=4)
for (i in 1:4)
  for (j in 1:4) {
    if (i==j) next
    zjizij <- (Z[i,j] - Z[j,i])^2
    sum_zjkzki <- sum((Z[j,c(-i,-j)] - Z[i,c(-i,-j)])^2)
    sum_zkjzik <- sum((Z[c(-i,-j),j] - Z[c(-i,-j),i])^2)
    d[i,j]<-sqrt(zjizij + sum_zjkzki + sum_zkjzik)
  }
dmax <- apply(d, 1, max)

se <- matrix(0, ncol=4, nrow=4)
for (i in 1:4)
  for (j in 1:4) {
    if (i==j) next
    se[i,j] <- (dmax[i] - d[j,i])/sum(dmax[i] - d[-i,i])
  }

# Naming
dimnames(se) <- dimnames(graph)
dimnames(d) <- dimnames(graph)

# Comparing
for (i in names(x)) {
  if (grepl("static", i)) {
    expect_equal(d, x[[i]]$d, tolerance=getOption("diffnet.tol"), scale=1)
    expect_equal(se, x[[i]]$SE, tolerance=getOption("diffnet.tol"), scale=1)
  } else {
    expect_equal(d, x[[i]][[1]]$d, tolerance=getOption("diffnet.tol"), scale=1)
    expect_equal(se, x[[i]][[1]]$SE, tolerance=getOption("diffnet.tol"), scale=1)
  }
}
