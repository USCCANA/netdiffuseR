context("Structural equivalence")

# ------------------------------------------------------------------------------
test_that("Computation", {
  # Basis graph
  rn <- LETTERS[1:4]
  graph <- matrix(c(0,0,1,0,1,0,1,0,0,1,0,1,0,0,0,0), ncol=4)

  dyngraph <- lapply(1:3,function(x) methods::as(graph, "dgCMatrix"))
  # Static graphs
  x <- list(
    `static matrix`=struct_equiv(graph),
    `static dgCMatrix`= struct_equiv(methods::as(graph, "dgCMatrix")),
    `dynamic array` = struct_equiv(array(graph, dim = c(4,4,3))),
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
  dimnames(se) <- list(1:4,1:4)
  dimnames(d) <- list(1:4,1:4)

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
})

# ------------------------------------------------------------------------------
test_that("Printing", {
  set.seed(1122)
  dn <- rdiffnet(100,5)

  ans <- struct_equiv(dn)

  expect_output(print(ans), "nodes : 100")
  expect_output(print(ans), "slices: 5")
})

# ------------------------------------------------------------------------------
test_that("By group", {
  # METHOD 1: Using the c.diffnet method:

  # Creating subsets by city
  cities <- unique(medInnovationsDiffNet[["city"]])

  diffnet <- medInnovationsDiffNet[medInnovationsDiffNet[["city"]] == cities[1]]
  diffnet[["expo_se"]] <- exposure(diffnet, alt.graph="se", valued=TRUE)

  for (v in cities[-1]) {
    diffnet_v <- medInnovationsDiffNet[medInnovationsDiffNet[["city"]] == v]
    diffnet_v[["expo_se"]] <- exposure(diffnet_v, alt.graph="se", valued=TRUE)
    diffnet <- c(diffnet, diffnet_v)
  }

  # We can set the original order (just in case) of the data
  diffnet <- diffnet[medInnovationsDiffNet$meta$ids]

  # Checking everything is equal
  expect_equal(
    summary(medInnovationsDiffNet, no.print=TRUE),
    summary(diffnet, no.print=TRUE)
    )


  # METHOD 2: Using the 'groupvar' argument
  # Further, we can compare this with using the groupvar
  diffnet[["expo_se2"]] <- exposure(diffnet, alt.graph="se",
                                    groupvar="city", valued=TRUE)

  # These should be equivalent
  expect_equivalent(
    diffnet[["expo_se", as.df=TRUE]],
    diffnet[["expo_se2", as.df=TRUE]]
    )

  # METHOD 3: Computing exposure, rbind and then adding it to the diffnet object
  expo_se3 <- NULL
  for (v in unique(cities))
    expo_se3 <- rbind(
      expo_se3,
      exposure(
        diffnet[diffnet[["city"]] == v],
        alt.graph = "se", valued=TRUE
      ))

  # Just to make sure, we sort the rows
  expo_se3 <- expo_se3[diffnet$meta$ids,]

  diffnet[["expo_se3"]] <- expo_se3

  # These should be equivalent
  expect_equivalent(
    diffnet[["expo_se", as.df=TRUE]],
    diffnet[["expo_se3", as.df=TRUE]]
  )

  # METHOD 4: Using the groupvar in struct_equiv
  se <- struct_equiv(diffnet, groupvar="city")
  se <- lapply(se, "[[", "SE")
  se <- lapply(se, function(x) {
    x <- 1/x
    x[!is.finite(x)] <- 0
    x
  })

  diffnet[["expo_se4"]] <- exposure(diffnet, alt.graph=se, valued=TRUE)

  # These should be equivalent
  expect_equivalent(
    diffnet[["expo_se", as.df=TRUE]],
    diffnet[["expo_se4", as.df=TRUE]]
  )

})

# ------------------------------------------------------------------------------
test_that("transformGraphBy", {
  set.seed(123)
  x <- rdiffnet(50, 5)
  x[["group"]] <- sample(1:3, nnodes(x), TRUE)

  # Baseline computation
  myfun <- function(x) struct_equiv(x)$SE
  ans0  <- Map(function(G) {
    # Empty matrix
    ans <- methods::new("dgCMatrix", Dim=c(nnodes(x),nnodes(x)), p=rep(0L,nnodes(x) + 1L))

    # Analizing data
    Group <- x[["group"]]
    GROUP <- unique(Group)

    for (i in GROUP) {
      index <- which(Group == i)
      ans[index,index] <- myfun(G[index,index])
    }

    ans

  }, G=x$graph)

  # Method
  ans1 <- Map(function(G) transformGraphBy(G, x[["group"]], myfun), G=x$graph)
  ans2 <- transformGraphBy(x, x[["group"]], myfun)$graph

  expect_equivalent(ans1, ans0)
  expect_equivalent(ans1, ans2)

})
