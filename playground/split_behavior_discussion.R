farm <- brfarmersDiffNet
str(farm)
farm$vertex.static.attrs$liveout
class(farm$vertex.dyn.attrs)
class(farm$vertex.static.attrs)

med_data <- medInnovationsDiffNet
str(med_data)
class(med_data$vertex.static.attrs$city)

k_fam <- kfamilyDiffNet
str(k_fam)
class(k_fam$vertex.static.attrs$study)

n=40; t=5
diffnet <- rdiffnet(40, 5, seed.p.adopt = .2)
X <- matrix(diffnet[["real_threshold"]], ncol=t, nrow=n, byrow = FALSE)
#ans0 <- exposure(diffnet, attrs=X)
net_1 <- rdiffnet(n, t, seed.nodes = 'random',
                  exposure.args = list(attrs = matrix(runif(n), nrow=n, ncol=t, byrow = FALSE)))
summary(net_1)
str(net_1)
class(net_1$toa)
class(net_1$adopt)
class(net_1$cumadopt)
class(net_1$vertex.static.attrs)
class(net_1$graph.attrs)
class(net_1$meta)

n=40; t=5
net_2 <- rdiffnet(n, t, seed.p.adopt = list(0.5,0.5),
                  exposure.args = list(attrs = X))
summary(net_2)
str(net_2)
class(net_2$toa)
class(net_2$adopt)
class(net_2$cumadopt)
class(net_2$vertex.static.attrs)
class(net_2$graph.attrs)
class(net_2$meta)

split_behaviors <- function(diffnet_obj) {

  diffnets <- rep(diffnet_obj, ncol(diffnet_obj$toa))
  diffnets_list <- list()

  #ver_static_att_nams <- colnames(diffnet_obj$vertex.static.attrs)

  for (q in 1:ncol(diffnet_obj$toa)) {

    for (i in seq_along(diffnet_obj)) {
      if (!is.null(diffnets[i]$toa)) {
        #print(diffnets[i]$toa)
        diffnets[i]$toa <- diffnet_obj$toa[, q, drop = FALSE]
      } else if (!is.null(diffnets[i]$adopt)) {
        diffnets[i]$adopt <- diffnet_obj$adopt[[q]]
      } else if (!is.null(diffnets[i]$cumadopt)) {
        diffnets[i]$cumadopt <- diffnet_obj$cumadopt[[q]]
      }# else if (!is.null(diffnets[i]$vertex.dyn.attrs)) {
      #  diffnets[i]$vertex.dyn.attrs <- setNames(data.frame(diffnet_obj$vertex.static.attrs[, q]), ver_static_att_nams[q])
      #}
    }

    # diffnets[2]$toa <- diffnet_obj$toa[, q, drop = FALSE]
    # diffnets[[q]]$adopt <- diffnet_obj$adopt[[q]]
    # diffnets[[q]]$cumadopt <- diffnet_obj$cumadopt[[q]]

    diffnets_list[[q]] <- diffnets[q*(1:length(diffnet_obj))]
  }

  return(diffnets_list)

  # for (q in ncol(net_2$toa)) {
  #
  #   #graph <- net_2$graph
  #   #
  #   diff_obj_2$toa <- net_2$toa[,q]
  #   #class(toa_slice)
  #   adopt_slice <- net_2$adopt[[q]]
  #   #class(adopt_slice)
  #   cumadopt_slice <- net_2$cumadopt[[q]]
  #   #class(cumadopt_slice)
  #   ver_static_att_slice <- setNames(data.frame(net_2$vertex.static.attrs[, q]), ver_static_att_nams[q])
  #   #class(ver_static_att_slice)
  #
  #   meta_slice$behavior <- strsplit(meta_slice$behavior, ", ")[[1]][q]
  #   #class(meta_slice)
  # }
}

###############################################################################

split_behaviors <- function(diffnet_obj) {

  diffnets <- replicate(ncol(diffnet_obj$toa), diffnet_obj, simplify = FALSE)

  for (q in 1:ncol(diffnet_obj$toa)) {

    diffnets[[q]]$toa <- as.integer(diffnet_obj$toa[, q, drop = FALSE])
    names(diffnets[[q]]$toa) <- rownames(diffnet_obj$toa)

    diffnets[[q]]$adopt <- diffnet_obj$adopt[[q]]

    diffnets[[q]]$cumadopt <- diffnet_obj$cumadopt[[q]]

  }
  return(diffnets)
}

test_that("toa, adopt, and cumadopt should be equal! (split_behaviors)", {
  set.seed(12131)
  n            <- 50
  t            <- 5
  graph        <- rgraph_ws(n, 4, p=.3)
  seed.nodes   <- c(1,5,7,10)
  thr          <- runif(n, .2,.4)

  # Generating identical networks
  net_single <- rdiffnet(seed.graph = graph, seed.nodes = seed.nodes, seed.p.adopt = 0.1,
                         t = t, rewire = FALSE, threshold.dist = thr)

  net_multiple <- rdiffnet(seed.graph = graph, seed.nodes = seed.nodes, seed.p.adopt = list(0.1, 0.1),
                           t = t, rewire = FALSE, threshold.dist = thr)

  net_single_from_multiple <- split_behaviors(net_multiple)
  net_single_from_multiple_1 <- net_single_from_multiple[[1]]

  expect_equal(net_single_from_multiple_1$toa, net_single$toa) # Error: names for current but not for target
  expect_equal(net_single_from_multiple_1$adopt, net_single$adopt)
  expect_equal(net_single_from_multiple_1$cumadopt, net_single$cumadopt)
})

# Let's check the plots.

plot_diffnet(net_single$graph, net_single$adopt)
plot_diffnet(net_single_from_multiple_1$graph, net_single_from_multiple_1$adopt)

plot_infectsuscep(net_single$graph, net_single$toa)
plot_infectsuscep(net_single_from_multiple_1$graph, net_single_from_multiple_1$toa)

set.seed(1234) # they are almost the same
plot_threshold(net_single$graph,
               exposure(net_single$graph, net_single$cumadopt),
               net_single$toa)
plot_threshold(net_single_from_multiple_1$graph,
               exposure(net_single_from_multiple_1$graph, net_single$cumadopt),
               net_single_from_multiple_1$toa)

################################################################################

set.seed(1234)
net_2 <- rdiffnet(50,5, seed.p.adopt = list(0.1, 0.1))
#str(net_2)
net_2_splitted <- split_behaviors(net_2)
net_1_from_2 <- net_2_splitted[[1]]
expect_s3_class(net_1_from_2, "diffnet")
#str(net_1_from_2)
#str_net_1_from_2 <- capture.output(str(net_1_from_2))

net_1 <- rdiffnet(50,5, seed.p.adopt = 0.1, seed.nodes = c(1,2,3,4,5))
net_1_1 <- rdiffnet(50,5, seed.p.adopt = 0.1, seed.nodes = c(1,2,3,4,5))
expect_equivalent(net_1$toa, net_1_1$toa)

#str(net_1)
#str_net_1 <- capture.output(str(net_1))

expect_equivalent(net_1_from_2$toa, net_1$toa)
expect_equal(net_1_from_2$adopt, net_1$adopt)
expect_equal(net_1_from_2$cumadopt, net_1$cumadopt)
#identical(str_net_1_from_2, str_net_1)
