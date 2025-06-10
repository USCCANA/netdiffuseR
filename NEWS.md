# Changes in netdiffuseR version 1.23.0 (2025-06-10)

* New methods for simulating multi-diffusion models.

* It is now possible to simulate diffusions with general disadoption functions.


# Changes in netdiffuseR version 1.22.7 (2024-09-18)

* Minor changes to testing (skip warnings).


# Changes in netdiffuseR version 1.22.5 (2022-11-30)

* Adressing [`roxygen2` issue #1491](https://github.com/r-lib/roxygen2/issues/1491)


# Changes in netdiffuseR version 1.22.5 (2022-11-30)

* Solved warning and errors reported by CRAN before the package was archived.

* New S3 generic functions `is_self`, `is_multiple`, `is_valued`, and
  `is_undirected` allow querying graph information for some methods.
  
* Fixed bug in `diag_expand`. Graphs with self ties were not transformed correctly
  (diagonals were excluded.)


# Changes in netdiffuseR version 1.22.4 (2022-09-16)

* Replaced `getMethod("t"...)` by `t` responding to changes in the
  `Matrix` package.
  

# Changes in netdiffuseR version 1.22.1 (2021-05-27)

* netdiffuseR has now a logo!

* Making updates after changes in knitr and Matrix.


# Changes in netdiffuseR version 1.22.0 (2020-05-17)

* Fixing a new issue regarding structural equivalence calculation. In the
  new version, the function has been fully ported to R, which should avoid
  problems related to the C++ code.

* As documented, `struct_equiv` now returns he euclidean distance matrix (it was
  not doing that).


# Changes in netdiffuseR version 1.21.0 (2020-02-10)

* Getting netdiffuseR back to CRAN. We have fixed issues that arise from the
  next big CRAN release, in particular, matrices are now arrays. That was
  creating issues throughout the package (now fixed).


# Changes in netdiffuseR version 1.20.2 (2019-03-25)

* Repaired broken link in bootnet.


# Changes in netdiffuseR version 1.20.1 (2019-03-22)

* This version has no user level visible changes.

## Other changes

* Changing `PI` macros in C++ code as requested by RcppCore.

* Setting 3.5 seed version for tests so that test won't break following message
  from CRAN.
  



# Changes in netdiffuseR version 1.20.0 (2018-06-06)

## New functions and features

* `diffreg` provides a wrapper of `glm` to run lagged regression models.

* Default colors for `plot_diffnet2`, `plot_infectsucept`, and others are now
  selected from the `viridis` R package, which provides perceptually uniform
  and colorblind proof colors.


## Bug fixes

* `plot_diffnet2` now has the correct scaling in nodes positions.

* `rdiffnet_multiple` calls `library(netdiffuseR)` when using multicore.


## Other changes

* `arrow.width` in `plot_threshold.default` now set to be equal to
  `nslices(graph)/80`.

* `curved` option passed to `plot_._threshold.default`.

* The c++ function `edges_arrow` now returns two different elements (the edge,
  and the arrow).



# Changes in netdiffuseR version 1.19.0 (2017-10-16)

## New functions and features

* `rdiffnet` now allows passing scalars for `threshold.dist`, more over, the user
  can also ask the function to just warn when there is no diffusion instead of
  returning with error.

* `plot.diffnet`, `plot_diffnet`, `plot.diffnet_mentor`, and `plot_diffnet2` use igraph for plotting. 
  Also, users can now pass "degree" for -vertex.size- (vertex.cex not used anymore),
  allowing to "automatically" scale vertices by (in/out/.)degree. Also, 
  plotting arguments like vertex.* or edge.* are standarized so these match igraph.
  
* `plot_diffnet` has a new parameter: `background`.

* `rdiffnet_multiple`, a wrapper of `rdiffnet`, allows performing simulations studies
  by running simulating multiple diffusion networks using `rdiffnet`.
  
* `exposure` has a new parameter: `lags`. By default lags = 0, returns a lagged
  exposure matrix.
  

## Bug fixes

* `igraph_to_diffnet` was failing with the graph had no weights.

* `drop_isolated` was not behaving well for diffnet objects.

* `vertex_covariate_dist` was incorrectly specified. Only the default p=2
  were OK. Now fixed and the tests/ folder includes a tests on this. 
  
* `plot_diffnet2` was not passing `color.ramp` to `drawColorKey`. Now fixed.

* `plot.diffnet_mentor` had a bug. Uncesary permutation of vertices was done,
  but it actually had no visible effect. Similar problem was corrected in
  `diffnet_to_igraph`, and other plot methods using igraph for plotting.


## Other changes

* Replacing some C++ functions by R functions in cases in which there
  was no decrease in performance.
  
* `plot_diffnet` function now has smaller margins, so looks more appealing.

* New examples in vignettes "netdiffuseR showcase: Medical Innovations", and
  "Simulating diffusion networks: Using the `rdiffnet` function".


# Changes in netdiffuseR version 1.18.1 (2017-07-22)

## New functions

* `network_to_diffnet`, `diffnet_to_network` coercion between 
  diffnet and network objects.
  
* `networkDynamic_to_diffnet`, `diffnet_to_networkDynamc` coercion
  between diffnet and networkDynamic objects.


## New features

* `new_diffnet` and `as_diffnet` now receive static networks as well.


## Bug fixes

* `diffnet_to_igraph` was copying over a single adjacency matrix,
  which did a difference in dynamic networks.
  
* `diffnet_to_igraph` was not considering loops correctly.

* In unexpected situations `egonet_attrs` was crashing.


# Changes in netdiffuseR version 1.18.0 (2017-07-16)

## New functions

* `bootnet` implements network bootstrapping based on Snijders and Borgatti (1999)

* `mentor_matching` implements Valente and Davis (1999) Mentor matching algorithm.
  including a plot method.

* `approx_geodesic` an alternative to `igraph::distances` and `sna::geodist`.
  computes geodesics up to a certain number of steps and returns a sparse
  matrix.
  
* `matrix_compare` Efficiently compares two sparse matrices looking only at the valued
  cells.
  
* `as_dgCMatrix` Coerce matrix-like objects into dgCMatrix objects (sparse matrices
  from the `Matrix` package).
  
* `fitbass` Fits the Bass Diffusion Model to an observed vector of cumulative
  adopters. The estimation is done via `stats::nls`.
  
* `netmatch` and `netmatch_prepare` (on development) implement matching estimators
  with network data.

  
## Bugs Fixes

* `dgr` returned with error when `self == TRUE`

* In some calls to `igraph::graph_from_adj...` sorting of vertices was not preserved.

* The `matrix` method in `egonet_attrs` was passing a list of vertices instead of
  the attributes. Fixed.
  
* `transformGraphBy` was returning with error when the time periods ranged other than
  1, 2, ...


## New Features

* `rgraph_er` is now significantly faster (orders of magnitude compare to
  previous versions). `rgraph_ba` is faster too.
  
* `moran` now returns the sd, expected and p-value.

* `exposure` now receives static graphs in `alt.graph` with a warning.

* `rewire_graph` now also uses QAP. This affects directly to `structu_test`.
  

# Changes in netdiffuseR version 1.17.0 (2016-11-10)

## New features and changes

* The title of the package is now _Analysis of Diffusion and Contagion Processes
  on Networks_.

* The function `struct_test` now allows other types of graphs. Before it only
  supported `diffnet` objects.

* The function `rewire_graph` gains a new argument for the algorithm "swap". Now
  to ensure aperiodicity in MCMC a chance of skiping a rewire has been included.

* The function `n_rewires` now has a default of 20 (before it was 100).
  This is based on Ray et al (2012) (more details in the manual).
  
* The function `rgraph_ba` gains a new argument, `self=TRUE`. By default behaves
  as before following Bollobas, but now can deviate to generate graphs with no
  autolinks.
  
* In `rgraph_ba`, the argument `eta` allows implementing De Almeida et al. (2013)
  Scale-free homophilic networks.
  
* The functions `exposure` and `dgr` are now pure R code (C++ functions were
  replace since there were no significant speed gains).
  
* `diffnet` class objects now have two new meta-values: name and behavior.

* Elements -graph-, -toa-, -adopt- and -cumadopt- in `diffnet` class objects
  have lost their dimnames (more efficient storage).
  
* `classify_adopters` now always includes Non-Adopters.
  
## New functions

* `vertex_covariate_dist` computes distances between vertices using both the graph
  and a matrix of length nxK.
  
* `vertex_mahalanobis_dist` computes mahalanobis distance between vertices (as above).

* `struct_test_asymp` an asymptotic approximation of `struct_test` (not recomended).

* `ego_variance` computes a pseudo variance at the ego level (aux function for
  `struct_test`).

* `transformGraphBy` applies a function that transforms a graph considering
  structural zeros given by groups. Similar to the idea of the -by- option
  in `struct_equiv`.
  
* `read_ucinet` read UCINET binary files (both header and graph file). Still
  work in progress.
  
* `plot.diffnet_degSeq` method allows visualizing degree sequence as log-log
  plots (default).
  
* `diag_expand` creates a single adjacency matrix from a dynamic graph.

* `summary.diffnet_adoptChange` method generates a summary table of the df
  generated by `select_egoalter`.
  
* `permute_graph` permutes the values of an adjacency matrix (Conditional 
  Uniform Graph).

* `rewire_qap` generates isomorphic graphs by "changing the labels".

## Bug fixes

* `^.diffnet` method was rasing to +1 power, e.g. `diffnet^2` was actually
  `diffnet^3`.
  
* `/.diffnet` was not working.

* `plot_diffnet` was computing the coordinates of the cells wrongly. Most of the
  time causing adding figures outside of the plotting area.
  
* The `c.diffnet_struct_test` method was not updating the `p.value`.

* The function `edgelist_to_adjmat` was not processing correctly undirected
  graphs when the edgelist represented a lower triangular matrix.
  
* The function `survey_to_diffnet` had an issue processin dyn graph attrs returning
  errors. Now fixed.
  
* The function `select_egoalter` returned error when `graph` was an array.

* The method `[[<-.diffnet` failed when replacing a dynamic attribute with
  a `NULL` value (e.g. `dn[["my_dyn_att"]] <- NULL`).


# Changes in netdiffuseR version 1.16.7 (2016-07-07)

## Bug fixes

* Fixed bug in `struct_equiv`: When `groupvar` was a list (dynamic attr), the
  function returned error.
  
* Fixed bug in `rewire_graph.array`: Returned error when `algorithm="swap"`

* Fixed bug in `rewire_graph`: The option `copy.first` was not been applied correctly.

* In `hist.diffnet_struct_test`: `...` now passed to `hist.default`.

* Fixed bug on `egonet_attrs`: The matrix method was returning with error.

## New features and changes

* `plot_infectsuscept` includes 2D kernel smoother via `MASS::kde2d`.

* `infection`, `susceptibility` and `threshold` now report `NA` for non-adopters
  or excluded variables.
  
* `egonets_attrs` now has new argument: `self.attrs` allows including ego's attributes
  as part of the outcome so it can be used by the user.
  
* `plot_diffnet` now uses `igraph::plot.igraph` for plotting instead

* `threshold` gains a new argument: `lags` now users can define threshold as
  exposure `lags` time periods prior to the time of adoption. By default is 0
  so its exposure at the time of adoption.
  

## New functions

* New method `c.diffnet_struct_test`: A wrapper of `boot:::c.boot`.

* `diffusionMap` computes the required matrix to be used with
  `image`-like functions mapping a vertex covariate given a graph structure.
  
* `n_rewires`: computes a suggested number of rewires per step in order to attain
  unbiased graph samples.
  
* `diffnetLapply`: Apply throught periods on diffnet objects.

* Several new methods for the class `diffnet`. Now users can apply `str`,
  `dimnames` (so `colnames` and `rownames`), `t`, `&`, `|`, `dim` and `%*%`.
  
* `drawColorKey`: Handy function to draw a color key in the current plot.

* `classify_adopters`: As in Valente (1995), depending on time of adoption, adopters
  are classified as early adopters, early  majority, late majority, and laggards.
  The function introduces a new class and has methods for `ftable` and `plot`.
  
* `rescale_vertex_igraph`: Helper function to fix the size of vertex when calling
  `plot.igraph` so that the size is proportional to the x-axis.


# Changes in netdiffuseR version 1.16.5 (2016-05-02)

## Bug fixes

* Bug fixed on `edgelist_to_adjmat`: Counting number of vertices is now done
  right after `recode`. (Reported by Tom)

* Fixed bug in `diffnet.attrs(..., as.df=TRUE)`. ids were wrongly retrieved.

* Fixed bugs for `rgraph_ba_cpp`: Degree of new vertices was not changing apropiately.
  This only was an issue when `m>1`.
  
* Fixed bugs for the `as_diffnet` method for arrays.
  
* Fixed bugs in `rewire_graph`. Indexing of the jth component (when rewiring) was
  not been made correctly (now it does). Also, when rewiring, the new endpoints
  were truncated to n-1 (now fixed).
  
* Fixed bugs for `as_diffnet`: When a dynamic graph was passed with slices names
  different from the time periods, the slices names were kept. Now these are
  replaced by `meta$pers`.


## New features and changes

* Support for `int64_t` in `RcppArmadillo` now allows for creating/reading
  adjacency matrices with more than 4 billion elements (big graphs).

* In `edgelist_to_adjmat` `use.incomplete` has been replaced by `keep.isolated`
  which makes more sense for naming. Incomplete cases on `times` or `weights` are
  still ignored (as these cannot be processed by the c++ 'engine'). (Reported by Tom)
  
* In `edgelist_to_adjmat` `times` has been replaced by `t0` and `t1`. So now
  the user can import graphs with spells.
  
* Added new elements to the `diffnet_struct_test` class: `p.value`, `t0`,
  `mean_t`, and `R`. All these were available before either to be computed
  or retrieved from the `boot` list at the class.
  
* New argument for `threshold`. Now, by default, threshold levels are not computed
  for adopters in the first time period as this can be a biased estimate. If
  the user wants to compute such, he/she can set `include_censored=TRUE`.
  
* Attributes in diffnet objects are now stored as data frames (instead of
  matrices). This affects the function `diffnet.attrs`, and `egonet_attrs` as
  these use attributes directly. (Requested by Tom)

* New features for the `rewire_graph` function. In particular, `p` can now be
  a vector of length `T`, so each slice can have different rewiring prob., and
  the new option `copy.first` which allows to recycle the first rewired slice
  (see details).
  
* New features for the `exposure` function. When `graph` is of class diffnet, 
  the function accepts `attrs` equal to the name of some the graph's attributes.
  Also, `alt.graph` can be specified as `se`, which will be replaced by the
  inverse of the structural equivalence. When `valued=FALSE` the function will
  switch it to `TRUE` and warn the user.
  
* New argument for `struct_equiv` and `exposure`, `groupvar`. This new option
  provides a convenient way of calculating structural equivalence and
  exposure clusterized by group. Specially useful when there are different
  communities in a graph. See examples in the manual.
  
* `as_diffnet` now has an internal function, `check_as_diffnet_attrs`, to check
  input attributes dimensions and coerce them into proper class/structure. Valid
  attributes are now documented in the function's manual.
  
* New arguments for `edges_coords`: `dev` and `ran` allow including device +
  margins aspec ratio and plotting area y/x limits for improved aspect ratio
  computation.

* New internal function `edges_arrow`: Computes the coordinates of a 4 points
  polygon allowing to draw pretty arrows considering aspec ratio of device,
  margins and y/x.

* Geodesic distances are now computed using `igraph::distances` instead of
  `sna::geodist` as it is more flexible and faster.
  
* New arguments for `plot_threshold`: `vertex.sides`, `vertex.lab.cex`, `vertex.lab.adj`,
  `vertex.lab.col`, `vertex.rot`, `jitter.factor`, and `jitter.amount` to give
  more control.
  
* New internal function `vertex_coords`: Creates polygons of any given number of
  sides considering aspec ratio of both x/y and device.
  
* New features for `rdiffnet`. `seed.graph` can be either a function that generates
  a random graph, a character string (as before) indicating the class of graph
  to generate, or any other class of graph (either static or dynamic) as specified
  in `netdiffuseR-graphs`. `seed.nodes` can now be a vector with indices pointing
  to the initial adopters.
  
* The rewiring algorithm for `rgraph_ws` has been replaced with a `rewire_ws` which
  has been implemented as it was presented in Watts and Strogatz (1998).

## New functions

* New function: `survey_to_diffnet`. This function allows importing network
  nomination data (in survey fashion) of both types, cross-section and panel
  formats (static network only varying adoption, or dynamic network varying
  attributes and network structure simultaneously).
  
* New function: `edgelist_to_diffnet`. Similar to `survey_to_diffnet`, this
  function reads diffusion networks from an adjacency matrix and a vertex
  attributes data frame. Both the attributes and the edgelist can be static
  or dynamic.

* New method: `as.array.diffnet`.

* New functions: `read_pajek` and `read_ucinet`. Still on development.

* New functions: `nvertices` and `nedges` return the number of vertices and
  edges that a graph has. This can be applied to any class of graph accepted
  by the package.

* New indexing methods via `[[.diffnet`, `[[<-.diffnet`, for network attributes
  and `[.diffnet` and `[<-.diffnet` for adjacency matrix. The function
  `diffnet.attrs<-` will be deprecated for the next CRAN release. The function
  `diffnet.subset.slices` is now not exported (internal use), so the user
  needs to use the `[.diffnet` method instead.

* New concatenating method `c.diffnet` for diffnet objects. This method allows
  'adding up' diffnet objects.

* New print method for `diffnet_se`, objects returned by `struct_equiv`.

* New function `diffnet_to_igraph`.

* New rewiring algorithm, `rewire_swap` implements the edge-switch algorithm in
  an efficient way. This preserves degree sequences.


# Changes in netdiffuseR version 1.16.2 (2016-02-18)

* First CRAN version.



