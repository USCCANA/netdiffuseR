# netdiffuseR 1.16.3.X (beta)

* Attributes in diffnet objects are now stored as data frames (instead of
  matrices). This affects the function `diffnet.attrs`,  
  `diffnet.attrs<-` and `egonet_attrs` as
  these use attributes directly. (Requested by Tom)
  
* In `edgelist_to_adjmat` `use.incomplete` has been replaced by `keep.isolated`
  which makes more sense for naming. Incomplete cases on `times` or `weights` are
  still ignored (as these cannot be processed by the c++ 'engine'). (Reported by Tom)
  
* Bug fixed on `edgelist_to_adjmat`: Counting number of vertices is now done
  right after `recode`. (Reported by Tom)
  
* New function: `survey_to_diffnet`. This function allows importing network
  nomination data (in survey fashion) of both types, cross-section and panel
  formats (static network only varying adoption, or dynamic network varying
  attributes and network structure simultaneously).

* New method: `as.array.diffnet`.

* Fixed bug in `diffnet.attrs(..., as.df=TRUE)`. ids were wrongly retrieved. 

* In `edgelist_to_adjmat` `times` has been replaced by `t0` and `t1`. So now
  the user can import graphs with spells.
  
* New functions: `read_pajek` and `read_ucinet`. Still on development.

* New functions: `nvertices` and `nedges` return the number of vertices and
  edges that a graph has. This can be applied to any class of graph accepted
  by the package.
  
* Added new elements to the `diffnet_struct_test` class: `p.value`, `t0`,
  `mean_t`, and `R`. All these were available before either to be computed
  or retrieved from the `boot` list at the class.
  
* New argument for `threshold`. Now, by default, threshold levels are not computed
  for adopters in the first time period as this can be a biased estimate. If
  the user wants to compute such, he/she can set `include_censored=TRUE`.
  
* New indexing methods via `[[.diffnet`, `[[<-.diffnet`, for network attributes
  and `[.diffnet` and `[<-.diffnet` for adjacency matrix.

# netdiffuseR 1.16.2 (CRAN release)

* First CRAN version.



