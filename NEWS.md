# netdiffuseR 1.16.2.XXX (beta)

* Attributes in diffnet objects are now stored as data frames (instead of
  matrices). This affects the function `diffnet.attrs`,  
  `diffnet.attrs<-` and `egonet_attrs` as
  these use attributes directly.
  
* In `edgelist_to_adjmat` `use.incomplete` has been replaced by `keep.isolated`
  which makes more sense for naming. Incomplete cases on `times` or `weights` are
  still ignored (as these cannot be processed by the c++ 'engine').
  
* New function: `survey_to_diffnet`. This function allows importing network
  nomination data (in survey fashion) of both types, cross-section and panel
  formats (static network only varying adoption, or dynamic network varying
  attributes and network structure simultaneously).

* New method: `as.array.diffnet`.

* Fixed bug in `diffnet.attrs(..., as.df=TRUE)`. ids were wrongly retrieved. 

* In `edgelist_to_adjmat` `times` has been replaced by `t0` and `t1`. So now
  the user can import graphs with spells.
  
* New functions: `read_pajek` and `read_ucinet`. Still on development.

# netdiffuseR 1.16.2 (CRAN release)

* First CRAN version.



