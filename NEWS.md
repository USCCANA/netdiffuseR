# netdiffuseR 1.16.2.XXX (beta)

* Attributes in diffnet objects are now stored as data frames (instead of
  matrices). This affects the function `diffnet.attrs`,  
  `diffnet.attrs<-` and `egonet_attrs` as
  these use attributes directly.
  
* In `edgelist_to_adjmat` `use.incomplete` has been replaced by `keep.isolated`
  which makes more sense for naming. Incomplete cases on `times` or `weights` are
  still ignored (as these cannot be processed by the c++ 'engine').
  
* New function: `survey_to_diffnet`.

* New method: `as.array.diffnet`.

* Fixed bug in `diffnet.attrs(..., as.df=TRUE)`. ids were wrongly retrieved. 

# netdiffuseR 1.16.2 (CRAN release)

* First CRAN version.



