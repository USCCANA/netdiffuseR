# Development road

The following document details expectations for the comming versions. There is no guarantee on the implementation of any of the points here.

# CRAN 1.16.4: Expected

## Polishing

- Review manuals: Description of arguments and name of the functions.
- ~~Aliases creation: `select_egoalter` `table_...`~~

## Testing

- `survey_to_diffnet`, `edgelist_to_adjmat`
- `diffnet_as_igraph`, `igraph_as_diffnet`


## Developing

- `diffnet_to_networkDynamic`, `networkDynamic_to_diffnet`. Need to reach out the author of Carter Butts: Is there any way to access networkDynamic objects formal definition??

- ~~`rewire_dgr_preserve`: A brief comparison on igraph and netdiffuseR rewiring algorithms shows that igraph has no significan speed improvement in small-medium graphs. This may be due to having a similar rewiring algorithm. The dgr preserve should work in a similar fashion and should incorporate the call to `std::<vector>remove` method.~~

- Raise the level of `ARMA_BIT_64INT`: Create a toy package in which, by simulations, we create big sparse matrices of class `dgCMatrix` do read/write operations both from R and from C++.

- `read_ucinet` and rename `read_ucinet` to `read_dl`. And write versions of it. It will be useful to ask Steven B. if is there any way of using UCINET to analyze dynamic data.

- `read_net` Talk with the developer to see is they have more information (aditional reference) on the data formatting for dynamic graphs and attributes. Recall UCINET's documents provided with spreadsheets pointing to data types and structure.

## Improving

-   ~~`rdiffnet`:~~
    Set initial adopters by providing a vector:
    - logical
    - integer (of length n)
    - integer (of random length)
    - character (ids, rownames)

-   What happens in `survey_to_diffnet` when a data.frame is provided, but se set to fill data and the checking done over the dimmensions of the data?

## Examples

- Include the survival analysis example: Maybe use the recidivism data mixed with social data? Find a paper where someone does that
- Include a probit/logit model
- Write down an example with the function `select_egoalter` and rename it.
