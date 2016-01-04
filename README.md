netdiffuseR: Network Analysis for Diffusion of Innovations
==========================================================

[![Build Status](https://travis-ci.org/USCCANA/netdiffuseR.svg?branch=master)](https://travis-ci.org/USCCANA/netdiffuseR) [![Build status](https://ci.appveyor.com/api/projects/status/6u48wgl1lqak2jum?svg=true)](https://ci.appveyor.com/project/gvegayon/netdiffuser) [![codecov.io](https://codecov.io/github/USCCANA/netdiffuseR/coverage.svg?branch=master)](https://codecov.io/github/USCCANA/netdiffuseR?branch=master)

This package contains functions useful for analyzing network data for diffusion of innovations applications.

The package was developed as part of the paper Thomas W. Valente, Stephanie R. Dyal, Kar-Hai Chu, Heather Wipfli, Kayo Fujimoto, *Diffusion of innovations theory applied to global tobacco control treaty ratification*, Social Science & Medicine, Volume 145, November 2015, Pages 89-97, ISSN 0277-9536 (available [here](http://www.sciencedirect.com/science/article/pii/S027795361530143X))

Functions include selection, susceptibility & infection, etc.

Installation
------------

Using the `devtools` package, you can install `netdiffuseR` dev version as follows

``` r
devtools::install_github('USCCANA/netdiffuseR', build_vignettes = TRUE)
```

You can skip building vignettes by setting `build_vignettes = FALSE` (so it is not required).

For the case of OSX users, there seems to be a problem when installing packages depending on `Rcpp`. This issue, developed [here](https://github.com/USCCANA/netdiffuseR/issues/3), can be solved by open the terminal and typing the following

``` sh
curl -O http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2
sudo tar fvxz gfortran-4.8.2-darwin13.tar.bz2 -C /
```

before installing the package through `devtools`.

Examples
--------

This example has been taken from the package's vignettes:

``` r
library(netdiffuseR)
```

### Infectiousness and Susceptibility

``` r
# Generating a random graph
set.seed(1234)
n <- 100
nper <- 20
graph <- rgraph_er(n,nper, p=.40, undirected = FALSE)
toa <- sample(c(1:(1+nper-1), NA), n, TRUE)
head(toa)
```

    ## [1] 20  7  9  1 17 14

``` r
# Creating a diffnet object
diffnet <- as_diffnet(graph, toa)
diffnet
```

    ## Dynamic network of class -diffnet-
    ##  # of nodes        : 100
    ##  # of time periods : 20
    ##  Adoption rate     : 0.95
    ##  Type              : directed

``` r
summary(diffnet)
```

    ## Diffusion network summary statistics
    ## -----------------------------------------------------------------------
    ##  Period  Adopters Cum Adopt. Cum Adopt. % Hazard Rate Density Moran's I 
    ## -------- -------- ---------- ------------ ----------- ------- --------- 
    ##        1        3          3         0.03           -    0.80     -0.01 
    ##        2        8         11         0.11       11.00    0.81     -0.01 
    ##        3        4         15         0.15       15.00    0.81     -0.01 
    ##        4        1         16         0.16       16.00    0.79     -0.01 
    ##        5        2         18         0.18       18.00    0.79     -0.01 
    ##        6        5         23         0.23       23.00    0.81     -0.01 
    ##        7        6         29         0.29       29.00    0.82     -0.02 
    ##        8        1         30         0.30       30.00    0.80     -0.01 
    ##        9        6         36         0.36       36.00    0.81     -0.01 
    ##       10        5         41         0.41       41.00    0.78     -0.01 
    ##       11        2         43         0.43       43.00    0.80     -0.01 
    ##       12        6         49         0.49       49.00    0.80     -0.01 
    ##       13        4         53         0.53       53.00    0.79     -0.01 
    ##       14        7         60         0.60       60.00    0.80     -0.01 
    ##       15        7         67         0.67       67.00    0.80     -0.01 
    ##       16        3         70         0.70       70.00    0.80     -0.01 
    ##       17       10         80         0.80       80.00    0.79     -0.01 
    ##       18        4         84         0.84       84.00    0.79     -0.01 
    ##       19        2         86         0.86       86.00    0.81     -0.01 
    ##       20        9         95         0.95       95.00    0.81     -0.01 
    ## -----------------------------------------------------------------------
    ##  Left censoring  : 0.03 (3)
    ##  Right centoring : 0.05 (5)

``` r
# Visualizing distribution of suscep/infect

out <- plot_infectsuscep(diffnet, K=1, logscale = TRUE)
```

    ## Warning in plot_infectsuscep.list(graph$graph, graph$toa, normalize, K, :
    ## When applying logscale some observations are missing.

![](README_files/figure-markdown_github/plot_infectsuscept-1.png)

### Threshold

``` r
# Generating a random graph
set.seed(123)
n <- 6
nper <- 5
toa  <- sample(c(2000:(2000+nper-1), NA), n, TRUE)
nper <- length(unique(toa))
graph <- rgraph_er(n,nper, p=.3, undirected = FALSE)

# Creating a diffnet object
diffnet <- as_diffnet(graph, toa, undirected = FALSE)

# Computing exposure
expos <- exposure(diffnet)

# Threshold with fixed vertex size
plot_threshold(diffnet)
```

![](README_files/figure-markdown_github/plot_threshold-1.png)

``` r
# Threshold with vertex size = avg degree
cex <- rowMeans(dgr(graph))
cex <- (cex - min(cex) + 1)/(max(cex) - min(cex) + 1)/2
plot_threshold(diffnet, vertex.cex = cex)
```

![](README_files/figure-markdown_github/plot_threshold-2.png)

### Diffusion process

``` r
plot_diffnet(diffnet)
```

    ## Loading required package: SparseM
    ## 
    ## Attaching package: 'SparseM'
    ## 
    ## The following object is masked from 'package:base':
    ## 
    ##     backsolve

![](README_files/figure-markdown_github/plot_diffnet-1.png)

To-do list
----------

-   Import/Export functions for interfacing other package's clases, in particular: `statnet` set (specially the packages `networkDynamic` and `ndtv`), `igraph` and `Rsiena`.
-   Populate the tests folder.
-   Use spells? (`select_egoalter` would use this)
