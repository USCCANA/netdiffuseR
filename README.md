
[![R-CMD-check](https://github.com/USCCANA/netdiffuseR/actions/workflows/r.yml/badge.svg)](https://github.com/USCCANA/netdiffuseR/actions/workflows/r.yml)
[![Build
status](https://ci.appveyor.com/api/projects/status/6u48wgl1lqak2jum?svg=true)](https://ci.appveyor.com/project/gvegayon/netdiffuser)
[![codecov.io](https://codecov.io/github/USCCANA/netdiffuseR/coverage.svg?branch=master)](https://codecov.io/github/USCCANA/netdiffuseR?branch=master)
[![](https://cranlogs.r-pkg.org/badges/netdiffuseR)](https://cran.r-project.org/package=netdiffuseR)
[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/netdiffuseR)](https://cran.r-project.org/package=netdiffuseR)
[![](https://cranlogs.r-pkg.org/badges/grand-total/netdiffuseR)](https://cran.r-project.org/package=netdiffuseR)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1039317.svg)](https://doi.org/10.5281/zenodo.1039317)
[![Dependencies](https://tinyverse.netlify.com/badge/netdiffuseR)](https://cran.r-project.org/package=netdiffuseR)
[![USC’s Department of Preventive
Medicine](https://raw.githubusercontent.com/USCbiostats/badges/master/tommy-uscprevmed-badge.svg)](https://preventivemedicine.usc.edu)

<h1>

netdiffuseR: Analysis of Diffusion and Contagion Processes on
Networks<img src="man/figures/netdiffuser-logo.svg" style="max-width:200px;width:50%;" align="right">

</h1>

This package contains functions useful for analyzing network data for
diffusion of innovations applications.

The package was developed as part of the paper Thomas W. Valente,
Stephanie R. Dyal, Kar-Hai Chu, Heather Wipfli, Kayo Fujimoto,
*Diffusion of innovations theory applied to global tobacco control
treaty ratification*, Social Science & Medicine, Volume 145, November
2015, Pages 89-97, ISSN 0277-9536 (available
[here](https://www.sciencedirect.com/science/article/pii/S027795361530143X))

From the description:

> Empirical statistical analysis, visualization and simulation of
> diffusion and contagion processes on networks. The package implements
> algorithms for calculating network diffusion statistics such as
> transmission rate, hazard rates, exposure models, network threshold
> levels, infectiousness (contagion), and susceptibility. The package is
> inspired by work published in Valente, et al., (2015); Valente (1995),
> Myers (2000), Iyengar and others (2011), Burt (1987); among others.

**Acknowledgements**: netdiffuseR was created with the support of grant
R01 CA157577 from the National Cancer Institute/National Institutes of
Health.

``` r
citation(package="netdiffuseR")
```

``` 

To cite netdiffuseR in publications use:

  Valente TW, Vega Yon GG. Diffusion/Contagion Processes on Social
  Networks. Health Education & Behavior. 2020;47(2):235-248.
  doi:10.1177/1090198120901497

Vega Yon G, Valente T (2021). _netdiffuseR: Analysis of Diffusion and
Contagion Processes on Networks_. doi: 10.5281/zenodo.1039317 (URL:
https://doi.org/10.5281/zenodo.1039317), R package version 1.22.2,
<URL: https://github.com/USCCANA/netdiffuseR>.

To see these entries in BibTeX format, use 'print(<citation>,
bibtex=TRUE)', 'toBibtex(.)', or set
'options(citation.bibtex.max=999)'.
```

## News

Changelog can be view [here](NEWS.md).

  - \[2016-06-02\] A video of the **netdiffuseR** workshop at SUNBELT
    2016 is now online on
    [youtube](https://www.youtube.com/playlist?list=PLT-GgRN1lFI4coHDqkRJm3flDw9e1gg2P),
    and the workshop materials can be found
    [here](https://github.com/USCCANA/netdiffuser-sunbelt2016/)
  - \[2016-04-11\] **netdiffuseR** will be on
    [useR\! 2016](https://user2016.r-project.org/) on as a presentation
    and on
    [IC2S2 2016](https://www.kellogg.northwestern.edu/news-events/conference/ic2s2/2016.aspx)
    in the posters session.
  - \[2016-03-16\] Next CRAN release scheduled for April 11th 2016
    (after the workshop).
  - \[2016-02-18\] **netdiffuseR** vers 1.16.2 is now on CRAN\!
  - \[2016-02-18\] We will be offering a workshop at the 2016 SUNBELT
    Conference (checkout the full list of workshops
    [here](http://sunbelt2016.insna.org/program/)).

## Installation

### CRAN version

To get the CRAN (stable) version of the package, simple type

``` r
install.packages("netdiffuseR")
```

### Bleeding edge version

If you want the latest (unstable) version of **netdiffuseR**, using the
`devtools` package, you can install `netdiffuseR` dev version as follows

``` r
devtools::install_github('USCCANA/netdiffuseR', build_vignettes = TRUE)
```

You can skip building vignettes by setting `build_vignettes = FALSE` (so
it is not required).

For the case of OSX users, there seems to be a problem when installing
packages depending on `Rcpp`. This issue, developed
[here](https://github.com/USCCANA/netdiffuseR/issues/3), can be solved
by open the terminal and typing the following

``` sh
curl -O http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2
sudo tar fvxz gfortran-4.8.2-darwin13.tar.bz2 -C /
```

before installing the package through `devtools`.

### Binary versions

For the case of windows and mac users, they can find binary versions of
the package [here](https://github.com/USCCANA/netdiffuseR/releases),
netdiffuseR\_1…zip, and netdiffuseR\_1…tgz respectively. They can
install this directly as follows (using the 1.16.3.29 version):

1.  Install dependencies from CRAN
    
    ``` r
    > install.packages(c("igraph", "Matrix", "SparseM", "RcppArmadillo", "sna"), dependencies=TRUE)
    ```

2.  Download the binary version and install it as follows:
    
    ``` r
    > install.packages("netdiffuseR_1.16.3.29.zip", repos=NULL)
    ```
    
    For windows users, and for Mac users:
    
    ``` r
    > install.packages("netdiffuseR_1.16.3.29.tgz", repos=NULL)
    ```

## Tutorials

Since starting netdiffuseR, we have done a couple of workshops at
Sunbelt and NASN. Here are the repositories:

  - Sunbelt 2018: <https://usccana.github.io/netdiffuser-sunbelt2018/>
    ([source code](https://github.com/USCCANA/netdiffuser-sunbelt2018))
  - NASN 2017: <https://usccana.github.io/netdiffuser-nasn2017/>
    ([source code](https://github.com/USCCANA/netdiffuser-nasn2017))
  - Sunbelt 2016: <https://github.com/USCCANA/netdiffuser-sunbelt2016>

## Presentations

  - ic2s2 2016 Evanston, IL:
    <https://github.com/USCCANA/netdiffuser-ic2s22016> (poster)
  - useR\! 2016 Stanford, CA:
    <https://github.com/USCCANA/netdiffuser-user2016> (slides)
  - useR\! 2016: <https://github.com/USCCANA/netdiffuser-user2016>

## Examples

This example has been taken from the package’s vignettes:

``` r
library(netdiffuseR)
```

    ## 
    ## Attaching package: 'netdiffuseR'

    ## The following object is masked from 'package:base':
    ## 
    ##     %*%

### Infectiousness and Susceptibility

``` r
# Generating a random graph
set.seed(1234)
n <- 100
nper <- 20
graph <- rgraph_er(n, nper, .5)
toa <- sample(c(1:(1+nper-1), NA), n, TRUE)
head(toa)
```

    ## [1] 16  3 14  3 13  5

``` r
# Creating a diffnet object
diffnet <- as_diffnet(graph, toa)
diffnet
```

    ## Dynamic network of class -diffnet-
    ##  Name               : Diffusion Network
    ##  Behavior           : Unspecified
    ##  # of nodes         : 100 (1, 2, 3, 4, 5, 6, 7, 8, ...)
    ##  # of time periods  : 20 (1 - 20)
    ##  Type               : directed
    ##  Final prevalence   : 0.95
    ##  Static attributes  : -
    ##  Dynamic attributes : -

``` r
summary(diffnet)
```

    ## Diffusion network summary statistics
    ## Name     : Diffusion Network
    ## Behavior : Unspecified
    ## -----------------------------------------------------------------------------
    ##  Period   Adopters   Cum Adopt. (%)   Hazard Rate   Density   Moran's I (sd)  
    ## -------- ---------- ---------------- ------------- --------- ---------------- 
    ##        1          8         8 (0.08)             -      0.50 -0.01 (0.00)     
    ##        2          3        11 (0.11)          0.03      0.50 -0.01 (0.00)     
    ##        3          6        17 (0.17)          0.07      0.51 -0.01 (0.00)     
    ##        4          3        20 (0.20)          0.04      0.49 -0.01 (0.00)     
    ##        5          9        29 (0.29)          0.11      0.50 -0.01 (0.00)     
    ##        6          5        34 (0.34)          0.07      0.50 -0.01 (0.00)     
    ##        7          2        36 (0.36)          0.03      0.51 -0.01 (0.00)     
    ##        8          3        39 (0.39)          0.05      0.50 -0.01 (0.00)     
    ##        9          5        44 (0.44)          0.08      0.50 -0.01 (0.00)     
    ##       10          1        45 (0.45)          0.02      0.49 -0.01 (0.00)     
    ##       11          3        48 (0.48)          0.05      0.50 -0.01 (0.00)     
    ##       12          6        54 (0.54)          0.12      0.50 -0.01 (0.00)     
    ##       13          8        62 (0.62)          0.17      0.50 -0.01 (0.00)     
    ##       14          9        71 (0.71)          0.24      0.50 -0.01 (0.00)     
    ##       15          5        76 (0.76)          0.17      0.50 -0.00 (0.00) **  
    ##       16          7        83 (0.83)          0.29      0.50 -0.01 (0.00)     
    ##       17          5        88 (0.88)          0.29      0.49 -0.00 (0.00) *** 
    ##       18          4        92 (0.92)          0.33      0.50 -0.01 (0.00)     
    ##       19          1        93 (0.93)          0.12      0.50 -0.01 (0.00)     
    ##       20          2        95 (0.95)          0.29      0.50 -0.01 (0.00)     
    ## -----------------------------------------------------------------------------
    ##  Left censoring  : 0.08 (8)
    ##  Right centoring : 0.05 (5)
    ##  # of nodes      : 100
    ## 
    ##  Moran's I was computed on contemporaneous autocorrelation using 1/geodesic
    ##  values. Significane levels  *** <= .01, ** <= .05, * <= .1.

``` r
# Visualizing distribution of suscep/infect
out <- plot_infectsuscep(diffnet, bins = 20,K=5, logscale = FALSE, h=.01)
```

![](README_files/figure-gfm/plot_infectsuscept-1.png)<!-- -->

``` r
out <- plot_infectsuscep(diffnet, bins = 20,K=5, logscale = TRUE,
                         exclude.zeros = TRUE, h=1)
```

    ## Warning in plot_infectsuscep.list(graph$graph, graph$toa, t0, normalize, : When
    ## applying logscale some observations are missing.

![](README_files/figure-gfm/plot_infectsuscept-2.png)<!-- -->

### Threshold

``` r
# Generating a random graph
set.seed(123)
diffnet <- rdiffnet(500, 20,
                    seed.nodes = "random",
                    rgraph.args = list(m=3),
                    threshold.dist = function(x) runif(1, .3, .7))
diffnet
```

    ## Dynamic network of class -diffnet-
    ##  Name               : A diffusion network
    ##  Behavior           : Random contagion
    ##  # of nodes         : 500 (1, 2, 3, 4, 5, 6, 7, 8, ...)
    ##  # of time periods  : 20 (1 - 20)
    ##  Type               : directed
    ##  Final prevalence   : 1.00
    ##  Static attributes  : real_threshold (1)
    ##  Dynamic attributes : -

``` r
# Threshold with fixed vertex size
plot_threshold(diffnet)
```

![](README_files/figure-gfm/BoringThreshold,%20plot_threshold-1.png)<!-- -->

Using more features

``` r
data("medInnovationsDiffNet")
set.seed(131)
plot_threshold(
  medInnovationsDiffNet,
  vertex.color     = viridisLite::inferno(4)[medInnovationsDiffNet[["city"]]],
  vertex.sides     = medInnovationsDiffNet[["city"]] + 2,
  sub = "Note: Vertices' sizes and shapes given by degree and city respectively",
  jitter.factor = c(1,1), jitter.amount = c(.25,.025)
)
```

    ## Warning in (function (graph, expo, toa, include_censored = FALSE, t0 =
    ## min(toa, : -vertex.sides- will be coerced to integer.

![](README_files/figure-gfm/NiceThreshold-1.png)<!-- -->

### Adoption rate

``` r
plot_adopters(diffnet)
```

![](README_files/figure-gfm/Adopters-1.png)<!-- -->

### Hazard rate

``` r
hazard_rate(diffnet)
```

![](README_files/figure-gfm/Hazard-1.png)<!-- -->

### Diffusion process

``` r
plot_diffnet(medInnovationsDiffNet, slices=c(1,9,8))
```

![](README_files/figure-gfm/plot_diffnet-1.png)<!-- -->

``` r
diffnet.toa(brfarmersDiffNet)[brfarmersDiffNet$toa >= 1965] <- NA
plot_diffnet2(brfarmersDiffNet, vertex.size = "indegree")
```

![](README_files/figure-gfm/plot_diffnet2-1.png)<!-- -->

``` r
set.seed(1231)

# Random scale-free diffusion network
x <- rdiffnet(1000, 4, seed.graph="scale-free", seed.p.adopt = .025,
                           rewire = FALSE, seed.nodes = "central",
                           rgraph.arg=list(self=FALSE, m=4),
                           threshold.dist = function(id) runif(1,.2,.4))

# Diffusion map (no random toa)
dm0 <- diffusionMap(x, kde2d.args=list(n=150, h=1), layout=igraph::layout_with_fr)

# Random
diffnet.toa(x) <- sample(x$toa, size = nnodes(x))

# Diffusion map (random toa)
dm1 <- diffusionMap(x, layout = dm0$coords, kde2d.args=list(n=150, h=.5))

oldpar <- par(no.readonly = TRUE)
col <- viridisLite::plasma(100)
par(mfrow=c(1,2), oma=c(1,0,0,0), cex=.8)
image(dm0, col=col, main="Non-random Times of Adoption\nAdoption from the core.")
image(dm1, col=col, main="Random Times of Adoption")
par(mfrow=c(1,1))
mtext("Both networks have the same distribution on times of adoption", 1,
      outer = TRUE)
```

![](README_files/figure-gfm/plot_diffnet2%20with%20map-1.png)<!-- -->

``` r
par(oldpar)
```

### Adopters classification

``` r
out <- classify(kfamilyDiffNet, include_censored = TRUE)
ftable(out)
```

    ##                thr Non-Adopters Very Low Thresh. Low Thresh. High Thresh. Very High Thresh.
    ## toa                                                                                        
    ## Non-Adopters               0.00             0.00        0.00         0.00              0.00
    ## Early Adopters             0.00            14.04        8.40         0.57              0.29
    ## Early Majority             0.00             5.64       11.65         5.54              2.58
    ## Late Majority              0.00             1.34        5.06         6.21              2.96
    ## Laggards                   0.00             1.53        0.00         0.00             34.19

``` r
# Plotting 
oldpar <- par(no.readonly = TRUE)
par(xpd=TRUE)
plot(out, color=viridisLite::inferno(5), las = 2, xlab="Time of Adoption",
     ylab="Threshold", main="")

# Adding key
legend("bottom", legend = levels(out$thr), fill=viridisLite::inferno(5), horiz = TRUE,
       cex=.6, bty="n", inset=c(0,-.1))
```

![](README_files/figure-gfm/mosaic-1.png)<!-- -->

``` r
par(oldpar)
```

### Session info

``` r
sessionInfo()
```

    ## R version 4.1.0 (2021-05-18)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 18.04.5 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/atlas/libblas.so.3.10.3
    ## LAPACK: /usr/lib/x86_64-linux-gnu/atlas/liblapack.so.3.10.3
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] netdiffuseR_1.22.2
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] igraph_1.2.6          Rcpp_1.0.6            knitr_1.33           
    ##  [4] magrittr_2.0.1        networkDynamic_0.10.1 network_1.16.1       
    ##  [7] MASS_7.3-54           viridisLite_0.4.0     lattice_0.20-44      
    ## [10] rlang_0.4.11          fansi_0.4.2           highr_0.9            
    ## [13] stringr_1.4.0         tools_4.1.0           MatchIt_4.1.0        
    ## [16] grid_4.1.0            rle_0.9.2             xfun_0.23            
    ## [19] utf8_1.2.1            sna_2.6               coda_0.19-4          
    ## [22] htmltools_0.5.1.1     ellipsis_0.3.2        yaml_2.2.1           
    ## [25] digest_0.6.27         tibble_3.1.2          lifecycle_1.0.0      
    ## [28] crayon_1.4.1          Matrix_1.3-3          vctrs_0.3.8          
    ## [31] evaluate_0.14         statnet.common_4.4.1  rmarkdown_2.8        
    ## [34] stringi_1.6.2         compiler_4.1.0        pillar_1.6.1         
    ## [37] backports_1.2.1       boot_1.3-28           SparseM_1.81         
    ## [40] pkgconfig_2.0.3

## To-do list

  - Import/Export functions for interfacing other package’s clases, in
    particular: `statnet` set (specially the packages `networkDynamic`
    and `ndtv`), ~~`igraph`~~ and `Rsiena`.
  - Populate the tests folder.
  - ~~Use spells? (`select_egoalter` would use this)~~
  - ~~Classify individuals by adoption category using early adopters,
    adopters, and laggards, and by threshold using very low, low, high
    and very high threshold (Valente 95’ p. 94).~~
  - ~~Double check all functions using adjacency matrix values.~~
  - ~~Remove dimnames from matrices and vectors. It is more efficient to
    use the ones stored in meta instead.~~
  - Implement the Bass model
  - ~~Include function to import survey data (as shown on the
    vignettes)~~
  - Exposure based on Mahalanobis distances and also Roger Leenders on
    weighting exposure (internal note).
  - (2016-03-30): use `xspline` for drawing polygons & edges.
  - ~~(2016-04-04): Add more options to `exposure`, namely, `self` (so
    removes diagonal or not\!).~~
  - (2016-04-19): animal behaviorists.
  - (2016-10-18): Review language throughout the manual (more than
    innovation).
  - (2016-10-18): Evaluate and eventually use a standard graph format
    (`network` for instance?).
  - (2016-10-18): Standarize graph plot methods (choose either
    statnet/igraph/own)
