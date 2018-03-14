<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Build Status](https://travis-ci.org/mpadge/spatialcluster.svg)](https://travis-ci.org/mpadge/spatialcluster) [![Project Status: Concept - Minimal or no implementation has been done yet.](http://www.repostatus.org/badges/0.1.0/concept.svg)](http://www.repostatus.org/#concept) [![codecov](https://codecov.io/gh/mpadge/spatialcluster/branch/master/graph/badge.svg)](https://codecov.io/gh/mpadge/spatialcluster)

spatialcluster
==============

R port of redcap spatial clustering, as described in [D. Guo's 2008 paper, "Regionalization with dynamically constrained agglomerative clutering and partitioning."](https://www.tandfonline.com/doi/abs/10.1080/13658810701674970) (pdf available [here](https://pdfs.semanticscholar.org/ead1/7df8aaa1aed0e433b3ae1ec1ec5c7e785b2b.pdf)).

Installation
------------

You can install spatialcluster from github with:

``` r
# install.packages("devtools")
devtools::install_github("mpadge/spatialcluster")
```

Usage
-----

    #> Loading spatialcluster

``` r
library (spatialcluster)
```

``` r
xy <- matrix (runif (100), ncol = 2)
edges <- scl_edges (xy)
# add some fake data to the edges
edges %<>% dplyr::mutate (d = runif (nrow (.))) %>%
   dplyr::arrange (desc (d))
edges
#> # A tibble: 274 x 3
#>     from    to     d
#>    <int> <int> <dbl>
#>  1    44    39 0.996
#>  2    28    32 0.993
#>  3     1     4 0.989
#>  4    37    34 0.989
#>  5    47    27 0.987
#>  6    17    24 0.980
#>  7    14    43 0.978
#>  8    32    19 0.977
#>  9    18    22 0.974
#> 10    32    28 0.970
#> # ... with 264 more rows
# get tree with component numbers
ncl <- 12 # desired number of clusters/components
tree <- scl_spantree (edges) %>%
    scl_cuttree (edges, ncl = ncl)
```

The clusters resulting from the cutting the spanning tree can then be viewed with:

``` r
g <- scl_plot (tree, xy)
```

![](README-plot-1.png)
