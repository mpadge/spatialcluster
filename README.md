<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Build Status](https://travis-ci.org/mpadge/spatialcluster.svg)](https://travis-ci.org/mpadge/spatialcluster) [![Project Status: Concept - Minimal or no implementation has been done yet.](http://www.repostatus.org/badges/0.1.0/concept.svg)](http://www.repostatus.org/#concept) [![codecov](https://codecov.io/gh/mpadge/spatialcluster/branch/master/graph/badge.svg)](https://codecov.io/gh/mpadge/spatialcluster)

spatialcluster
==============

R port of redcap spatially-constrained clustering, as described in [D. Guo's 2008 paper, "Regionalization with dynamically constrained agglomerative clutering and partitioning."](https://www.tandfonline.com/doi/abs/10.1080/13658810701674970) (pdf available [here](https://pdfs.semanticscholar.org/ead1/7df8aaa1aed0e433b3ae1ec1ec5c7e785b2b.pdf)). ''Spatially-constrained'' implies that the points which are to be clustered have defined spatial coordinates, while the clustering itself is based on additional data, generally in the form of a distance or similarity matrix quantifying relationships between those points.

Installation
------------

You can install spatialcluster from github with:

``` r
# install.packages("devtools")
devtools::install_github("mpadge/spatialcluster")
```

Usage
-----

The single function `scl_cluster()` requires three arguments:

1.  A rectangular matrix of coordinates of points to be clustered (`n` rows; at least 2 columns);
2.  An `n`-by-`n` square matrix quantifying relationships between those points;
3.  A single value (`ncl`) specifying the desired number of clusters.

Usage can be demonstrated with some simple fake data:

``` r
n <- 20
xy <- matrix (runif (2 * n), ncol = 2)
dmat <- matrix (runif (n ^ 2), ncol = n)
```

The load the package and call the function:

``` r
library (spatialcluster)
scl <- scl_cluster (xy, dmat, ncl = 4)
plot (scl)
```

![](README-plot-1.png)
