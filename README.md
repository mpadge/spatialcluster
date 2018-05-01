<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Build Status](https://travis-ci.org/mpadge/spatialcluster.svg)](https://travis-ci.org/mpadge/spatialcluster) [![Project Status: WIP](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip) [![codecov](https://codecov.io/gh/mpadge/spatialcluster/branch/master/graph/badge.svg)](https://codecov.io/gh/mpadge/spatialcluster)

spatialcluster
==============

Clustering algorithms are borne of historical computational necessity. Approximations were sought to boost efficiency and thereby enable clusters to be computed at all. This package includes an algorithm to *exactly* compute spatially constrained clusters. The algorithm is as efficient as can be, but is of course slower than approximate algorithms. It is nevertheless exact.

''Spatially-constrained'' means that the data from which clusters are to be formed also map on to spatial reference points, and the constraint is that clusters must be spatially contiguous. The actual clustering data are generally in the form of a distance or similarity matrix quantifying relationships between some collection of objects or points, as with most clustering procedures, while the spatial data are simply coordinates of those objects or points.

For comparison, this package also includes the REDCAP collection of efficient yet approximate algorithms described in [D. Guo's 2008 paper, "Regionalization with dynamically constrained agglomerative clustering and partitioning."](https://www.tandfonline.com/doi/abs/10.1080/13658810701674970) (pdf available [here](https://pdfs.semanticscholar.org/ead1/7df8aaa1aed0e433b3ae1ec1ec5c7e785b2b.pdf)).

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
n <- 100
xy <- matrix (runif (2 * n), ncol = 2)
dmat <- matrix (runif (n ^ 2), ncol = n)
```

The load the package and call the function:

``` r
library (spatialcluster)
scl <- scl_cluster (xy, dmat, ncl = 8, linkage = "single")
plot (scl)
```

![](/docs/figs/README-plot-single-1.png)

``` r
scl <- scl_cluster (xy, dmat, ncl = 8, linkage = "average")
plot (scl)
```

![](/docs/figs/README-plot-average-1.png)

``` r
scl <- scl_cluster (xy, dmat, ncl = 8, linkage = "complete")
plot (scl)
```

![](/docs/figs/README-plot-complete-1.png)

This example illustrates the universal danger in all clustering algorithms: they can not fail to produce results, even when the data fed to them are definitely devoid of any information as in this example. Clustering algorithms should only be applied to reflect a very specific hypothesis for why data should be clustered in the first place; spatial clustering algorithms should only be applied to reflect two very specific hypothesis for (i) why data should be clustered at all, and (ii) why those clusters should manifest a spatial pattern.
