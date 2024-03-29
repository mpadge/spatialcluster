% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/redcap.R
\name{scl_redcap}
\alias{scl_redcap}
\title{scl_redcap}
\usage{
scl_redcap(
  xy,
  dmat,
  ncl,
  full_order = TRUE,
  linkage = "single",
  shortest = TRUE,
  nnbs = 6L,
  iterate_ncl = FALSE,
  quiet = FALSE
)
}
\arguments{
\item{xy}{Rectangular structure (matrix, data.frame, tibble), containing
coordinates of points to be clustered.}

\item{dmat}{Square structure (matrix, data.frame, tibble) containing
distances or equivalent metrics between all points in \code{xy}. If \code{xy}
has \code{n} rows, then \code{dat} must have \code{n} rows and \code{n}
columns.}

\item{ncl}{Desired number of clusters. See description of `ncl_iterate`
parameter for conditions under which actual number may be less than this
value.}

\item{full_order}{If \code{FALSE}, build spanning trees from first-order
relationships only, otherwise build from full-order relationships (see Note).}

\item{linkage}{One of \code{"single"}, \code{"average"}, or
\code{"complete"}; see Note.}

\item{shortest}{If \code{TRUE}, the \code{dmat} is interpreted as distances
such that lower values are preferentially selected; if \code{FALSE}, then
higher values of \code{dmat} are interpreted to indicate stronger
relationships, as is the case for example with covariances.}

\item{nnbs}{Number of nearest neighbours to be used in calculating clustering
trees. Triangulation will be used if \code{nnbs <= 0}.}

\item{iterate_ncl}{Actual numbers of clusters found may be less than the
specified value of `ncl`, because clusters formed from < 3 edges are removed.
If `iterate_ncl = FALSE` (the default), the value is returned with whatever
number of actual clusters is found. Setting this parameter to `TRUE` forces
the algorithm to iterate until the exact number of clusters has been found.
For large data sets, this may result in considerable longer calculation
times.}

\item{quiet}{If `FALSE` (default), display progress information on screen.}
}
\value{
A object of class \code{scl} with \code{tree} containing the
clustering scheme, and \code{xy} the original coordinate data of the
clustered points. An additional component, \code{tree_rest}, enables the tree
to be re-cut to a different number of clusters via \link{scl_recluster},
rather than calculating clusters anew.
}
\description{
Cluster spatial data with REDCAP (REgionalization with Dynamically
Constrained Agglomerative clustering and Partitioning) routines.
}
\note{
Please refer to the original REDCAP paper ('Regionalization with
dynamically constrained agglomerative clustering and partitioning (REDCAP)',
by D. Guo (2008), Int.J.Geo.Inf.Sci 22:801-823) for details of the
\code{full_order} and \code{linkage} parameters. This paper clearly
demonstrates the general inferiority of spanning trees constructed from
first-order relationships. It is therefore strongly recommended that the
default \code{full_order = TRUE} be used at all times.
}
\examples{
n <- 100
xy <- matrix (runif (2 * n), ncol = 2)
dmat <- matrix (runif (n ^ 2), ncol = n)
scl <- scl_redcap (xy, dmat, ncl = 4)
# Those clusters will by default be constructed by connecting edges with the
# lowest (\code{shortest}) values of \code{dmat}, and will differ from
scl <- scl_redcap (xy, dmat, ncl = 4, shortest = FALSE)
# using 'full_order = FALSE' constructs clusters from first-order
# relationships only; not recommended, but possible nevertheless:
scl <- scl_redcap (xy, dmat, ncl = 4, full_order = FALSE)
}
\seealso{
Other clustering_fns: 
\code{\link{scl_full}()},
\code{\link{scl_recluster}()}
}
\concept{clustering_fns}
