#' scl_exact
#'
#' Exact spatially-constrained clustering.
#'
#' @param xy Rectangular structure (matrix, data.frame, tibble), containing
#' coordinates of points to be clustered.
#' @param dmat Square structure (matrix, data.frame, tibble) containing
#' distances or equivalent metrics betwen all points in \code{xy}. If \code{xy}
#' has \code{n} rows, then \code{dat} must have \code{n} rows and \code{n}
#' columns.
#' @param ncl Desired number of clusters
#'
#' @return A object of class \code{scl} with \code{tree} containing the
#' clustering scheme, and \code{xy} the original coordinate data of the
#' clustered points. An additional component, \code{tree_rest}, enables the tree
#' to be re-cut to a different number of clusters via \link{scl_recluster},
#' rather than calculating clusters anew.
#'
#'
#' @export
scl_exact <- function (xy, dmat, ncl)
{
    xy <- scl_tbl (xy)
    edges <- scl_edges_nn (xy, dmat, shortest = TRUE)
    # cluster numbers can be joined with edges through either from or to:
    cl <- tibble::tibble (from = seq (nrow (xy)),
                          cl = rcpp_exact_initial (edges) + 1)
    edges <- dplyr::left_join (edges, cl, by = c ("from"))
    # edges then include initial cluster allocation

    edges <- rcpp_exact_merge (edges, ncl = ncl)
}
