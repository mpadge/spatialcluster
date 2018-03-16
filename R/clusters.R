#' scl_cluster
#'
#' Cluster spatial data
#'
#' @param xy Rectangular structure (matrix, data.frame, tibble), containing
#' coordinates of points to be clustered.
#' @param dmat Square structure (matrix, data.frame, tibble) containing
#' distances or equivalent metrics betwen all points in \code{xy}. If \code{xy}
#' has \code{n} rows, then \code{dat} must have \code{n} rows and \code{n}
#' columns.
#' @param ncl Desired number of clusters
#' @param shortest If \code{TRUE}, the \code{dmat} is interpreted as distances
#' such that lower values are preferentially selected; if \code{FALSE}, then
#' higher values of \code{dmat} are interpreted to indicate stronger
#' relationships, as is the case for example with covariances.
#'
#' @return A object of class \code{scl} with \code{tree} containing the
#' clustering scheme, and \code{xy} the original coordinate data of the
#' clustered points.
#' @export
#' @examples
#' \dontrun{
#' xy <- matrix (runif (100), ncol = 2)
#' edges <- scl_edges (xy)
#' # add some fake data to the edges
#' edges %<>% dplyr::mutate (d = runif (nrow (.))) %>%
#'    dplyr::arrange (desc (d))
#' tree <- scl_spantree (edges) # plain tree; no components
#' ncl <- 8 # desired number of clusters/components
#' tree <- scl_cuttree (tree, edges, ncl = ncl) # tree with component numbers
#' }
scl_cluster <- function (xy, dmat, ncl, shortest = TRUE)
{
    xy <- scl_tbl (xy)
    edges <- scl_edges (xy, dmat, shortest)
    tree_full <- scl_spantree (edges)
    trees <- scl_cuttree (tree_full, edges, ncl)

    tree <- trees$tree_in
    tree$id <- seq (nrow (tree))
    cmps <- rcpp_get_component_vector (tree)
    tree <- tibble::tibble (id = as.numeric (cmps$edge_id),
                    comp = cmps$edge_component) %>%
        dplyr::arrange (id) %>%
        dplyr::left_join (tree, ., by = "id") %>%
        dplyr::select (from, to, d, comp)

    structure (list (xy = xy, tree = tree, tree_rest = trees$tree_out),
               class = "scl")
}
