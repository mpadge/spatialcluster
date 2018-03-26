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
#' @param full_order If \code{FALSE}, build spanning trees from first-order
#' relationships only, otherwise build from full-order relationships (see Note).
#' @param linkage One of \code{"single"}, \code{"average"}, or
#' \code{"complete"}; see Note.
#' @param shortest If \code{TRUE}, the \code{dmat} is interpreted as distances
#' such that lower values are preferentially selected; if \code{FALSE}, then
#' higher values of \code{dmat} are interpreted to indicate stronger
#' relationships, as is the case for example with covariances.
#'
#' @return A object of class \code{scl} with \code{tree} containing the
#' clustering scheme, and \code{xy} the original coordinate data of the
#' clustered points. An additional component, \code{tree_rest}, enables the tree
#' to be re-cut to a different number of clusters via \link{scl_recluster},
#' rather than calculating clusters anew.
#'
#' @note Please refer to the original REDCAP paper ('Regionalization with
#' dynamically constrained agglomerative clustering and partitioning (REDCAP)',
#' by D. Guo (2008), Int.J.Geo.Inf.Sci 22:801-823) for details of the
#' \code{full_order} and \code{linkage} parameters. This paper clearly
#' demonstrates the general inferiority of spanning trees constructed from
#' first-order relationships. It is therefore strongly recommended that the default
#' \code{full_order = TRUE} be used at all times.
#'
#' @export
#' @examples
#' n <- 100
#' xy <- matrix (runif (2 * n), ncol = 2)
#' dmat <- matrix (runif (n ^ 2), ncol = n)
#' scl <- scl_cluster (xy, dmat, ncl = 4)
#' # Those clusters will by default be constructed by connecting edges with the
#' # lowest (\code{shortest}) values of \code{dmat}, and will differ from
#' scl <- scl_cluster (xy, dmat, ncl = 4, shortest = FALSE)
#' # using 'full_order = FALSE' constructs clusters from first-order
#' # relationships only; not recommended, but possible nevertheless:
#' scl <- scl_cluster (xy, dmat, ncl = 4, full_order = FALSE)
scl_cluster <- function (xy, dmat, ncl, full_order = TRUE, linkage = "single",
                         shortest = TRUE)
{
    linkage <- scl_linkage_type (linkage)

    if (is (xy, "scl"))
    {
        message ("scl_cluster is for initial cluster construction; ",
                 "passing to scl_recluster")
        scl_recluster (xy, ncl = ncl)
    } else
    {
        xy <- scl_tbl (xy)
        edges_nn <- scl_edges_nn (xy, dmat, shortest)
        if (!full_order)
        {
            tree_full <- scl_spantree_O1 (edges_nn)
        } else
        {
            if (linkage == "single")
            {
                edges_all <- scl_edges_all (xy, dmat, shortest)
                tree_full <- scl_spantree_slk (edges_all, edges_nn)
            } else if (linkage == "average")
            {
                tree_full <- scl_spantree_alk (edges_nn)
            } else
            {
                stop ("Only single linkage implemented at the moment.")
            }
        }
        trees <- scl_cuttree (tree_full, edges_nn, ncl)

        tree <- tree_components (trees$tree_in)

        # meta-data:
        pars <- list (ncl = ncl,
                      cl_order = c ("single", "full") [match (full_order,
                                                              c (FALSE, TRUE))],
                      linkage = linkage)

        structure (list (xy = xy, tree = tree,
                         tree_rest = trees$tree_out,
                         pars = pars),
                   class = "scl")
    }
}

tree_components <- function (tree)
{
    tree$id <- seq (nrow (tree))
    cmps <- rcpp_get_component_vector (tree)
    tree <- tibble::tibble (id = as.numeric (cmps$edge_id),
                            comp = cmps$edge_component) %>%
        dplyr::arrange (id) %>%
        dplyr::left_join (tree, ., by = "id") %>%
        dplyr::select (from, to, d, comp)

    return (tree)
}

#' scl_reccluster
#'
#' Re-cut a spatial cluster tree (\code{scl}) at a different number of clusters.
#'
#' @param scl An \code{scl} object returned from \link{scl_cluster}.
#' @param ncl Number of clusters or components into which tree is to be re-cut
#' @param shortest If \code{TRUE}, the \code{dmat} is interpreted as distances
#' such that lower values are preferentially selected; if \code{FALSE}, then
#' higher values of \code{dmat} are interpreted to indicate stronger
#' relationships, as is the case for example with covariances.
#'
#' @return Modified \code{scl} object in which \code{tree} is re-cut into
#' \code{ncl} clusters.
#' @export
#' @examples
#' n <- 100
#' xy <- matrix (runif (2 * n), ncol = 2)
#' dmat <- matrix (runif (n ^ 2), ncol = n)
#' scl <- scl_cluster (xy, dmat, ncl = 4)
#' plot (scl)
#' scl <- scl_recluster (scl, ncl = 5)
#' plot (scl)
scl_recluster <- function (scl, ncl, shortest = TRUE)
{
    if (!is (scl, "scl"))
        stop ("scl_recluster can only be applied to 'scl' objects ",
              "returned from scl_cluster")

    tree_full <- scl$tree %>% dplyr::select (from, to, d) %>%
        dplyr::bind_rows (scl$rest)
    if (shortest)
        tree_full %<>% dplyr::arrange (d)
    else
        tree_full %<>% dplyr::arrange (dplyr::desc (d))

    # cut tree:
    n <- nrow (tree_full)
    tree_rest <- tree_full [1:(n - 1), ]
    tree <- tree_components (tree_full [ncl:n, ])

    structure (list (xy = scl$xy,
                     tree = tree,
                     tree_rest = tree_rest,
                     ncl = ncl),
               class = "scl")
}
