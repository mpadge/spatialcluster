#' scl_redcap
#'
#' Cluster spatial data with REDCAP (REgionalization with Dynamically
#' Constrained Agglomerative clustering and Partitioning) routines.
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
#' @param nnbs Number of nearest neighbours to be used in calculating clustering
#' trees. Triangulation will be used if \code{nnbs <= 0}.
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
#' first-order relationships. It is therefore strongly recommended that the
#' default \code{full_order = TRUE} be used at all times.
#'
#' @export
#' @examples
#' n <- 100
#' xy <- matrix (runif (2 * n), ncol = 2)
#' dmat <- matrix (runif (n ^ 2), ncol = n)
#' scl <- scl_redcap (xy, dmat, ncl = 4)
#' # Those clusters will by default be constructed by connecting edges with the
#' # lowest (\code{shortest}) values of \code{dmat}, and will differ from
#' scl <- scl_redcap (xy, dmat, ncl = 4, shortest = FALSE)
#' # using 'full_order = FALSE' constructs clusters from first-order
#' # relationships only; not recommended, but possible nevertheless:
#' scl <- scl_redcap (xy, dmat, ncl = 4, full_order = FALSE)
scl_redcap <- function (xy,
                        dmat,
                        ncl,
                        full_order = TRUE,
                        linkage = "single",
                        shortest = TRUE,
                        nnbs = 6) {
    linkage <- scl_linkage_type (linkage)

    if (methods::is (xy, "scl")) {
        if (!identical (xy$pars$method, "redcap"))
            stop ("scl_redcap can pass to scl_recluster only for scl objects",
                  " previously generated with scl_redcap")

        message ("scl_redcap is for initial cluster construction; ",
                 "passing to scl_recluster")
        scl_recluster_redcap (xy, ncl = ncl, shortest = shortest)
    } else {
        xy <- scl_tbl (xy)
        if (nnbs <= 0)
            edges_nn <- scl_edges_tri (xy, dmat, shortest)
        else
            edges_nn <- scl_edges_nn (xy, dmat, nnbs, shortest)
        if (!full_order) {
            tree_full <- scl_spantree_ord1 (edges_nn)
        } else {
            if (linkage == "average") {
                tree_full <- scl_spantree_alk (edges_nn, shortest)
            } else {
                edges_all <- scl_edges_all (xy, dmat, shortest)
                if (linkage == "single") {
                    tree_full <- scl_spantree_slk (edges_all, edges_nn,
                                                   shortest)
                } else if (linkage == "complete") {
                    tree_full <- scl_spantree_clk (edges_all, edges_nn,
                                                   shortest)
                } else {
                    stop ("linkage must be one of ",
                          "(single, average, complete)")
                }
            }

        }

        tree <- scl_cuttree (tree_full, edges_nn, ncl, shortest)

        # meta-data:
        clo <- c ("single", "full") [match (full_order, c (FALSE, TRUE))]
        pars <- list (method = "redcap",
                      ncl = ncl,
                      cl_order = clo,
                      linkage = linkage)

        structure (list (tree = tree,
                         nodes = dplyr::bind_cols (tree_nodes (tree), xy),
                         pars = pars),
                   class = "scl")
    }
}

# Match cluster numbers in edge tree to actual nodes
tree_nodes <- function (tree) {
    node <- NULL # suppress no visible binding note
    res <- tibble::tibble (node = c (tree$from, tree$to),
                           cluster = rep (tree$cluster, 2)) %>%
        dplyr::distinct () %>%
        dplyr::arrange (node) %>%
        dplyr::filter (!is.na (cluster))
    # remove clusters with < 3 members:
    res$cluster [res$cluster %in% which (table (res$cluster) < 3)] <- NA
    return (res)
}

#' scl_reccluster
#'
#' Re-cut a spatial cluster tree (\code{scl}) at a different number of clusters.
#'
#' @param scl An \code{scl} object returned from \link{scl_redcap}.
#' @inheritParams scl_redcap
#'
#' @return Modified \code{scl} object in which \code{tree} is re-cut into
#' \code{ncl} clusters.
#' @export
#' @examples
#' n <- 100
#' xy <- matrix (runif (2 * n), ncol = 2)
#' dmat <- matrix (runif (n ^ 2), ncol = n)
#' scl <- scl_redcap (xy, dmat, ncl = 4)
#' plot (scl)
#' scl <- scl_recluster (scl, ncl = 5)
#' plot (scl)
scl_recluster <- function (scl, ncl, shortest = TRUE) {
    if (!methods::is (scl, "scl"))
        stop ("scl_recluster can only be applied to 'scl' objects ",
              "returned from scl_redcap")
    else if (identical (scl$pars$method, "redcap"))
        scl_recluster_redcap (scl = scl, ncl = ncl, shortest = shortest)
    else if (identical (scl$pars$method, "full"))
        scl_recluster_full (scl = scl, ncl = ncl)
}

scl_recluster_redcap <- function (scl, ncl, shortest = TRUE) {
    from <- to <- d <- NULL # no visible binding messages
    tree_full <- scl$tree %>% dplyr::select (from, to, d)
    if (shortest)
        tree_full %<>% dplyr::arrange (d)
    else
        tree_full %<>% dplyr::arrange (dplyr::desc (d))

    tree_full$cluster <- rcpp_cut_tree (tree_full, ncl,
                                        shortest = shortest) + 1

    pars <- scl$pars
    pars$ncl <- ncl

    structure (list (tree = tree_full,
                     nodes = dplyr::bind_cols (tree_nodes (tree_full),
                                               scl$nodes [, c ("x", "y")]),
                     pars = pars),
               class = "scl")
}
