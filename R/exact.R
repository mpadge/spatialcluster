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
#' @param linkage One of \code{"single"}, \code{"average"}, or
#' \code{"complete"}; see \link{scl_redcap} for details.
#'
#' @return A object of class \code{scl} with \code{tree} containing the
#' clustering scheme, and \code{xy} the original coordinate data of the
#' clustered points. An additional component, \code{tree_rest}, enables the tree
#' to be re-cut to a different number of clusters via \link{scl_recluster},
#' rather than calculating clusters anew.
#'
#' @export
#' @examples
#' n <- 100
#' xy <- matrix (runif (2 * n), ncol = 2)
#' dmat <- matrix (runif (n ^ 2), ncol = n)
#' scl <- scl_exact (xy, dmat, ncl = 4)
scl_exact <- function (xy, dmat, ncl, linkage = "single")
{
    linkage <- scl_linkage_type (linkage)

    xy <- scl_tbl (xy)
    edges <- scl_edges_nn (xy, dmat, shortest = TRUE)
    # cluster numbers can be joined with edges through either from or to:
    cl <- rcpp_exact_initial (edges) + 1

    # make 3 vectors of cluster numbers:
    #   1. cl = cluster number for intra-cluster edges only;
    #   2. cl_from = Number of origin cluster for inter-cluster edges only; and
    #   3. cl_to = Number of destination cluster for inter-cluster edges only.
    from_cl <- cl [edges$from]
    to_cl <- cl [edges$to]
    indx <- which (from_cl == to_cl)
    cl_in <- cl_join_from <- cl_join_to <- rep (NA, nrow (edges))
    cl_in [indx] <- from_cl [indx]
    indx <- which (from_cl != to_cl)
    cl_join_from [indx] <- from_cl [indx]
    cl_join_to [indx] <- to_cl [indx]

    edges$cl <- cl_in - 1 # convert back to C++ 0-indexed values
    edges$cl_from <- cl_join_from - 1
    edges$cl_to <- cl_join_to - 1
    edges$cl [is.na (edges$cl)] <- -1
    edges$cl_from [is.na (edges$cl_from)] <- -1
    edges$cl_to [is.na (edges$cl_to)] <- -1

    merges <- rcpp_exact_merge (edges, linkage = linkage) %>%
        data.frame ()

    merges <- tibble::tibble (from = as.integer (merges$from),
                              to = as.integer (merges$to),
                              dist = merges$dist)
    # exact_cluster_nodes just auto-merges the tree to the specified number, but
    # some of these may be clusters with only 2 members. These are excluded here
    # by iterating until the desired number is achieved in which each cluster
    # has >= 3 members:
    num_nodes <- 0
    ncl_trial <- ncl
    while (num_nodes < ncl)
    {
        nodes <- exact_cluster_nodes (edges, merges, ncl_trial)
        num_nodes <- length (which (table (nodes$cluster) > 2))
        ncl_trial <- ncl_trial + 1
        if (ncl_trial >= nrow (nodes))
            break
    }
    n <- which (table (nodes$cluster) == 2)
    nodes$cluster [nodes$cluster %in% n] <- NA

    pars <- list (method = "exact",
                  ncl = ncl,
                  linkage = linkage)


    structure (list (merges = merges,
                     ord = order_merges (merges),
                     nodes = dplyr::bind_cols (nodes, xy),
                     pars = pars),
               class = "scl")
}

#' order_merges
#'
#' Order merges so they can be plotted as dendrogram
#' @param merges output from rccp_exact_merge
#' @noRd
order_merges <- function (merges)
{
    merges <- as.matrix (merges)
    nodes <- merges [nrow (merges), c ("from", "to")]
    for (i in rev (seq (nrow (merges))) [-1])
    {
        ii <- which (nodes == merges [i, 2])
        n1 <- n2 <- NULL
        if (ii > 1)
            n1 <- nodes [1:(ii - 1)]
        if (ii <= length (nodes))
            n2 <- nodes [ii:length (nodes)]
        nodes <- c (n1, merges [i, 1], n2)
    }
    return (as.numeric (nodes))
}

#' exact_cluster_nodes
#'
#' Transform edge and merge data into rectangle of nodes and cluster IDs
#' @noRd
exact_cluster_nodes <- function (edges, merges, ncl)
{
    edges$cl [edges$cl < 0] <- NA
    ncl_exact <- length (unique (edges$cl))
    merge_tree <- merges [1:(ncl_exact - ncl - 1), ]
    for (i in seq (nrow (merge_tree)))
        edges$cl [edges$cl == merge_tree$from [i]] <- merge_tree$to [i]

    node <- cluster <- NULL # rm undefined variable note
    nodes <- tibble::tibble (node = c (edges$from, edges$to),
                             cluster = rep (edges$cl, 2)) %>%
        dplyr::distinct () %>%
        dplyr::arrange (node) %>%
        dplyr::filter (!is.na (cluster))

    # re-order cluster numbers by frequencies
    nt <- sort (table (nodes$cluster), decreasing = TRUE)
    nodes$cluster <- match (nodes$cluster, names (nt))

    return (nodes)
}
