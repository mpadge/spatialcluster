#' scl_full
#'
#' Full spatially-constrained clustering.
#'
#' @inheritParams scl_redcap
#'
#' @family clustering_fns
#' @export
#' @examples
#' n <- 100
#' xy <- matrix (runif (2 * n), ncol = 2)
#' dmat <- matrix (runif (n ^ 2), ncol = n)
#' scl <- scl_full (xy, dmat, ncl = 4)
scl_full <- function (xy,
                      dmat,
                      ncl,
                      linkage = "single",
                      shortest = TRUE,
                      nnbs = 6L) {
    linkage <- scl_linkage_type (linkage)
    if (methods::is (xy, "scl")) {
        message ("scl_full is for initial cluster construction; ",
                 "passing to scl_recluster")
        scl_recluster_full (xy, ncl = ncl)
    } else {
        xy <- scl_tbl (xy)
        edges_all <- scl_edges_all (xy, dmat, shortest)

        if (nnbs <= 0) {
            edges <- scl_edges_tri (xy, dmat, shortest = shortest)
        } else {
            edges <- scl_edges_nn (
                xy,
                dmat,
                edges_all = edges_all,
                nnbs = nnbs,
                shortest = shortest
            )
        }

        # cluster numbers can be joined with edges through either from or to:
        cl <- as.integer (rcpp_full_initial (edges, shortest) + 1)

        # make 3 vectors of cluster numbers:
        #   1. cl = cluster number for intra-cluster edges only;
        #   2. cl_from = Num of origin cluster for inter-cluster edges only; and
        #   3. cl_to = Num of destination cluster for inter-cluster edges only.
        from_cl <- cl [edges$from]
        to_cl <- cl [edges$to]
        indx <- which (from_cl == to_cl)
        cl_in <- cl_join_from <- cl_join_to <- rep (NA, nrow (edges))
        cl_in [indx] <- from_cl [indx]
        indx <- which (from_cl != to_cl)
        cl_join_from [indx] <- from_cl [indx]
        cl_join_to [indx] <- to_cl [indx]

        edges$cluster <- cl_in - 1L # convert back to C++ 0-indexed values
        edges$cl_from <- cl_join_from - 1L
        edges$cl_to <- cl_join_to - 1L
        edges$cluster [is.na (edges$cluster)] <- -1L
        edges$cl_from [is.na (edges$cl_from)] <- -1L
        edges$cl_to [is.na (edges$cl_to)] <- -1L

        merges <- rcpp_full_merge (edges, linkage = linkage,
                                    shortest = shortest) %>%
            data.frame ()

        merges <- tibble::tibble (from = as.integer (merges$from),
                                  to = as.integer (merges$to),
                                  dist = merges$dist)
        # full_cluster_nodes just auto-merges the tree to the specified number,
        # but some of these may be clusters with only 2 members. These are
        # excluded here by iterating until the desired number is achieved in
        # which each cluster has >= 3 members:
        num_nodes <- 0
        ncl_trial <- ncl
        while (num_nodes < ncl) {
            nodes <- full_cluster_nodes (edges, merges, ncl_trial)
            num_nodes <- length (which (table (nodes$cluster) > 2))
            ncl_trial <- ncl_trial + 1
            if (ncl_trial >= nrow (nodes))
                break
        }
        n <- which (table (nodes$cluster) <= 2)
        nodes$cluster [nodes$cluster %in% n] <- NA

        # tree at that point has initial cluster numbers which must be
        # re-aligned with clusters from the nodal merges:
        tree <- edges %>% dplyr::select (from, to, d, cluster)
        tree$cluster <- tree$cl_fr <-
                nodes$cluster [match (tree$from, nodes$node)]
        tree$cl_to <- nodes$cluster [match (tree$to, nodes$node)]
        tree$cluster [tree$cl_fr != tree$cl_to] <- NA

        pars <- list (method = "full",
                      ncl = ncl,
                      linkage = linkage)

        res <- structure (list (tree = dplyr::select (tree,
                                                      c (from, to, d, cluster)),
                                merges = merges,
                                ord = order_merges (merges),
                                nodes = dplyr::bind_cols (nodes, xy),
                                pars = pars),
                          class = "scl")

        res <- scl_statistics (res)

        return (res)
    }
}

#' order_merges
#'
#' Order merges so they can be plotted as dendrogram
#' @param merges output from rccp_full_merge
#' @noRd
order_merges <- function (merges) {
    merges <- as.matrix (merges)
    nodes <- merges [nrow (merges), c ("from", "to")]
    for (i in rev (seq (nrow (merges))) [-1]) {
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

#' full_cluster_nodes
#'
#' Transform edge and merge data into rectangle of nodes and cluster IDs
#' @noRd
full_cluster_nodes <- function (edges, merges, ncl) {
    edges$cluster [edges$cluster < 0] <- NA
    ncl_full <- length (unique (edges$cluster))
    merge_tree <- merges [1:(ncl_full - ncl - 1), ]
    for (i in seq (nrow (merge_tree)))
        edges$cluster [edges$cluster == merge_tree$from [i]] <-
            merge_tree$to [i]

    node <- cluster <- NULL # rm undefined variable note
    nodes <- tibble::tibble (node = c (edges$from, edges$to),
                             cluster = rep (edges$cluster, 2)) %>%
        dplyr::distinct () %>%
        dplyr::arrange (node) %>%
        dplyr::filter (!is.na (cluster))

    # nodes can still be in multiple clusters, so these are set to NA
    dup_nodes <- unique (nodes$node [which (duplicated (nodes$node))])
    nodes$cluster [nodes$node %in% dup_nodes] <- NA_integer_
    nodes <- nodes [which (!duplicated (nodes)), ]

    # re-order cluster numbers by frequencies
    nt <- sort (table (nodes$cluster), decreasing = TRUE)
    nodes$cluster <- match (nodes$cluster, names (nt))

    return (nodes)
}

#' scl_recluster_full
#'
#' @noRd
scl_recluster_full <- function (scl, ncl = ncl) {
    xy <- scl$nodes %>% dplyr::select (x, y)
    num_nodes <- 0
    ncl_trial <- ncl
    while (num_nodes < ncl) {
        scl$nodes <- full_cluster_nodes (scl$tree, scl$merges, ncl_trial)
        num_nodes <- length (which (table (scl$nodes$cluster) > 2))
        ncl_trial <- ncl_trial + 1
        if (ncl_trial >= nrow (scl$nodes))
            break
    }
    n <- which (table (scl$nodes$cluster) == 2)
    scl$nodes$cluster [scl$nodes$cluster %in% n] <- NA

    scl$nodes <- dplyr::bind_cols (scl$nodes, xy)

    return (scl)
}
