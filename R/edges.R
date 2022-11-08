#' scl_edges_tri
#'
#' Generate triangulated nearest-neighbour edges between a set of input points
#'
#' @inheritParams scl_redcap
#' @noRd
scl_edges_tri <- function (xy, dmat, shortest = TRUE) {
    nbs <- dplyr::select (xy, c (x, y)) %>%
        tripack::tri.mesh () %>%
        tripack::neighbours ()

    n <- length (nbs)
    edges <- lapply (seq (n), function (i)
                     cbind (i, nbs [[i]])) %>%
            do.call (rbind, .)
    edges <- tibble::tibble (from = edges [, 1],
                             to = edges [, 2])

    append_dist_to_edges (edges, dmat, shortest)
}

#' scl_edges_nn
#'
#' Generate distance-based nearest-neighbour edges between a set of input
#' points, ensuring that all edges connect to a single component.
#' @param edges_all `data.frame` of all edges returned from
#' \link{scl_edges_all}.
#' @param nnbs Number of nearest neighbours
#'
#' @inheritParams scl_redcap
#' @noRd
scl_edges_nn <- function (xy, dmat, edges_all, nnbs, shortest = TRUE) {

    # Initially contruct with nnbs + 1, because the `d` matrix includes
    # self-distances of zero, which are subsequently removed.
    nnbs <- nnbs + 1
    d <- apply (as.matrix (stats::dist (xy)), 2, function (i)
                order (i, decreasing = !shortest) [seq_len (nnbs)])

    edges <- tibble::tibble (from = rep (as.integer (colnames (d)),
                                         each = nnbs),
                             to = as.vector (d))
    # rm self-edges:
    edges <- edges [which (edges$from != edges$to), ]

    # then ensure that the minimal spanning tree is included, to ensure all
    # nearest neighbour edges are connected in a single component.
    mst <- scl_spantree_ord1 (edges_all) [, c ("from", "to")]
    # duplicate all of those:
    mst <- rbind (mst, tibble::tibble (from = mst$to, to = mst$from))

    edges <- rbind (edges, mst)
    edges <- edges [which (!duplicated (edges)), ]

    edges <- append_dist_to_edges (edges, dmat, shortest)

    return (edges)
}

append_dist_to_edges <- function (edges, dmat, shortest) {
    index <- (edges$to - 1) * nrow (dmat) + edges$from
    edges$d <- dmat [index]

    if (shortest)
        edges %<>% dplyr::arrange (d) # lowest-to-highest
    else
        edges %<>% dplyr::arrange (dplyr::desc (d))

    return (edges)
}

#' scl_edges_all
#'
#' Generate full set of edges between a set of input points
#'
#' @inheritParams scl_redcap
#' @noRd
scl_edges_all <- function (xy, dmat, shortest = TRUE) {
    n <- nrow (dmat)
    edges <- tibble::tibble (from = rep (seq (n), times = n),
                             to = rep (seq (n), each = n),
                             d = as.vector (dmat))
    edges <- na.omit (edges)

    if (shortest)
        edges %<>% dplyr::arrange (d) # lowest-to-highest
    else
        edges %<>% dplyr::arrange (dplyr::desc (d))

    return (edges)
}
