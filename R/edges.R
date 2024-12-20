#' scl_edges_tri
#'
#' Generate triangulated nearest-neighbour edges between a set of input points
#'
#' @inheritParams scl_redcap
#' @examples
#' set.seed (1)
#' n <- 100
#' xy <- matrix (runif (2 * n), ncol = 2)
#' edges <- scl_edges_tri (xy)
#' @noRd
scl_edges_tri <- function (xy, shortest = TRUE) {

    if (!inherits (xy, "data.frame")) {
        xy <- as.data.frame (xy)
    }
    if (ncol (xy) > 2) {
        xy <- dplyr::select (xy, c (x, y))
    }
    names (xy) <- c ("x", "y")
    nbs <- tripack::tri.mesh (xy) |>
        tripack::neighbours ()

    n <- length (nbs)
    edges <- lapply (
        seq (n),
        function (i) cbind (i, nbs [[i]])
    )
    edges <- do.call (rbind, edges)

    edges <- tibble::tibble (
        from = edges [, 1],
        to = edges [, 2]
    )

    dxy <- as.matrix (stats::dist (xy))

    append_dist_to_edges (edges, dxy, shortest)
}

#' scl_edges_nn
#'
#' Generate distance-based nearest-neighbour edges between a set of input
#' points, ensuring that all edges connect to a single component. The minimal
#' spanning tree is constructed from **spatial** distances, not from distances
#' given in `dmat`.
#' @param nnbs Number of nearest neighbours
#' @inheritParams scl_redcap
#'
#' @return A `tibble` of `from` and `to` vertex indices for the minimal spanning
#' tree edges, along with corresponding spatial distances calculated from 'xy'.
#'
#' @examples
#' set.seed (1)
#' n <- 100
#' xy <- matrix (runif (2 * n), ncol = 2)
#' scl_edges_nn (xy, nnbs = 3L)
#' @noRd
scl_edges_nn <- function (xy, nnbs, shortest = TRUE) {

    # Initially contruct with nnbs + 1, because the `d` matrix includes
    # self-distances of zero, which are subsequently removed.
    nnbs <- nnbs + 1
    d <- apply (
        as.matrix (stats::dist (xy)),
        2,
        function (i) order (i, decreasing = !shortest) [seq_len (nnbs)]
    )

    edges <- tibble::tibble (
        from = rep (as.integer (colnames (d)),
            each = nnbs
        ),
        to = as.vector (d)
    )
    # rm self-edges:
    edges <- edges [which (edges$from != edges$to), ]

    # then ensure that the minimal spanning tree is included, to ensure all
    # nearest neighbour edges are connected in a single component. The distances
    # used for this MST are spatial distances, not from `dmat`.
    dxy <- as.matrix (stats::dist (xy))

    n <- nrow (xy)
    edges_all <- tibble::tibble (
        from = rep (seq_len (n), n),
        to = rep (seq_len (n), each = n),
        d = as.vector (dxy)
    ) |>
        dplyr::arrange (d, from, to) |>
        dplyr::filter (from != to)

    mst <- scl_spantree_ord1 (edges_all)
    # duplicate all of those:
    mst <- rbind (
        mst [, c ("from", "to")],
        tibble::tibble (from = mst$to, to = mst$from)
    )

    edges <- rbind (edges, mst)
    edges <- edges [which (!duplicated (edges)), ]

    # Then append final distances from `dxy` to the return value:
    edges <- append_dist_to_edges (edges, dxy, shortest)

    return (edges)
}

append_dist_to_edges <- function (edges, dmat, shortest) {
    index <- (edges$to - 1) * nrow (dmat) + edges$from
    edges$d <- dmat [index]

    if (shortest) {
        edges <- dplyr::arrange (edges, d)
    } # lowest-to-highest
    else {
        edges <- dplyr::arrange (edges, dplyr::desc (d))
    }

    return (edges)
}

#' scl_edges_all
#'
#' Generate full set of edges between a set of input points
#'
#' @param dmat Either a spatial distance matrix generated from 'xy', or a
#' separate distance matrix passed as the 'dmat' parameter to the main
#' \link{scl_redcap} function.
#' @inheritParams scl_redcap
#' @noRd
scl_edges_all <- function (xy, dmat, shortest = TRUE) {

    n <- nrow (dmat)
    edges <- tibble::tibble (
        from = rep (seq (n), times = n),
        to = rep (seq (n), each = n),
        d = as.vector (dmat)
    )
    edges <- na.omit (edges)

    if (shortest) {
        edges <- dplyr::arrange (edges, d)
    } # lowest-to-highest
    else {
        edges <- dplyr::arrange (edges, dplyr::desc (d))
    }

    return (edges)
}
