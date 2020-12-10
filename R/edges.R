#' scl_edges_tri
#'
#' Generate triangulated nearest-neighbour edges between a set of input points
#'
#' @inheritParams scl_redcap
#' @noRd
scl_edges_tri <- function (xy, dmat, shortest = TRUE)
{
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
#' Generate distance-based nearest-neighbour edges between a set of input points
#' @param nnbs Number of nearest neighbours
#'
#' @inheritParams scl_redcap
#' @noRd
scl_edges_nn <- function (xy, dmat, nnbs, shortest = TRUE)
{
    d <- apply (as.matrix (dist (xy)), 2, function (i)
                order (i, decreasing = !shortest) [1:nnbs])
    edges <- tibble::tibble (from = rep (as.integer (colnames (d)),
                                         each = nnbs),
                             to = as.vector (d))

    append_dist_to_edges (edges, dmat, shortest)
}

append_dist_to_edges <- function (edges, dmat, shortest)
{
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
scl_edges_all <- function (xy, dmat, shortest = TRUE)
{
    n <- nrow (dmat)
    edges <- cbind (seq (n), rep (seq (n), each = n), as.vector (dmat)) %>%
        tibble::as_tibble () %>%
        dplyr::rename (from = V1, to = V2, d = V3) %>%
        na.omit ()

    if (shortest)
        edges %<>% dplyr::arrange (d) # lowest-to-highest
    else
        edges %<>% dplyr::arrange (dplyr::desc (d))

    return (edges)
}
