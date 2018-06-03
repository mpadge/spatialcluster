#' scl_edges_nn
#'
#' Generate triangulated nearest-neighbour edges between a set of input points
#'
#' @inheritParams scl_redcap
#' @noRd
scl_edges_nn <- function (xy, dmat, distances = TRUE)
{
    nbs <- dplyr::select (xy, c (x, y)) %>%
        tripack::tri.mesh () %>%
        tripack::neighbours ()

    n <- length (nbs)
    edges <- lapply (seq (n), function (i)
                     cbind (i, nbs [[i]])) %>%
            do.call (rbind, .) %>%
            tibble::as.tibble () %>%
            dplyr::rename (from = i, to = V1)

    index <- (edges$to - 1) * nrow (dmat) + edges$from
    edges$d <- dmat [index]

    if (distances)
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
scl_edges_all <- function (xy, dmat, distances = TRUE)
{
    n <- nrow (dmat)
    edges <- cbind (seq (n), rep (seq (n), each = n), as.vector (dmat)) %>%
        tibble::as.tibble () %>%
        dplyr::rename (from = V1, to = V2, d = V3) %>%
        na.omit ()

    if (distances)
        edges %<>% dplyr::arrange (d) # lowest-to-highest
    else
        edges %<>% dplyr::arrange (dplyr::desc (d))

    return (edges)
}
