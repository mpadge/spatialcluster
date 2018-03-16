#' scl_edges
#'
#' Generate triangulated edges between a set of input points
#'
#' @inheritParams scl_cluster
#' @noRd
scl_edges <- function (xy, dmat, shortest = TRUE)
{
    nbs <- dplyr::select (xy, c (x, y)) %>%
        tripack::tri.mesh () %>%
        tripack::neighbours ()

    n <- length (nbs)
    edges <- lapply (seq (n), function (i)
                     cbind (i, nbs [[i]])) %>%
            do.call (rbind, .) %>%
            as.tibble () %>%
            rename (from = i, to = V1)

    index <- (edges$to - 1) * nrow (dmat) + edges$from
    edges$d <- dmat [index]

    if (shortest)
        edges %<>% dplyr::arrange (d)
    else
        edges %<>% dplyr::arrange (dplyr::desc (d))

    return (edges)
}
