#' scl_edges
#'
#' Generate triangulated edges between a set of input points
#'
#' @param xy A matrix, data.frame, or other rectangular structure containing x
#' and y coordinates for triangulation
#'
#' @return A \pkg{tibble} of two columns indexing into \code{xy}, where each row
#' denotes the start and end of an edge
#'
#' @export
scl_edges <- function (xy)
{
    nbs <- tripack::tri.mesh (xy) %>%
        tripack::neighbours ()

    n <- length (nbs)
    edges <- lapply (seq (n), function (i)
                     cbind (i, nbs [[i]])) %>%
            do.call (rbind, .) %>%
            as.tibble () %>%
            rename (from = i, to = V1)

    return (edges)
}
