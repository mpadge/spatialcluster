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
#' @examples
#' # xy can be matrix, data.frame, or tibble:
#' xy <- matrix (runif (100), ncol = 2) %>%
#'     tibble::as.tibble () %>%
#'     dplyr::rename (x = V1, y = V2)
#' edges <- scl_edges (xy)
#' # plot edges:
#' \dontrun{
#' plot (xy, pch = 19)
#' with (edges, segments (xy$x [from], xy$y [from], xy$x [to], xy$y [to]))
#' }
scl_edges <- function (xy)
{
    if (!(inherits (xy, "data.frame") | inherits (xy, "matrix")))
        stop ("xy must be a rectangular structure")

    xy <- as.data.frame (xy)
    # The order of x and y columns in irrelevant For triangulation, so:
    names (xy) <- c ("x", "y")
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
