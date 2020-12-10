#' scl_tbl
#'
#' Convert anything to a tibble
#'
#' @param xy A rectangular object containing the coordinates
#' @return A tibble-ified version of the input, with coordinate columns
#' identified and re-labelled "x" and "y"
#' @noRd
scl_tbl <- function (xy)
{
    if (!inherits (xy, "data.frame"))
    {
        if (!is.numeric (xy))
            stop ("coordinates must be numeric")
        if (is.vector (xy))
            stop ("coordinates require at least 2 columns")
        xy <- data.frame (xy)
    }
    xi <- grep ("^x|^lon", names (xy), ignore.case = TRUE)
    yi <- grep ("^y|^lat", names (xy), ignore.case = TRUE)
    if (length (xi) == 1 & length (yi) == 1)
    {
        names (xy) [xi] <- "x"
        names (xy) [yi] <- "y"
    } else if (ncol (xy) == 2)
    {
        colnames (xy) <- c ("x", "y")
    } else
    {
        stop ("Cannot determine unambiguous coordinate columns")
    }
    tibble::as_tibble (xy)
}

#' scl_linkage_type
#'
#' Convert \code{linkage} string arg to matching type
#' @param linkage Type of linkage
#' @return Strict match to one of three options
#' @noRd
scl_linkage_type <- function (linkage)
{
    linkages <- c ("single", "average", "complete", "full")
    i <- grep (linkage, linkages, ignore.case = TRUE)
    if (length (i) == 0)
        stop ("linkage must be one of (single, average, complete, full)")

    return (linkages [i])
}
