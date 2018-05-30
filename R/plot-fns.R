#' scl_hulls
#'
#' Calculate convex hulls around redcap clusters, mostly cribbed from
#' osmplotr/R/add-osm-groups.R
#'
#' @param tree Spanning tree obtained from \link{scl_redcap}
#' @param xy Matrix of spatial coordinates of points indexed by \code{tree}.
#' @return tibble of (id, x, y), where the coordinates trace the convex hulls
#' for each cluster id
#' @noRd
scl_hulls <- function (nodes)
{
    ncl <- length (unique (nodes$cl [!is.na (nodes$cl)]))
    bdry <- list ()
    for (i in seq (ncl))
    {
        if (length (which (nodes$cl == i)) > 1)
        {
            xyi <- nodes %>%
                dplyr::filter (cl == i)
            xy2 <- spatstat::ppp (xyi$x, xyi$y,
                                  xrange = range (xyi$x),
                                  yrange = range (xyi$y))
            ch <- spatstat::convexhull (xy2)
            bdry [[i]] <- cbind (i, ch$bdry[[1]]$x, ch$bdry[[1]]$y)
        }
    }
    bdry <- data.frame (do.call (rbind, bdry))
    names (bdry) <- c ("id", "x", "y")
    return (bdry)
}

#' scl_ahulls
#'
#' Calculate alpha hulls around clusters via the \pkg{alphahull} package
#'
#' @param tree Spanning tree obtained from \link{scl_redcap}
#' @param xy Matrix of spatial coordinates of points indexed by \code{tree}.
#' @param alpha Parameter used to create alpha hulls
#' @return tibble of (id, x, y), where the coordinates trace the convex hulls
#' for each cluster id
#' @noRd
scl_ahulls <- function (nodes, alpha = 0.1)
{
    ncl <- length (unique (nodes$cl [!is.na (nodes$cl)]))
    bdry <- list ()
    for (i in seq (ncl))
    {
        if (length (which (nodes$cl == i)) > 2)
        {
            xyi <- nodes %>%
                dplyr::filter (cl == i) %>%
                dplyr::select (x, y)

            a <- alphahull::ashape (xyi, alpha = alpha)$edges %>%
                data.frame ()

            xy <- rbind (data.frame (ind = a$ind1, x = a$x1, y = a$y1),
                         data.frame (ind = a$ind2, x = a$x2, y = a$y2)) %>%
                    unique () %>%
                    dplyr::arrange (ind)
            inds <- data.frame (ind1 = a$ind1, ind2 = a$ind2)
            # Then just have to wrap those around xy:
            # TODO: Find a better way to do this!
            ind_seq <- as.numeric (inds [1, ])
            inds <- inds [-1, ]
            while (nrow (inds) > 0)
            {
                j <- which (inds$ind1 == tail (ind_seq, n = 1))
                if (length (j) > 0)
                {
                    ind_seq <- c (ind_seq, inds [j, 2])
                } else
                {
                    j <- which (inds$ind2 == tail (ind_seq, n = 1))
                    ind_seq <- c (ind_seq, inds [j, 1])
                }
                inds <- inds [-j, , drop = FALSE] #nolint
            }
            xy <- xy [match (ind_seq, xy$ind), ]
            bdry [[i]] <- cbind (i, xy$x, xy$y)
        }
    }
    bdry <- data.frame (do.call (rbind, bdry))
    names (bdry) <- c ("id", "x", "y")
    return (bdry)
}

#' plot.scl
#' @method plot scl
#' @param x object to be plotted
#' @param convex Should hull be convex? If not, the \code{ashape} routine from
#' the \pkg{alphahull} package is used to generate non-convex hulls, generated
#' with the \code{hull_alpha} parameter
#' @param hull_alpha alpha value of non-convex hulls (see ?alphashape::ashape
#' for details).
#' @param ... ignored here
#' @export
#' @examples
#' n <- 100
#' xy <- matrix (runif (2 * n), ncol = 2)
#' dmat <- matrix (runif (n ^ 2), ncol = n)
#' scl <- scl_redcap (xy, dmat, ncl = 4)
#' plot (scl)
#' # Connect clusters according to highest (\code{shortest = FALSE}) values of
#' # \code{dmat}:
#' scl <- scl_redcap (xy, dmat, ncl = 4, shortest = FALSE, full_order = FALSE)
#' plot (scl)
plot.scl <- function (x, ..., convex = TRUE, hull_alpha = 0.1)
{
    if (convex)
        hulls <- scl_hulls (x$nodes)
    else
        hulls <- scl_ahulls (x$nodes, alpha = hull_alpha)

    nc <- length (unique (x$nodes$cl [!is.na (x$nodes$cl)]))

    # clnum in cl_cols is + 1 because xy below increases cluster numbers by 1 to
    # allocate cl_num == 1 to unassigned points
    cl_cols <- rainbow (nc) %>%
        tibble::as.tibble () %>%
        dplyr::mutate (cl = seq (nc) + 1) %>%
        dplyr::rename (col = value)

    xy <- x$nodes %>%
        dplyr::mutate (cl = ifelse (is.na (cl), 1, cl + 1)) %>%
        dplyr::left_join (cl_cols, by = "cl") %>%
        dplyr::mutate (col = ifelse (is.na (col), "#333333FF", col))

    y <- id <- NULL # suppress no visible binding warnings
    hull_aes <- ggplot2::aes (x = x, y = y, group = id)
    hull_width <- 0.5
    g <- ggplot2::ggplot (xy, ggplot2::aes (x = x, y = y)) +
        ggplot2::geom_point (size = 5, color = xy$col,
                             show.legend = FALSE) +
        ggplot2::geom_polygon (data = hulls,
                               mapping = hull_aes,
                               colour = cl_cols$col [hulls$id],
                               fill = cl_cols$col [hulls$id],
                               alpha = 0.1,
                               size = hull_width) +
        ggthemes::theme_solarized ()

    print (g)
    invisible (g)
}
