#' scl_hulls
#'
#' Calculate convex hulls around clusters, mostly cribbed from
#' osmplotr/R/add-osm-groups.R
#'
#' @param tree Spanning tree obtained from \link{scl_redcap}
#' @param xy Matrix of spatial coordinates of points indexed by \code{tree}.
#' @return tibble of (id, x, y), where the coordinates trace the convex hulls
#' for each cluster id
#' @noRd
scl_hulls <- function (tree, xy)
{
    xy <- as.matrix (xy)
    ncl <- length (unique (tree$clnum))
    bdry <- list ()
    for (i in seq (ncl))
    {
        if (length (which (tree$clnum == i)) > 1)
        {
            xyi <- tree %>%
                dplyr::filter (clnum == i) %>%
                dplyr::select (from, to) %>%
                unlist () %>%
                unique () %>%
                sort () %>%
                xy [., ]
            xy2 <- spatstat::ppp (xyi [, 1], xyi [, 2],
                                  xrange = range (xyi [, 1]),
                                  yrange = range (xyi [, 2]))
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
scl_ahulls <- function (tree, xy, alpha = 0.1)
{
    xymat <- as.matrix (xy)
    ncl <- length (unique (tree$clnum))
    bdry <- list ()
    for (i in seq (ncl))
    {
        if (length (which (tree$clnum == i)) > 2)
        {
            xyi <- tree %>%
                dplyr::filter (clnum == i) %>%
                dplyr::select (from, to) %>%
                unlist () %>%
                unique () %>%
                sort () %>%
                xymat [., ]

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
#' # \coce{dmat}:
#' scl <- scl_redcap (xy, dmat, ncl = 4, shortest = FALSE, full_order = FALSE)
#' plot (scl)
plot.scl <- function (x, ..., convex = TRUE, hull_alpha = 0.1)
{
    if (convex)
        hulls <- scl_hulls (x$tree, x$xy)
    else
        hulls <- scl_ahulls (x$tree, x$xy, alpha = hull_alpha)

    nc <- length (unique (x$tree$clnum))

    # clnum in cl_cols is + 1 because xy below increases cluster numbers by 1 to
    # allocate cl_num == 1 to unassigned points
    cl_cols <- rainbow (nc) %>%
        tibble::as.tibble () %>%
        dplyr::mutate (clnum = seq (nc) + 1) %>%
        dplyr::rename (col = value)

    edge2vert <- dplyr::bind_rows (dplyr::select (x$tree, c (from, clnum)) %>%
                                       dplyr::rename (v = from),
                                   dplyr::select (x$tree, c (to, clnum)) %>%
                                       dplyr::rename (v = to)) %>%
                dplyr::arrange (v) %>%
                unique ()
    xy <- x$xy
    xy %<>% dplyr::mutate (v = seq (nrow (xy))) %>%
        dplyr::left_join (edge2vert, by = "v") %>%
        dplyr::mutate (clnum = ifelse (is.na (clnum), 1, clnum + 1)) %>%
        dplyr::left_join (cl_cols, by = "clnum") %>%
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
