#' scl_hulls
#'
#' Calculate convex hulls around clusters, mostly cribbed from
#' osmplotr/R/add-osm-groups.R
#'
#' @param tree Spanning tree obtained from \link{scl_cluster}
#' @param xy Matrix of spatial coordinates of points indexed by \code{tree}.
#' @return tibble of (id, x, y), where the coordinates trace the convex hulls
#' for each cluster id
#' @noRd
scl_hulls <- function (tree, xy)
{
    xy <- as.matrix (xy)
    ncomp <- length (unique (tree$comp))
    bdry <- list ()
    for (i in seq (ncomp))
    {
        if (length (which (tree$comp == i)) > 1)
        {
            xyi <- tree %>%
                dplyr::filter (comp == i) %>%
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

#' scl_plot
#'
#' plot cluster groups
#' 
#' @param tree Spanning tree obtained from \link{scl_cluster}
#' @param xy Matrix of spatial coordinates of points indexed by \code{tree}.
#' @return (Invisible) \link{ggplot2} object containing the plot
#' @export
#' @examples
#' \dontrun{
#' xy <- matrix (runif (100), ncol = 2)
#' edges <- scl_edges (xy)
#' # add some fake data to the edges
#' edges %<>% dplyr::mutate (d = runif (nrow (.))) %>%
#'    dplyr::arrange (desc (d))
#' # get tree with component numbers
#' ncl <- 12 # desired number of clusters/components
#' tree <- scl_spantree (edges) %>%
#'     scl_cuttree (edges, ncl = ncl)
#' xy <- tibble::tibble (x = xy [, 1], y = xy [, 2])
#' g <- scl_plot (tree, xy)
#' }
scl_plot <- function (tree, xy)
{
    hulls <- scl_hulls (tree, xy)
    nc <- length (unique (tree$comp))

    if (is (xy, "tbl_df"))
        xy <- tibble::as.tibble (xy)
    if (ncol (xy) == 2)
        xy %<>% dplyr::mutate (v = seq (nrow (xy)))
    names (xy) [1:2] <- c ("x", "y")

    ggthemes::solarized_pal (accent = "blue") (8)
    cl_cols <- rainbow (nc) %>%
        tibble::as.tibble () %>%
        dplyr::mutate (comp = seq (nc)) %>%
        dplyr::rename (col = value)

    tree <- dplyr::left_join (tree, cl_cols, by = "comp")
    edge2vert <- dplyr::bind_rows (dplyr::select (tree, c (from, comp)) %>%
                                       dplyr::rename (v = from),
                                   dplyr::select (tree, c (to, comp)) %>%
                                       dplyr::rename (v = to)) %>%
                dplyr::arrange (v) %>%
                unique ()
    xy %<>% dplyr::mutate (v = seq (nrow (xy))) %>%
        dplyr::left_join (edge2vert, by = "v") %>%
        dplyr::mutate (comp = ifelse (is.na (comp), 1, comp + 1)) %>%
        dplyr::left_join (cl_cols, by = "comp") %>%
        dplyr::mutate (col = ifelse (is.na (col), "#222222", col))

    hull_aes <- ggplot2::aes (x = x, y = y, group = id)
    hull_width <- 0.5
    g <- ggplot2::ggplot (xy, ggplot2::aes (x = x,
                                            y = y,
                                            colour = col)) +
        ggplot2::geom_point (size = 5, show.legend = FALSE) +
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
