#' scl_spantree_O1
#'
#' Generate a spanning tree from first-order relationships expressed via a set
#' of edges
#'
#' @param edges A set of edges resulting from \link{scl_edges}, which are sorted
#' in ascending order according to user-specified data. The only aspect of that
#' data which affect tree construction is this order, so only the set of
#' \code{edges} are needed here
#'
#' @return A tree
#' @noRd
scl_spantree_O1 <- function (edges)
{
    n <- edges %>%
        dplyr::select (from, to) %>%
        unlist () %>%
        max ()
    tree <- tibble::tibble (from = integer(), to = integer())
    clusters <- tibble::tibble (id = seq (n), clnum = seq (n))

    # make the minimal tree:
    for (e in seq (nrow (edges)))
    {
        clf <- clusters$clnum [edges$from [e]]
        clt <- clusters$clnum [edges$to [e]]
        if (clf != clt)
        {
            cli <- min (c (clf, clt))
            clj <- max (c (clf, clt))
            clusters %<>% dplyr::mutate (clnum = replace (clnum,
                                                          clnum == clj,
                                                          cli))
            tree %<>% dplyr::bind_rows (tibble::tibble (from = edges$from [e],
                                                        to = edges$to [e]))
        }
        if (length (unique (clusters$clnum)) == 1)
            break
    }
    return (tree)
}

#' scl_spantree_slk
#'
#' Generate a spanning tree from full-order, single linkage clustering (ALK)
#' relationships expressed via a set of edges
#'
#' @param edges_all A set of ALL edges resulting from \link{scl_edges_all},
#' which are sorted in ascending order according to user-specified data.
#' @param edges_nn A equivalent set of nearest neighbour edges only, resulting
#' from \link{scl_edges_nn}.
#'
#' @return A tree
#' @noRd
scl_spantree_slk <- function (edges_all, edges_nn)
{
    clusters <- rcpp_slk (edges_all, edges_nn) + 1
    tibble::tibble (from = edges_nn$from [clusters],
                    to = edges_nn$to [clusters])
}

#' scl_cuttree
#'
#' Cut a tree generated with \link{scl_spantree} into a specified number of
#' clusters or components
#'
#' @param tree result of \link{scl_spantree}
#' @param edges A set of edges resulting from \link{scl_edges}, but with
#' additional data specifying edge weights, distances, or desired properties
#' from which to construct the tree
#' @param ncl Number of clusters or components into which tree is to be cut
#'
#' @return Two trees, \code{tree_in} containing the desired cut portion, and
#' \code{tree_out} the excised portion, retained here to enable later
#' reconstruction of full-tree for subsequent re-cutting.
#' @noRd
scl_cuttree <- function (tree, edges, ncl)
{
    tree %<>% left_join (edges, by = c ("from", "to"))
    n <- nrow (tree)

    # define component as > 2 members, and cut tree until that is attained, or
    # until tree is only calculated from < half the points
    ncomps <- 1
    ncli <- ncl - 1
    ncmax <- 0
    while (ncomps < ncl & ncli < (n / 2))
    {
        ncli <- ncli + 1
        cmp <- tree_components (tree [ncli:nrow (tree), ])$comp
        ncomps <- length (which (table (cmp) > 2))

        if (ncomps > ncmax)
        {
            ncmax <- ncomps
            ncl_max <- ncli
        }
    }
    if (ncli >= (n / 2))
    {
        message ("Only able to cut tree into maximum of ", ncmax,
                 " components")
        ncli <- ncl_max
    }

    list (tree_in = tree [ncli:n, ], tree_out = tree [1:(ncli - 1),],
          ncl = ncli)
}
