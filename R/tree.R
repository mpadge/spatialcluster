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

#' scl_spantree_alk
#'
#' Generate a spanning tree from full-order, average linkage clustering (ALK)
#' relationships expressed via a set of edges
#'
#' @param edges A set of edges resulting from \link{scl_edges}, which are sorted
#' in ascending order according to user-specified data. The only aspect of that
#' data which affect tree construction is this order, so only the set of
#' \code{edges} are needed here
#'
#' @return A tree
#' @noRd
scl_spantree_alk <- function (edges)
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
    list (tree_in = tree [ncl:n, ], tree_out = tree [1:(n - 1),])
}
