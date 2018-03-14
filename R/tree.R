#' scl_spantree
#'
#' Generate a spanning tree from a set of edges
#'
#' @param edges A set of edges resulting from \link{scl_edges}, but with
#' additional data specifying edge weights, distances, or desired properties
#' from which to construct the tree
#'
#' @return A tree
#' @export
scl_spantree <- function (edges)
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
