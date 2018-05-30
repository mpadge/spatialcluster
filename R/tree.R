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
    clusters <- tibble::tibble (id = seq (n), cluster = seq (n))

    # make the minimal tree:
    for (e in seq (nrow (edges)))
    {
        clf <- clusters$cluster [edges$from [e]]
        clt <- clusters$cluster [edges$to [e]]
        if (clf != clt)
        {
            cli <- min (c (clf, clt))
            clj <- max (c (clf, clt))
            clusters %<>% dplyr::mutate (cluster = replace (cluster,
                                                            cluster == clj,
                                                            cli))
            tree %<>% dplyr::bind_rows (tibble::tibble (from = edges$from [e],
                                                        to = edges$to [e]))
        }
        if (length (unique (clusters$cluster)) == 1)
            break
    }
    return (tree)
}

#' scl_spantree_slk
#'
#' Generate a spanning tree from full-order, single linkage clustering (SLK)
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

#' scl_spantree_alk
#'
#' Generate a spanning tree from full-order, average linkage clustering (ALK)
#' relationships expressed via a set of edges
#'
#' @inheritParams scl_spantree_slk
#' @noRd
scl_spantree_alk <- function (edges)
{
    clusters <- rcpp_alk (edges) + 1
    tibble::tibble (from = edges$from [clusters],
                    to = edges$to [clusters])
}

#' scl_spantree_clk
#'
#' Generate a spanning tree from full-order, complete linkage clustering (CLK)
#' relationships expressed via a set of edges
#'
#' @inheritParams scl_spantree_slk
#' @noRd
scl_spantree_clk <- function (edges_all, edges_nn)
{
    clusters <- rcpp_clk (edges_all, edges_nn) + 1
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
#' @return Modified version of \code{tree}, including an additional column
#' specifying the cluster number of each edge, with NA for edges that lie
#' between clusters.
#'
#' @note The \code{rcpp_cut_tree} routine in \code{src/cuttree} includes
#' \code{constexpr MIN_CLUSTER_SIZE = 3}.
#'
#' @noRd
scl_cuttree <- function (tree, edges, ncl)
{
    tree %<>%
        dplyr::left_join (edges, by = c ("from", "to")) %>%
        dplyr::mutate (cluster = rcpp_cut_tree (., ncl = ncl) + 1)

    return (tree)
}
