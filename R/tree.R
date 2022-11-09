#' scl_spantree_ord1
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
scl_spantree_ord1 <- function (edges) {

    tree <- rcpp_mst (edges) %>%
        dplyr::arrange (from, to) %>%
        tibble::tibble ()

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
#' from \link{scl_edges_tri} or \link{scl_edges_nn}.
#'
#' @return A tree
#' @noRd
scl_spantree_slk <- function (edges_all, edges_nn, shortest, quiet = FALSE) {

    clusters <- rcpp_slk (edges_all, edges_nn,
        shortest = shortest, quiet = quiet) + 1

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
scl_spantree_alk <- function (edges, shortest, quiet = FALSE) {

    clusters <- rcpp_alk (edges, shortest = shortest, quiet = quiet) + 1
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
scl_spantree_clk <- function (edges_all, edges_nn, shortest, quiet = FALSE) {

    clusters <- rcpp_clk (edges_all, edges_nn,
        shortest = shortest, quiet = quiet) + 1

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
#' @inheritParams scl_redcap
#'
#' @return Modified version of \code{tree}, including an additional column
#' specifying the cluster number of each edge, with NA for edges that lie
#' between clusters.
#'
#' @note The \code{rcpp_cut_tree} routine in \code{src/cuttree} includes
#' \code{constexpr MIN_CLUSTER_SIZE = 3}.
#'
#' @noRd
scl_cuttree <- function (tree, edges, ncl, shortest,
                         iterate_ncl = FALSE, quiet = FALSE) {

    num_clusters <- 0
    ncl_trial <- ncl

    quiet <- !(!quiet & nrow (tree) > 100)

    while (num_clusters < ncl) {

        if (num_clusters > 0 && !quiet) {
            message ("Not enough clusters found; re-starting search.")
        }

        tree_temp <- tree %>%
            dplyr::left_join (edges, by = c ("from", "to")) %>%
            dplyr::mutate (cluster = rcpp_cut_tree (
                .,
                ncl = ncl_trial,
                shortest = shortest,
                quiet = quiet
            ) + 1)
        num_clusters <- length (which (table (tree_temp$cluster) > 2))
        if (!quiet) {
            message ("Total clusters found with > 2 members: ", num_clusters)
        }
        ncl_trial <- ncl_trial + 1
        if (ncl_trial >= nrow (tree) || iterate_ncl)
            break
    }

    return (tree_temp)
}
