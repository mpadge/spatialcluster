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
#' @examples
#' xy <- matrix (runif (100), ncol = 2)
#' edges <- scl_edges (xy)
#' # add some fake data to the edges
#' edges %<>% dplyr::mutate (d = runif (nrow (.))) %>%
#'    dplyr::arrange (desc (d))
#' tree <- scl_spantree (edges)
#' \dontrun{
#' # plot the tree
#' plot (xy, pch = 19)
#' with (edges, segments (xy [from, 1], xy [from, 2], xy [to, 1], xy [to, 2],
#'                        col = "gray", lty = 2))
#' with (tree, segments (xy [from, 1], xy [from, 2], xy [to, 1], xy [to, 2],
#'                       col = "red", lwd = 4))
#' }
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
#' @return A tree
#' @export
#' @examples
#' xy <- matrix (runif (100), ncol = 2)
#' edges <- scl_edges (xy)
#' # add some fake data to the edges
#' edges %<>% dplyr::mutate (d = runif (nrow (.))) %>%
#'    dplyr::arrange (desc (d))
#' tree <- scl_spantree (edges)
#' ncl <- 8 # desired number of clusters/components
#' tree <- scl_cuttree (tree, edges, ncl = ncl)
scl_cuttree <- function (tree, edges, ncl)
{
    tree %<>% left_join (edges, by = c ("from", "to"))
    n <- nrow (tree)
    tree <- tree [ncl:n, ]
    scl_components (tree)
}

#' scl_components
#'
#' Get component vector of tree edges
#'
#' @param tree result of \link{scl_spantree}
#'
#' @return Modified version of \code{tree}, with additional \code{comp} column
#' enumerating the component numbers
#' @export
#' @examples
#' xy <- matrix (runif (100), ncol = 2)
#' edges <- scl_edges (xy)
#' # add some fake data to the edges
#' edges %<>% dplyr::mutate (d = runif (nrow (.))) %>%
#'    dplyr::arrange (desc (d))
#' tree <- scl_spantree (edges) # plain tree; no components
#' ncl <- 8 # desired number of clusters/components
#' tree <- scl_cuttree (tree, edges, ncl = ncl) # tree with component numbers
scl_components <- function (tree)
{
    tree$id <- seq (nrow (tree))
    cmps <- rcpp_get_component_vector (tree)
    tibble::tibble (id = as.numeric (cmps$edge_id),
                    comp = cmps$edge_component) %>%
        dplyr::arrange (id) %>%
        dplyr::left_join (tree, ., by = "id") %>%
        dplyr::select (from, to, comp)
}
