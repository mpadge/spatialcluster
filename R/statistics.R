#' scl_statistics
#'
#' @param scl Output of either \link{scl_redcap} or \link{scl_full}.
#'
#' @return A modified version of the input object with statistics appended.
#'
#' @noRd
scl_statistics <- function (scl) {

    tree <- scl$tree |>
        dplyr::mutate (tf = paste0 (to, "-", from))
    edges_in <- scl$tree [which (scl$tree$cluster >= 0), ] |>
        dplyr::mutate (tf = paste0 (to, "-", from))
    tree <- tree [which (!tree$tf %in% edges_in$tf), ]

    # t.test (edges_in$d, tree$d, alternative = "greater", var.equal = TRUE)
    tt_global <- stats::t.test (
        edges_in$d,
        tree$d,
        alternative = "less",
        var.equal = TRUE
    )
    tt_global <- c (tt_global$statistic, tt_global$parameter, tt_global$p.value)
    names (tt_global) <- c ("statistic", "parameter", "p.value")

    tt_cl <- vapply (
        sort (unique (scl$tree$cluster)),
        function (i) {
            index <- which (edges_in$cluster == i)
            if (length (index) <= 3) {
                res <- rep (NA, 3L)
            } else {
                tt <- stats::t.test (edges_in$d [index], tree$d,
                    alternative = "less",
                    var.equal = TRUE
                )
                res <- c (tt$statistic, tt$parameter, tt$p.value)
            }
            names (res) <- c ("statistic", "parameter", "p.value")
            return (res)
        },
        numeric (3)
    ) |>
        t ()

    scl$statistics <- list (tt_global = tt_global, tt_clusters = tt_cl)

    return (scl)
}
