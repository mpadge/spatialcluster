#' spatialcluster.
#'
#' @name spatialcluster
#' @docType package
#' @importFrom dplyr arrange desc filter left_join mutate rename
#' @importFrom ggplot2 aes ggplot geom_point geom_polygon
#' @importFrom ggthemes solarized_pal
#' @importFrom graphics plot
#' @importFrom grDevices rainbow
#' @importFrom magrittr %>% %<>%
#' @importFrom methods is
#' @importFrom spatstat convexhull ppp
#' @importFrom stats na.omit
#' @importFrom tibble tibble as.tibble
#' @importFrom tripack neighbours tri.mesh
#' @useDynLib spatialcluster, .registration = TRUE
NULL

#' Pipe operator
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL
