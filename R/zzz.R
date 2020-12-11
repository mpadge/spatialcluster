.onLoad <- function (libname, pkgname) { # nolint

    # make data set names global to avoid CHECK notes
    utils::globalVariables (".")
    utils::globalVariables ("i")
    utils::globalVariables ("x")
    utils::globalVariables ("y")
    utils::globalVariables ("d")
    utils::globalVariables ("id")
    utils::globalVariables ("v")
    utils::globalVariables ("V1")
    utils::globalVariables ("V2")
    utils::globalVariables ("V3")
    utils::globalVariables ("from")
    utils::globalVariables ("to")
    utils::globalVariables ("clnum")
    utils::globalVariables ("cluster")
    utils::globalVariables ("merged")
    utils::globalVariables ("comp")
    utils::globalVariables ("value")
    utils::globalVariables ("na.omit")
    utils::globalVariables ("xfr")
    utils::globalVariables ("xto")
    utils::globalVariables ("yfr")
    utils::globalVariables ("yto")
    utils::globalVariables ("ind")

    invisible ()
}
