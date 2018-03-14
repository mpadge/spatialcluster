.onLoad <- function (libname, pkgname)
{
    # make data set names global to avoid CHECK notes
    utils::globalVariables (".")
    utils::globalVariables ("i")
    utils::globalVariables ("x")
    utils::globalVariables ("y")
    utils::globalVariables ("d")
    utils::globalVariables ("id")
    utils::globalVariables ("v")
    utils::globalVariables ("V1")
    utils::globalVariables ("from")
    utils::globalVariables ("to")
    utils::globalVariables ("clnum")
    utils::globalVariables ("comp")
    utils::globalVariables ("value")

    invisible ()
}
