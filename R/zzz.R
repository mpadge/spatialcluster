.onLoad <- function (libname, pkgname)
{
    # make data set names global to avoid CHECK notes
    utils::globalVariables (".")
    utils::globalVariables ("i")
    utils::globalVariables ("V1")
    utils::globalVariables ("from")
    utils::globalVariables ("to")
    utils::globalVariables ("clnum")

    invisible ()
}
