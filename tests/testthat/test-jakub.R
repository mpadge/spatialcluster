context("issue 22")

test_that("doesthisfail", {
    set.seed(2020-12-14)
    n <- 1000
    xy <- matrix (runif (2 * n), ncol = 2)
    dmat <- matrix (runif (n ^ 2), ncol = n)

    # works
    #scl <- scl_redcap (xy, dmat, ncl = 8, linkage = "single")
    #expect_is (scl, "scl")

    # R crashes (both inside and outside of RStudio)
    #scl2 <- scl_redcap (xy, dmat, ncl = 50, linkage = "single")
    #expect_is (scl, "scl")
})
