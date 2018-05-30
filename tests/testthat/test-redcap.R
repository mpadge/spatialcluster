context("redcap ")

test_that("structure", {
              n <- 100
              xy <- matrix (runif (2 * n), ncol = 2)
              dmat <- matrix (runif (n ^ 2), ncol = n)
              scl <- scl_redcap (xy, dmat, ncl = 4)
              expect_is (scl, "scl")
              expect_true (scl$pars$ncl >= 4)
              expect_true (all (names (scl) %in%
                                c ("xy", "tree", "nodes", "pars")))
              expect_true (nrow (scl$tree) < n)
})

test_that("methods", {
              n <- 100
              xy <- matrix (runif (2 * n), ncol = 2)
              dmat <- matrix (runif (n ^ 2), ncol = n)
              scl1 <- scl_redcap (xy, dmat, ncl = 4, shortest = FALSE)
              scl2 <- scl_redcap (xy, dmat, ncl = 4, shortest = TRUE)
              expect_true (!identical (scl1, scl2))

              # these are the default pars:
              scl3 <- scl_redcap (xy, dmat, ncl = 4, full_order = TRUE)
              expect_true (!identical (scl1, scl2))

              scl4 <- scl_redcap (xy, dmat, ncl = 4, linkage = "single")
              expect_true (identical (scl3, scl4))
              scl5 <- scl_redcap (xy, dmat, ncl = 4, linkage = "average")
              expect_false (identical (scl4, scl5))
              scl6 <- scl_redcap (xy, dmat, ncl = 4, linkage = "complete")
              expect_false (identical (scl6, scl5))
              expect_false (identical (scl6, scl4))
              scl7 <- scl_redcap (xy, dmat, ncl = 4, full_order = FALSE)
              expect_false (identical (scl7, scl4))
              expect_false (identical (scl7, scl5))
              expect_false (identical (scl7, scl6))
              expect_error (scl8 <- scl_redcap (xy, dmat, ncl = 4,
                                                linkage = "blah"),
                            "linkage must be one of")
})

test_that("recluster", {
              n <- 100
              xy <- matrix (runif (2 * n), ncol = 2)
              dmat <- matrix (runif (n ^ 2), ncol = n)
              scl <- scl_redcap (xy, dmat, ncl = 4)
              scl2 <- scl_recluster (scl, ncl = 3)
              expect_message ( scl3 <- scl_redcap (scl, ncl = 3),
                              "scl_redcap is for initial cluster")
              expect_identical (scl2, scl3)
              expect_true (!identical (scl, scl2))
})
