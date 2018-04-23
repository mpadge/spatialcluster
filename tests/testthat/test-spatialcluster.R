context("tree")
test_that("scl structure", {
              n <- 100
              xy <- matrix (runif (2 * n), ncol = 2)
              dmat <- matrix (runif (n ^ 2), ncol = n)
              scl <- scl_cluster (xy, dmat, ncl = 4)
              expect_is (scl, "scl")
              expect_true (scl$pars$ncl >= 4)
              expect_true (all (names (scl) %in%
                                c ("xy", "tree", "tree_rest", "pars")))
              expect_true (nrow (scl$tree) < n)
})

test_that("scl methods", {
              n <- 100
              xy <- matrix (runif (2 * n), ncol = 2)
              dmat <- matrix (runif (n ^ 2), ncol = n)
              scl1 <- scl_cluster (xy, dmat, ncl = 4, shortest = FALSE)
              scl2 <- scl_cluster (xy, dmat, ncl = 4, shortest = TRUE)
              expect_true (!identical (scl1, scl2))

              # these are the default pars:
              scl3 <- scl_cluster (xy, dmat, ncl = 4, full_order = TRUE)
              expect_true (!identical (scl1, scl2))
              scl4 <- scl_cluster (xy, dmat, ncl = 4, linkage = "single")
              expect_true (identical (scl3, scl4))
              scl5 <- scl_cluster (xy, dmat, ncl = 4, linkage = "average")
              expect_false (identical (scl4, scl5))
              scl6 <- scl_cluster (xy, dmat, ncl = 4, linkage = "complete")
              expect_false (identical (scl6, scl5))
              expect_false (identical (scl6, scl4))
})

test_that("recluster", {
              n <- 100
              xy <- matrix (runif (2 * n), ncol = 2)
              dmat <- matrix (runif (n ^ 2), ncol = n)
              scl <- scl_cluster (xy, dmat, ncl = 4)
              scl2 <- scl_recluster (scl, ncl = 3)
              expect_message ( scl3 <- scl_cluster (scl, ncl = 3),
                              "scl_cluster is for initial cluster")
              expect_identical (scl2, scl3)
              expect_true (!identical (scl, scl2))
})

test_that("plot", {
              n <- 100
              xy <- matrix (runif (2 * n), ncol = 2)
              dmat <- matrix (runif (n ^ 2), ncol = n)
              scl <- scl_cluster (xy, dmat, ncl = 4)
              g <- plot (scl)
              expect_is (g, "ggplot")
              g2 <- plot (scl, tree = TRUE)
              expect_true (!identical (g, g2))
})
