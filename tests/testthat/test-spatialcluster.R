context("tree")
test_that("cluster", {
              n <- 20
              xy <- matrix (runif (2 * n), ncol = 2)
              dmat <- matrix (runif (n ^ 2), ncol = n)
              scl <- scl_cluster (xy, dmat, ncl = 4)
              expect_is (scl, "scl")
              expect_equal (scl$ncl, 4)
              expect_true (all (names (scl) %in%
                                c ("xy", "tree", "tree_rest", "ncl")))
              expect_true (nrow (scl$tree) < n)
              scl2 <- scl_cluster (xy, dmat, ncl = 4, shortest = FALSE)
              expect_true (!identical (scl, scl2))
})

test_that("recluster", {
              n <- 20
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
              n <- 20
              xy <- matrix (runif (2 * n), ncol = 2)
              dmat <- matrix (runif (n ^ 2), ncol = n)
              scl <- scl_cluster (xy, dmat, ncl = 4)
              g <- plot (scl)
              expect_is (g, "ggplot")
})
