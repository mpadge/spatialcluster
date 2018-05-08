context("exact")

test_that("structure", {
              n <- 100
              xy <- matrix (runif (2 * n), ncol = 2)
              dmat <- matrix (runif (n ^ 2), ncol = n)
              ncl <- 4
              scl <- scl_exact (xy, dmat, ncl = ncl)
              expect_is (scl, "scl_exact")
              expect_true (scl$ncl == ncl)
              expect_true (all (names (scl) %in%
                                c ("merges", "ord", "nodes", "ncl")))
              expect_true (length (unique (scl$nodes$cluster)) == ncl)
})
