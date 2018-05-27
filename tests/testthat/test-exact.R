context("exact")

test_that("structure", {
              n <- 100
              xy <- matrix (runif (2 * n), ncol = 2)
              dmat <- matrix (runif (n ^ 2), ncol = n)
              ncl <- 4
              scl <- scl_exact (xy, dmat, ncl = ncl)
              expect_is (scl, "scl_exact")
              expect_true (scl$pars$ncl == ncl)
              expect_true (all (names (scl) %in%
                                c ("merges", "ord", "nodes", "pars")))
              expect_true (length (unique (scl$nodes$cluster)) == ncl)
})

test_that("methods", {
              n <- 100
              xy <- matrix (runif (2 * n), ncol = 2)
              dmat <- matrix (runif (n ^ 2), ncol = n)
              ncl <- 8
              scl1 <- scl_exact (xy, dmat, ncl = ncl, method = "single")
              scl2 <- scl_exact (xy, dmat, ncl = ncl, method = "average")
              expect_true (!identical (scl1, scl2))
              expect_equal (length (unique (scl1$nodes$cluster)), ncl)
              expect_equal (length (unique (scl2$nodes$cluster)), ncl)
})
