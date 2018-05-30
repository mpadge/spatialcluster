context("exact")

test_that("structure", {
              n <- 100
              xy <- matrix (runif (2 * n), ncol = 2)
              dmat <- matrix (runif (n ^ 2), ncol = n)
              ncl <- 4
              scl <- scl_exact (xy, dmat, ncl = ncl)
              expect_is (scl, "scl")
              expect_true (scl$pars$ncl == ncl)
              expect_true (all (names (scl) %in%
                                c ("merges", "ord", "nodes", "pars")))
              cl <- scl$nodes$cluster [!is.na (scl$nodes$cluster)]
              expect_true (length (unique (cl)) == ncl)
})

test_that("methods", {
              n <- 100
              xy <- matrix (runif (2 * n), ncol = 2)
              dmat <- matrix (runif (n ^ 2), ncol = n)
              ncl <- 8
              scl1 <- scl_exact (xy, dmat, ncl = ncl, linkage = "single")
              scl2 <- scl_exact (xy, dmat, ncl = ncl, linkage = "average")
              expect_true (!identical (scl1, scl2))
              cl1 <- scl1$nodes$cluster [!is.na (scl1$nodes$cluster)]
              expect_equal (length (unique (cl1)), ncl)
              cl2 <- scl2$nodes$cluster [!is.na (scl2$nodes$cluster)]
              expect_equal (length (unique (cl2)), ncl)
})
