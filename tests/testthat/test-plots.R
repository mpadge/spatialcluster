context("plot")

test_that("redcap plot", {
              n <- 100
              xy <- matrix (runif (2 * n), ncol = 2)
              dmat <- matrix (runif (n ^ 2), ncol = n)
              scl <- scl_redcap (xy, dmat, ncl = 4)
              g1 <- plot (scl)
              expect_is (g1, "ggplot")
              g2 <- plot (scl, convex = FALSE)
              # The hulls are then contained in
              h1 <- g1$layers [[2]]$data
              h2 <- g2$layers [[2]]$data
              expect_true (!identical (h1, h2))
})
