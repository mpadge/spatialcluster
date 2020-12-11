context("plot")

test_that("redcap plot", {
              set.seed (1)
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

test_that("full plot", {
              set.seed (1)
              n <- 100
              xy <- matrix (runif (2 * n), ncol = 2)
              dmat <- matrix (runif (n ^ 2), ncol = n)
              scl <- scl_full (xy, dmat, ncl = 4)
              g1 <- plot (scl)
              expect_is (g1, "ggplot")
              g2 <- plot (scl, convex = FALSE)
              # The hulls are then contained in
              h1 <- g1$layers [[2]]$data
              h2 <- g2$layers [[2]]$data
              expect_true (!identical (h1, h2))
})

test_that("plot_merges", {
              set.seed (1)
              n <- 100
              xy <- matrix (runif (2 * n), ncol = 2)
              dmat <- matrix (runif (n ^ 2), ncol = n)
              scl <- scl_full (xy, dmat, ncl = 8)
              expect_silent (plot_merges (scl))
              graphics.off ()
              expect_silent (plot_merges (scl, root_tree = TRUE))
              graphics.off ()
              scl <- scl_redcap (xy, dmat, ncl = 8)
              expect_error (plot_merges (scl),
                            "plot_merges can only be applied to scl objects")
})
