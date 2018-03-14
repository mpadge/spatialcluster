context("tree")
test_that("edges", {
              xy <- matrix (runif (100), ncol = 2)
              edges1 <- scl_edges (xy)
              expect_true (nrow (edges1) > nrow (xy))

              xy <- data.frame (xy)
              edges2 <- scl_edges (xy)
              expect_identical (edges1, edges2)

              xy <- tibble::as.tibble (xy)
              edges3 <- scl_edges (xy)
              expect_identical (edges1, edges3)

              names (xy) <- c ("y", "x") # reverse column names
              edges4 <- scl_edges (xy)
              expect_identical (edges1, edges4)
})

test_that("tree", {
              xy <- matrix (runif (100), ncol = 2)
              edges <- scl_edges (xy)
              edges %<>% dplyr::mutate (d = runif (nrow (.))) %>%
                  dplyr::arrange (desc (d))
              tree <- scl_spantree (edges)
              expect_true (nrow (tree) < nrow (edges))

              # tree with components
              treec <- scl_cuttree (tree, edges, ncl = 8)
              expect_identical (names (tree), c ("from", "to"))
              expect_identical (names (treec), c ("from", "to", "comp"))
              expect_true (length (unique (treec$comp)) <= 8) # can be less
              expect_true (max (treec$comp) <= 8)
})

test_that("plot", {
              xy <- matrix (runif (100), ncol = 2)
              edges <- scl_edges (xy)
              edges %<>% dplyr::mutate (d = runif (nrow (.))) %>%
                  dplyr::arrange (desc (d))
              ncl <- 12 # desired number of clusters/components
              tree <- scl_spantree (edges) %>%
                  scl_cuttree (edges, ncl = ncl)
              xy <- tibble::tibble (x = xy [, 1], y = xy [, 2])
              g <- scl_plot (tree, xy)
              expect_is (g, "ggplot")
})
