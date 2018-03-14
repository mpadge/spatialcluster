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
              edges %<>% dplyr::mutate (d = runif (nrow (.)),
                                        id = seq (nrow (.))) %>%
                  dplyr::arrange (desc (d))
              tree <- scl_spantree (edges)
              expect_true (nrow (tree) < nrow (edges))
})
