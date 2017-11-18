context("Simulation")
set.seed(1)

test_that("simulate gives valid output for chain network", {
  test_mat <- matrix(c(-1.0539190, 1.2555341, 0.4210145, -1.3346208, 0.3021105,
                       -0.4912558, 1.6943262, 0.1882838, -0.1825210, 0.5992027), nrow = 5, ncol = 2)

  expect_equal(SimulateFromModel(CreateModel(1, 5, 2, ChainNetwork)),
               test_mat, tolerance = 1.5e-5)
})


