context("Simulation")
set.seed(1)

test_that("simulate gives valid output for chain network", {
  test_mat <- matrix(c(
    0.6972238, -2.4047035, -0.2672252, 0.2590463, -0.8504294,
    -1.4957939, 1.7819398, 0.5975326, -1.8941850, 0.4287759
  ), nrow = 5, ncol = 2)

  expect_equal(
    SimulateFromModel(CreateModel(1, 5, 2, ChainNetwork)),
    test_mat, tolerance = 1.5e-5
  )
})
