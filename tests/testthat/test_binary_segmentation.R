context("Binary segmentation")

test_dat <- system.file("test_data", "test_data.Rdata", package = "hdcd")
load(test_dat)

test_that("binary segmentation finds same changepoints as in original version", {

  # tree <- BinarySegmentation(x = test_data, delta = 0.1, lambda = 0.01, method = "glasso")

  # expect_equal(PruneTreeGamma(tree, gamma = 0.05, 1)[[1]][[1]],
  #            c(42, 87, 130, 174, 216, 259))

  tree <- BinarySegmentation(x = test_data, delta = 0.1, lambda = 0.01, method = "summed_regression")

  expect_equal(
    PruneTreeGamma(tree, gamma = 0.05)[[1]][[1]],
    c(44, 103, 134, 174, 216, 259)
  )

  tree <- BinarySegmentation(x = test_data, delta = 0.1, lambda = 0.01, method = "ratio_regression")

  expect_equal(
    PruneTreeGamma(tree, gamma = 0.05)[[1]][[1]],
    c(41, 87, 129, 174, 217, 259)
  )

  tree <- BinarySegmentation(x = test_data, delta = 0.1, lambda = 0.1, method = "nodewise_regression", node = 10)

  expect_equal(
    PruneTreeGamma(tree, gamma = 0.05)[[1]][[1]],
    c(35, 86, 136, 170, 211, 258)
  )
})


test_that("split function gets split correct", {

  # loss_fun <- SegmentLoss(n_obs = nrow(test_data), lambda = 0.1, const = 0.05, penalize_diagonal = F, method = "glasso")

  # expect_equal(FindBestSplit(x = test_data, n_obs = nrow(test_data), delta = 0.1, use_ternary_search = F,
  #                           SegmentLossFUN = loss_fun)[["opt_split"]],
  #             44)

  loss_fun <- SegmentLoss(n_obs = NROW(test_data), lambda = 0.1, const = 0.05, penalize_diagonal = F, method = "summed_regression")

  expect_equal(
    FindBestSplit(x = test_data, n_obs = NROW(test_data), delta = 0.1, SegmentLossFUN = loss_fun)[["opt_split"]],
    174
  )

  loss_fun <- SegmentLoss(n_obs = NROW(test_data), lambda = 0.1, const = 0.05, penalize_diagonal = F, method = "ratio_regression")

  expect_equal(
    FindBestSplit(x = test_data, n_obs = NROW(test_data), delta = 0.1, SegmentLossFUN = loss_fun)[["opt_split"]],
    174
  )

  loss_fun <- SegmentLoss(n_obs = NROW(test_data), lambda = 0.1, const = 0.05, penalize_diagonal = F, method = "nodewise_regression", node = 10)

  expect_equal(
    FindBestSplit(x = test_data, n_obs = NROW(test_data), delta = 0.1, SegmentLossFUN = loss_fun)[["opt_split"]],
    258
  )
})


test_that("section search", {

  min_points <- 10
  stepsize <- 0.5

  x <- test_data

  SegmentLossFUN <- SegmentLoss(n_obs = NROW(x), lambda = 0.05, method = "summed_regression")

  rec <- SectionSearch()


  result <- rec(
    1:NROW(x), left = 1, right = NROW(x), x = x,
    SegmentLossFUN = SegmentLossFUN, RecFUN = rec, stepsize = stepsize,
    min_points = min_points
  )


  expect_equal(
    result$opt_split,
    97
  )

  expect_equal(
    result$loss,
    2.080238099
  )
})

