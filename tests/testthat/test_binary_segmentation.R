context("Binary segmentation")

test_dat <- system.file("data", "test_data.Rdata", package="hdcd")
load(test_dat)

test_that("binary segmentation finds same changepoints as in original version", {

  tree <- BinarySegmentation(x = test_data, delta = 0.1, lambda = 0.01, method = "glasso")

  expect_equal(PruneTreeGamma(tree, gamma_max = 0.05, 1)[[1]][[1]],
               c(42, 87, 130, 174, 216, 259))

  tree <- BinarySegmentation(x = test_data, delta = 0.1, lambda = 0.01, method = "summed_regression")

  expect_equal(PruneTreeGamma(tree, gamma_max = 0.05, 1)[[1]][[1]],
               c(44, 103, 134, 174, 216, 259))

  tree <- BinarySegmentation(x = test_data, delta = 0.1, lambda = 0.01, method = "ratio_regression")

  expect_equal(PruneTreeGamma(tree, gamma_max = 0.05, 1)[[1]][[1]],
               c(41, 87, 129, 174, 217, 259))

  tree <- BinarySegmentation(x=test_data, delta = 0.1, lambda = 0.1, penalize_diagonal = F,
                             use_ternary_search = F, method = "nodewise_regression", p = 10)

  expect_equal(PruneTreeGamma(tree, gamma_max = 0.05, 1)[[1]][[1]],
               c(52, 86, 137, 170, 225, 258))

})


test_that("split function gets split correct", {

  loss_fun <- SegmentLoss(n_obs = nrow(test_data), lambda = 0.1, const = 0.05, penalize_diagonal = F, method = "glasso")

  expect_equal(FindBestSplit(x = test_data, n_obs = nrow(test_data), delta = 0.1, use_ternary_search = F,
                             SegmentLossFUN = loss_fun)[["opt_split"]],
               44)

  loss_fun <- SegmentLoss(n_obs = nrow(test_data), lambda = 0.1, const = 0.05, penalize_diagonal = F, method = "summed_regression")

  expect_equal(FindBestSplit(x = test_data, n_obs = nrow(test_data), delta = 0.1, use_ternary_search = F,
                             SegmentLossFUN = loss_fun)[["opt_split"]],
               174)

  loss_fun <- SegmentLoss(n_obs = nrow(test_data), lambda = 0.1, const = 0.05, penalize_diagonal = F, method = "ratio_regression")

  expect_equal(FindBestSplit(x = test_data, n_obs = nrow(test_data), delta = 0.1, use_ternary_search = F,
                             SegmentLossFUN = loss_fun)[["opt_split"]],
               174)

  loss_fun <- SegmentLoss(n_obs = nrow(test_data), lambda = 0.1, const = 0.05, penalize_diagonal = F, method = "nodewise_regression", p = 10)

  expect_equal(FindBestSplit(x = test_data, n_obs = nrow(test_data), delta = 0.1, use_ternary_search = F,
                             SegmentLossFUN = loss_fun)[["opt_split"]],
               225)
})

test_that("ternary search finds changepoints", {

  tree <- BinarySegmentation(x = test_data, delta = 0.1, lambda = 0.1, method = "glasso", use_ternary_search = T)

  expect_equal(PruneTreeGamma(tree, gamma_max = 0.05, 1)[[1]][[1]],
               c(44, 97, 130, 174, 216, 259))

  tree <- BinarySegmentation(x = test_data, delta = 0.1, lambda = 0.1, method = "summed_regression", use_ternary_search = T)

  expect_equal(PruneTreeGamma(tree, gamma_max = 0.05, 1)[[1]][[1]],
               c(44, 86, 130, 174, 217, 260))

  tree <- BinarySegmentation(x = test_data, delta = 0.1, lambda = 0.1, method = "ratio_regression", use_ternary_search = T)

  expect_equal(PruneTreeGamma(tree, gamma_max = 0.05, 1)[[1]][[1]],
               c(44, 86, 130, 174, 216, 260))

  tree <- BinarySegmentation(x=test_data, delta = 0.1, lambda = 0.01, penalize_diagonal = F,
                             use_ternary_search = T, method = "nodewise_regression", p = 10)

  expect_equal(PruneTreeGamma(tree, gamma_max = 0.05, 1)[[1]][[1]],
               c(52, 88, 135, 171, 212, 263))

})

test_that("ternary search with more than 3 intervals", {


  tree <- BinarySegmentation(x = test_data, delta = 0.1, lambda = 0.1, method = "summed_regression", use_ternary_search = T, intervals = 10)

  expect_equal(PruneTreeGamma(tree, gamma_max = 0.05, 1)[[1]][[1]],
               c(44, 86, 130, 174, 217, 260))

  tree <- BinarySegmentation(x = test_data, delta = 0.1, lambda = 0.1, method = "summed_regression", use_ternary_search = T, intervals = 5)

  expect_equal(PruneTreeGamma(tree, gamma_max = 0.05, 1)[[1]][[1]],
               c(44, 86, 130, 174, 217, 260))

})


