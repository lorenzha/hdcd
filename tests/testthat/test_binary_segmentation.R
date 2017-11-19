context("Binary segmentation")

test_dat <- system.file("data", "test_data.Rdata", package="hdcd")
load(test_dat)

test_that("binary segmentation finds same changepoints as in original version", {

  tree <- BinarySegmentation(x = test_data, delta = 0.1, lambda = 0.01, method = "glasso")

  expect_equal(PruneTreeGamma(tree, gamma_max = 0.05, 1)[[1]][[1]],
               c(44, 87, 130, 174, 216, 259))

  tree <- BinarySegmentation(x = test_data, delta = 0.1, lambda = 0.01, method = "summed_regression")

  expect_equal(PruneTreeGamma(tree, gamma_max = 0.05, 1)[[1]][[1]],
               c(87, 172, 258))

  tree <- BinarySegmentation(x = test_data, delta = 0.1, lambda = 0.01, method = "ratio_regression")

  expect_equal(PruneTreeGamma(tree, gamma_max = 0.05, 1)[[1]][[1]],
               c(45, 87, 130, 172, 210, 258))

  tree <- BinarySegmentation(x=test_data, delta = 0.1, lambda = 0.1, penalize_diagonal = F,
                             use_ternary_search = F, method = "nodewise_regression", p = 10)

  expect_equal(PruneTreeGamma(tree, gamma_max = 0.05, 1)[[1]][[1]],
               c(44, 261))

})


test_that("split function gets split correct", {

  loss_fun <- SegmentLoss(n_obs = nrow(test_data), lambda = 0.1, const = 0.05, penalize_diagonal = F, method = "glasso")

  expect_equal(FindBestSplit(x = test_data, n_obs = nrow(test_data), delta = 0.1, use_ternary_search = F,
                             SegmentLossFUN = loss_fun)[["opt_split"]],
               45)

  loss_fun <- SegmentLoss(n_obs = nrow(test_data), lambda = 0.1, const = 0.05, penalize_diagonal = F, method = "summed_regression")

  expect_equal(FindBestSplit(x = test_data, n_obs = nrow(test_data), delta = 0.1, use_ternary_search = F,
                             SegmentLossFUN = loss_fun)[["opt_split"]],
               131)

  loss_fun <- SegmentLoss(n_obs = nrow(test_data), lambda = 0.1, const = 0.05, penalize_diagonal = F, method = "ratio_regression")

  expect_equal(FindBestSplit(x = test_data, n_obs = nrow(test_data), delta = 0.1, use_ternary_search = F,
                             SegmentLossFUN = loss_fun)[["opt_split"]],
               45)

  loss_fun <- SegmentLoss(n_obs = nrow(test_data), lambda = 0.1, const = 0.05, penalize_diagonal = F, method = "nodewise_regression", p = 10)

  expect_equal(FindBestSplit(x = test_data, n_obs = nrow(test_data), delta = 0.1, use_ternary_search = F,
                             SegmentLossFUN = loss_fun)[["opt_split"]],
               44)
})

test_that("ternary search finds changepoints", {

  tree <- BinarySegmentation(x = test_data, delta = 0.1, lambda = 0.1, method = "glasso", use_ternary_search = T)

  expect_equal(PruneTreeGamma(tree, gamma_max = 0.05, 1)[[1]][[1]],
               c(45, 87, 131, 173, 216, 259))

  tree <- BinarySegmentation(x = test_data, delta = 0.1, lambda = 0.1, method = "summed_regression", use_ternary_search = T)

  expect_equal(PruneTreeGamma(tree, gamma_max = 0.05, 1)[[1]][[1]],
               c(44, 88, 173, 216, 259))

  tree <- BinarySegmentation(x = test_data, delta = 0.1, lambda = 0.1, method = "ratio_regression", use_ternary_search = T)

  expect_equal(PruneTreeGamma(tree, gamma_max = 0.05, 1)[[1]][[1]],
               c(44, 88, 130, 172, 216, 259))

  tree <- BinarySegmentation(x=test_data, delta = 0.1, lambda = 0.01, penalize_diagonal = F,
                             use_ternary_search = T, method = "nodewise_regression", p = 10)

  expect_equal(PruneTreeGamma(tree, gamma_max = 0.05, 1)[[1]][[1]],
               c(78, 179, 257))

})

test_that("ternary search with more than 3 intervals", {


  tree <- BinarySegmentation(x = test_data, delta = 0.1, lambda = 0.1, method = "summed_regression", use_ternary_search = T, intervals = 10)

  expect_equal(PruneTreeGamma(tree, gamma_max = 0.05, 1)[[1]][[1]],
               c(45, 130, 172, 259))

  tree <- BinarySegmentation(x = test_data, delta = 0.1, lambda = 0.1, method = "summed_regression", use_ternary_search = T, intervals = 5)

  expect_equal(PruneTreeGamma(tree, gamma_max = 0.05, 1)[[1]][[1]],
               c(45, 130, 172, 259))

})


