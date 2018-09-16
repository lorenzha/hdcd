context("test-section-search.R")

x	<- as.matrix(rep(c(0, 1, -1), times = c(80, 10, 10)))


test_that("OBS search using k_sigma finds right split in noiseless case forward", {
  expect_equal({
    GetChangePointsFromLeafs(
      BinarySegmentation(x, 0.05, lambda = 0, optimizer = "section_search",
                          FUN = InitSquaredLoss, control = list(k_sigma = 0.01)))},
    c(81, 91))
})

test_that("OBS search using k_sigma finds right split in noiseless case backward", {
  expect_equal({
    GetChangePointsFromLeafs(
      BinarySegmentation(x[100:1, , drop = F], 0.05, lambda = 0, optimizer = "section_search",
                         FUN = InitSquaredLoss, control = list(k_sigma = 0.01)))},
    c(11, 21))
})

test_that("OBS search using k_sigma does not find right split in noiseless case when k_sigma = 0", {
  expect_equal({
    GetChangePointsFromLeafs(
      BinarySegmentation(x, 0.05, lambda = 0, optimizer = "section_search",
                         FUN = InitSquaredLoss, control = list(k_sigma = 0)))},
    numeric(0))
})

