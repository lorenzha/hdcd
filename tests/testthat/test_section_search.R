context("test-section-search.R")

test_dat <- system.file("test_data", "test_data.Rdata", package = "hdcd")
load(test_dat)

test_that("multivariate section search", {
  expect_equal_to_reference({
    rec <- SectionSearch(test_data, split_candidates = 10:290,
                         n_obs = 300, start = 1, end = 300,
                         SegmentLossFUN = SegmentLoss(n_obs = 300,
                                                      lambda = 0.01,
                                                      method = "summed_regression"))
    rec(left = 1, right = 290, RecFUN = rec)
  }, test_path("reference_objects/section-search"))
})

test_that("univariate section search", {
  expect_equal_to_reference({
    rec <- SectionSearch(test_data[,1], split_candidates = 10:290,
                         n_obs = 300, start = 1, end = 300,
                         SegmentLossFUN = InitSquaredLoss(test_data[,1]),
                         k_sigma = 0)
    rec(left = 1, right = 290, RecFUN = rec)
  }, test_path("reference_objects/section-search-k_sigma"))
})

test_that("univariate section search with k_sigma > 0", {
  expect_equal_to_reference({
    rec <- SectionSearch(test_data[,1], split_candidates = 10:290,
                         n_obs = 300, start = 1, end = 300,
                         SegmentLossFUN = InitSquaredLoss(test_data[,1]),
                         k_sigma = 0.00001)
    rec(left = 1, right = 290, RecFUN = rec)
  }, test_path("reference_objects/section-search-k_sigma"))
})
