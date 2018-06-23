#' BinarySegmentation
#'
#' Uses the binary segmentation algorithmn in order to build a binary tree. The tree can then be pruned in order to obtain
#' a changepoint estimate.
#'
#' @param x A n times p data matrix.
#' @param delta Numeric value between 0 and 0.5. Tuning param which determines the minimal segment size
#' proportional to the size of the dataset and hence an upper bound for the number of changepoints.
#' @param lambda Positive numeric value. This is the regularization parameter in the single Lasso fits.
#' @param method Which estimator should be used? Possible choices are
#' \itemize{
#'   \item nodewise_regression: Nodewise regression is based on a single node that needs to be specified with an additional parameter p.
#'   \item summed_regression: Summed nodewise regression sums up the residual variances of nodewise regression over all nodes.
#'   \item ratio_regression: Likelihood ratio based regression sums the pseudo-profile-likelihood over all nodes.
#' }
#' @param penalize_diagonal Boolean, should the diagonal elements of the precision matrix be penalized?
#' @param optimizer Which search technique should be used for performing individual splits in the binary segmentation alogrithm? Possible choices are
#' \itemize{
#'   \item line_search: Exhaustive linear search. All datapoints are evaluated and the maximum is returned.
#'   \item ternary_search: Iteratively cuts the search space according by a fixed ratio and approximately finds a local maximum.
#'   \item section_search: Iteratively cuts the search space according by a flexible ratio as determined by stepsize in control parameter and approximately finds a local maximum.
#' }
#' @param control A list with parameters that is accessed by the optimizer.
#' \itemize{
#'   \item stepsize: Numeric value between 0 and 0.5. Used by section search.
#'   \item intervals: Integer value larger than 3. Used by ternary search.
#' }
#' @param standardize Boolean. If TRUE the penalty parameter \eqn{\lambda} will be the standard deviation for every dimension in the single Lasso fits.
#' @param threshold The threshold for halting the iteration in glasso or glmnet. In the former it controls the absolute change of single parameters in the latter it controls the total objective value.
#' @param verbose Boolean. If TRUE additional information will be printed.
#' @param ... Supply additional arguments for a specific method (e.g. p for nodewise_regression)
#'
#' @return An object of class \strong{bs_tree}.
#' @export
#'
#' @examples
#' dat <- SimulateFromModel(CreateModel(n_segments = 2,n = 100,p = 30, ChainNetwork))
#' res <- BinarySegmentation(dat, delta = 0.1, lambda = 0.01, method = "summed_regression")
#' print(res)
BinarySegmentation <- function(x, delta, lambda,
                               method = c("nodewise_regression", "summed_regression", "ratio_regression"),
                               penalize_diagonal = F,
                               optimizer = c("line_search", "ternary_search", "section_search"),
                               control = NULL,
                               standardize = T,
                               threshold = 1e-7,
                               verbose = FALSE,
                               ...) {
  SegmentLossFUN <- SegmentLoss(
    n_obs = nrow(x), lambda = lambda, penalize_diagonal = penalize_diagonal,
    method = method, standardize = standardize, threshold = threshold, ...
  )

  tree <- data.tree::Node$new("bs_tree", start = 1, end = nrow(x))
  class(tree) <- c("bs_tree", class(tree))


  Rec <- function(x, n_obs, delta, SegmentLossFUN, node, optimizer) {
    n_selected_obs <- nrow(x)

    if (verbose) print(tree)

    if (n_selected_obs / n_obs >= 2 * delta) { # check whether segment is still long enough

      res <- FindBestSplit(x, delta, n_obs, optimizer, control, SegmentLossFUN)

      node$min_loss <- min(res[["loss"]])
      node$loss <- res[["loss"]]
      node$segment_loss <- res[["segment_loss"]]
      split_point <- res[["opt_split"]]

      if (is.na(split_point)) {
        return(NA)
      } else {
        start <- node$start

        child_left <- node$AddChild(
          as.character(start), start = start,
          end = start + split_point - 1
        )
        Rec(
          x[1:(split_point - 1), , drop = F], n_obs, delta,
          SegmentLossFUN, child_left, optimizer
        )

        child_right <- node$AddChild(
          as.character(start + split_point - 1), start = start + split_point - 1,
          end = start + n_selected_obs - 1
        )
        Rec(
          x[split_point:n_selected_obs, , drop = F], n_obs, delta,
          SegmentLossFUN, child_right, optimizer
        )
      }
    } else {
      return(NA)
    }
  }
  Rec(
    x = x, n_obs = nrow(x), delta = delta, SegmentLossFUN = SegmentLossFUN, node = tree,
    optimizer = optimizer
  )
  tree
}

#' FindBestSplit
#'
#' Uses the SegmentLossFUN function for each possible split on the given segment
#' considering delta and returns the splitting index with the lowest loss
#' as well as the loss for each candidate.
#'
#' @inheritParams BinarySegmentation
#' @param SegmentLossFUN A loss function as created by \code{\link{SegmentLoss}}
#'
FindBestSplit <- function(x, delta, n_obs, optimizer = c("line_search", "ternary_search", "section_search"), control, SegmentLossFUN) {
  opt <- match.arg(optimizer)

  obs_count <- nrow(x)
  min_seg_length <- ceiling(delta * n_obs)

  if (obs_count < 2 * min_seg_length || obs_count < 4) {
    return(list(opt_split = NA, loss = NA))
  }

  split_candidates <- seq(
    max(3, min_seg_length + 1),
    min(obs_count - 1, obs_count - min_seg_length + 1), 1
  )

  segment_loss <- SegmentLossFUN(x)

  switch(opt,
    "line_search" = {
      loss <- sapply(split_candidates, function(y) SplitLoss(x, y, SegmentLossFUN = SegmentLossFUN))
      result <- list(opt_split = split_candidates[which.min(loss)], loss = loss)
    },
    "ternary_search" = {
      intervals <- control[["intervals"]]
      if (is.null(intervals) || !all.equal(intervals, as.integer(intervals)))
        intervals <- 3 # set default value if necessary
      result <- TernarySearch(
        split_candidates, 1, length(split_candidates), x, SegmentLossFUN, intervals
        )
    },
    "section_search" = {
      rec <- SectionSearch()
      min_points <- control[["min_points"]]
      stepsize <- control[["stepsize"]]
      if (is.null(stepsize) || stepsize <= 0)
        stepsize <- 0.1 # set default value if necessary
      result <- rec(
        split_candidates, left = 1, right = length(split_candidates), x = x,
        SegmentLossFUN = SegmentLossFUN, RecFUN = rec, stepsize = stepsize,
        min_points = min_points
      )
    }
  )

  if (round(min(result[["loss"]]), 15) >= round(segment_loss, 15)) {
    list(opt_split = NA, loss = result[["loss"]], segment_loss = segment_loss)
  } else {
    list(opt_split = result[["opt_split"]], loss = result[["loss"]], segment_loss = segment_loss)
  }
}
