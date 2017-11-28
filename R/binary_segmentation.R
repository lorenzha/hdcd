#' BinarySegmentation
#'
#' Uses the binary segmentation algorithmn in order to build a binary tree of the data sequence
#' and given tuning parameters lambda and delta recursively. The tree can then be pruned in order to obtain
#' a changepoint estimate.
#'
#' @inheritParams TernarySearch
#' @param x A n times p data matrix.
#' @param delta Value between 0 and 0.5. Tuning param which determines the minimal segment size
#' proportional to the size of the dataset and hence an upper bound for the number of changepoints.
#' @param lambda Sparsity penality parameter in single lasso fits.
#' @param method Which method should be used for fitting the model? See defaults for possible choices.
#' @param penalize_diagonal Boolean, should the diagonal elements of the precision matrix be penalized?
#' @param threshold The threshold for halting the iteration in glasso or glmnet. In the former it controls the change of single parameters
#' in the latter it controls the total objective value.
#' @param use_ternary_search Use ternary search in each level of the recursion to find a local optimum (EXPERIMENTAL).
#'
#' @return An object of class \strong{bs_tree}
#' @export
#'
#' @examples
#' dat <- SimulateFromModel(CreateModel(n_segments = 2,n = 100,p = 30, ChainNetwork))
#' res <- BinarySegmentation(dat, delta = 0.1, lambda = 0.01, method = "summed_regression")
#' print(res)
BinarySegmentation <- function(x, delta, lambda,
                               method = c("nodewise_regression", "summed_regression", "ratio_regression"),
                               threshold = 1e-7,
                               penalize_diagonal = F,
                               use_ternary_search = F,
                               intervals = 3,
                               ...) {
  mth <- match.arg(method)
  SegmentLossFUN <- SegmentLoss(n_obs = nrow(x), lambda = lambda, penalize_diagonal = penalize_diagonal,
                                method = mth, threshold = threshold, ...)

  tree <- data.tree::Node$new("bs_tree", start = 1, end = nrow(x))

  Rec <- function(x, n_obs, delta, SegmentLossFUN, node = tree, use_ternary_search, intervals) {
    n_selected_obs <- nrow(x)

    if (n_selected_obs / n_obs >= 2 * delta) { # check whether segment is still long enough

      res <- FindBestSplit(x, delta, n_obs, use_ternary_search, SegmentLossFUN, intervals)

      node$min_loss <- min(res[["loss"]])
      node$loss     <- res[["loss"]]
      node$segment_loss <- res[["segment_loss"]]
      split_point <- res[["opt_split"]]

      if (is.na(split_point)) {
        return(NA)
      } else {
        start <- node$start

        child_left <- node$AddChild(as.character(start), start = start,
                                    end = start + split_point - 1)
        Rec(x[1:(split_point - 1),, drop = F], n_obs, delta,
            SegmentLossFUN, child_left, use_ternary_search, intervals)

        child_right <- node$AddChild(as.character(start + split_point - 1), start = start + split_point - 1,
                                     end = start + n_selected_obs - 1)
        Rec(x[split_point:n_selected_obs,, drop = F], n_obs, delta,
            SegmentLossFUN, child_right, use_ternary_search, intervals)
      }
    } else {
      return(NA)
    }
  }
  Rec(x = x, n_obs = nrow(x), delta = delta, SegmentLossFUN = SegmentLossFUN, node = tree,
      use_ternary_search = use_ternary_search, intervals = intervals)
  class(tree) <- c("bs_tree", class(tree))
  tree
}

#' FindBestSplit
#'
#' Uses the SegmentLossFUN function for each possible split on the given segment
#' considering delta and returns the splitting index with the lowest loss
#' as well as the loss for each candidate.
#'
#' @inheritParams BinarySegmentation
#' @inheritParams TernarySearch
#' @param SegmentLossFUN A loss function as created by \code{\link{SegmentLoss}}
FindBestSplit <- function(x, delta, n_obs, use_ternary_search, SegmentLossFUN, intervals) {
  obs_count <- nrow(x)
  min_seg_length <- ceiling(delta * n_obs)

  if (obs_count < 2 * min_seg_length || obs_count < 4)
    return(list(opt_split = NA, loss = NA))

  split_candidates <- seq(
    max(3, min_seg_length + 1),
    min(obs_count - 1, obs_count - min_seg_length + 1), 1
  )

  segment_loss <- SegmentLossFUN(x)

  if (use_ternary_search) {
    result <- TernarySearch(split_candidates, 1, length(split_candidates), x, SegmentLossFUN, intervals)
  } else {
    loss <- sapply(split_candidates, function(y) SplitLoss(x, y, SegmentLossFUN = SegmentLossFUN))
    result <- list(opt_split = split_candidates[which.min(loss)], loss = loss)
  }

  if (round(min(result[["loss"]]), 15) >= round(segment_loss, 15))
    list(opt_split = NA, loss = result[["loss"]], segment_loss = segment_loss)
    else
    list(opt_split = result[["opt_split"]], loss = result[["loss"]], segment_loss = segment_loss)
}

#' TernarySearch
#'
#' Use a ternary search technique to find a local minimum for the next split recursively.
#'
#' Besides the standard ternary approach where the segment is divided into three intervals and evaluated at two points one
#' can also define a more granular grid until the number of intervals is equal to the number
#' of observations in which case all observations will be evaluated and we are in the case of searching for the minimum exactly.
#'
#' @param split_candidates A vector of indices of possible split points
#' @param x An n times p data matrix
#' @param left start index on the left
#' @param right start index on the right
#' @param SegmentLossFUN A loss function as created by \code{\link{SegmentLoss}}
#' @param intervals Number of intervals of the search grid. Only used by ternary search.
TernarySearch <- function(split_candidates, left, right, x, SegmentLossFUN, intervals = 3) {

  # Stopping condition for recursion
  if (abs(left - right) < intervals) {

    # Check all remaining points for the minimum
    inds <- split_candidates[left:right]
    loss <- sapply(inds, function(y) SplitLoss(x, y, SegmentLossFUN = SegmentLossFUN))

    return(list(opt_split = inds[which.min(loss)], loss = min(loss)))
  }

  # Caclulate the new grid points
  grid_points <- left + ceiling((right - left) * seq(1, intervals - 1) / intervals)

  # Evaluate loss for splitting at those points
  grid_loss <- sapply(grid_points, function(y) SplitLoss(x, split_candidates[y], SegmentLossFUN = SegmentLossFUN))

  boundaries <- c(left, grid_points, right)

  min_loss_ind <- which.min(grid_loss) + 1  # Add 1 because we padded the boundary vector
  # Discard outer segments with higher loss and do recursion for remaining variables

  TernarySearch(split_candidates = split_candidates, left = boundaries[min_loss_ind - 1],
                right = boundaries[min_loss_ind + 1], x = x, SegmentLossFUN = SegmentLossFUN,
                intervals = intervals)
}
