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


#' SectionSearch
#'
#' Implements a variant of golden section search where the stepsize can be choosen
#' freely.
#'
#' Hence its is possible to trade robustness against computational effciency. This function is
#' a closure that implements a cache for already calculated values and returns a function that can
#' be used to find the global minimum numerically.
#'
#'
#' @return
#' Returns a function with arguments split_candidates, left, mid, right, x, SegmentLossFUN, RecFUN, stepsize
#' where RecFun should always be set to the object name of the function has been assigned so the function can
#' call itself recursively.
SectionSearch <- function(){

  # Implement cache for storing already calculated values
  cache <- NULL

  cache_reset <- function() {
    cache <<- new.env(TRUE, emptyenv())
  }

  cache_set <- function(key, value) {
    assign(key, value, envir = cache)
  }

  cache_get <- function(key) {
    get(key, envir = cache, inherits = FALSE)
  }

  cache_has_key <- function(key) {
    exists(key, envir = cache, inherits = FALSE)
  }

  # Initialize the cache
  cache_reset()

  function(split_candidates, left, mid, right, x, SegmentLossFUN, RecFUN, stepsize = 0.1) {
    step <- (right - left) * stepsize

    # If no mid point is supplied start randomly left or right
    if (missing(mid))
      mid <- if (runif(1) <= 0.5) ceiling(right - step) else floor(left + step)

    # Stopping condition for recursion
    if (step < 1) {

      # Check all remaining points for the minimum
      inds <- split_candidates[left:right]
      loss <- sapply(inds, function(y) SplitLoss(x, y, SegmentLossFUN = SegmentLossFUN))

      return(list(opt_split = inds[which.min(loss)], loss = min(loss)))
    }

    # check if given index has already been computed before if not compute it and store in cache
    f <- function(ind){
      key <- as.character(ind)
      if (cache_has_key(key))
        cache_get(key)
      else {
        cache_set(key, SplitLoss(x = x, split_point = split_candidates[ind],
                                 SegmentLossFUN = SegmentLossFUN))
        cache_get(key)
      }
    }

    f_mid   <- f(mid)

    if (mid - left < right - mid){
      new <- ceiling(right - step)
      f_new <- f(new)
      if (f_new > f_mid)
        RecFUN(split_candidates, left = left, mid = mid, right = new, x = x,
               SegmentLossFUN = SegmentLossFUN, RecFUN = RecFUN, stepsize = stepsize)
      else
        RecFUN(split_candidates, left = mid, mid = new, right = right, x = x,
               SegmentLossFUN = SegmentLossFUN, RecFUN = RecFUN, stepsize = stepsize)
    }
    else{
      new <- floor(left + step)
      f_new <- f(new)
      if (f_new > f_mid)
        RecFUN(split_candidates, left = new, mid = mid, right = right, x = x,
               SegmentLossFUN = SegmentLossFUN, RecFUN = RecFUN, stepsize = stepsize)
      else
        RecFUN(split_candidates, left = left, mid = new, right = mid, x = x,
               SegmentLossFUN = SegmentLossFUN, RecFUN = RecFUN, stepsize = stepsize)
    }
  }
}


