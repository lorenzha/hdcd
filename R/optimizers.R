
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
SectionSearch <- function() {

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

  function(split_candidates, left, mid, right, x, SegmentLossFUN, RecFUN, min_points = 3, stepsize = 0.5) {

    # If no mid point is supplied start randomly left or right
    if (missing(mid)) {
      step <- (right - left) * stepsize
      mid <- if (runif(1) <= 0.5) floor(right - step) else ceiling(left + step)
    }

    # Stopping condition for recursion
    if (abs(left - right) < min_points) {

      # Check all remaining points for the minimum
      inds <- split_candidates[left:right]
      loss <- sapply(inds, function(y) SplitLoss(x, y, SegmentLossFUN = SegmentLossFUN))

      return(list(opt_split = inds[which.min(loss)], loss = min(loss)))
    }

    # check if given index has already been computed before if not compute it and store in cache
    f <- function(ind) {
      key <- as.character(ind)
      if (cache_has_key(key)) {
        cache_get(key)
      } else {
        cache_set(key, SplitLoss(
          x = x, split_point = split_candidates[ind],
          SegmentLossFUN = SegmentLossFUN
        ))
        cache_get(key)
      }
    }

    f_mid <- f(mid)

    if (mid - left < right - mid) {
      step <- (right - mid) * stepsize
      new <- floor(right - step)
      f_new <- f(new)
      if (f_new > f_mid) {
        RecFUN(
          split_candidates, left = left, mid = mid, right = new, x = x,
          SegmentLossFUN = SegmentLossFUN, RecFUN = RecFUN, stepsize = stepsize
        )
      } else {
        RecFUN(
          split_candidates, left = mid, mid = new, right = right, x = x,
          SegmentLossFUN = SegmentLossFUN, RecFUN = RecFUN, stepsize = stepsize
        )
      }
    }
    else {
      step <- (mid - left) * stepsize
      new <- ceiling(left + step)
      f_new <- f(new)
      if (f_new > f_mid) {
        RecFUN(
          split_candidates, left = new, mid = mid, right = right, x = x,
          SegmentLossFUN = SegmentLossFUN, RecFUN = RecFUN, stepsize = stepsize
        )
      } else {
        RecFUN(
          split_candidates, left = left, mid = new, right = mid, x = x,
          SegmentLossFUN = SegmentLossFUN, RecFUN = RecFUN, stepsize = stepsize
        )
      }
    }
  }
}
