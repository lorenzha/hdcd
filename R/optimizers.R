
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

  # Initialize alternating splits for equal segment length
  cache_set("left", TRUE)

  function(split_candidates, left, mid, right, x, SegmentLossFUN, RecFUN,
           start, end, min_points = 3, stepsize = 0.5) {

    # If no mid point is supplied start on cache status
    if (missing(mid)) {
      step <- (right - left) * stepsize
      if (cache_get("left")){
        mid <- ceiling(left + step)
        cache_set("left", FALSE)
      }  else {
        mid <- floor(right - step)
        cache_set("left", TRUE)
      }
    }

    # Stopping condition for recursion
    if (abs(left - right) < min_points) {

      # Check all remaining points for the minimum
      inds <- split_candidates[left:right]
      loss <- sapply(inds, function(y) SplitLoss(x, y, SegmentLossFUN = SegmentLossFUN,
                                                 start = start, end = end))

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
          SegmentLossFUN = SegmentLossFUN, start = start, end = end
        ))
        cache_get(key)
      }
    }

    f_mid <- f(mid)

    if (mid - left == right - mid){
      dir_left <- cache_get("left")
      cache_set("left", !dir_left)
    }
    else if (mid - left < right - mid){
      dir_left <- FALSE
    } else {
      dir_left <- TRUE
    }

    if (!dir_left) {
      step <- (right - mid) * stepsize
      new <- floor(right - step)
      f_new <- f(new)
      if (f_new > f_mid) {
        RecFUN(
          split_candidates, left = left, mid = mid, right = new, x = x,
          SegmentLossFUN = SegmentLossFUN, RecFUN = RecFUN, start = start,
          end = end, min_points = min_points, stepsize = stepsize
        )
      } else {
        RecFUN(
          split_candidates, left = mid, mid = new, right = right, x = x,
          SegmentLossFUN = SegmentLossFUN, RecFUN = RecFUN, start = start,
          end = end, min_points = min_points, stepsize = stepsize
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
          SegmentLossFUN = SegmentLossFUN, RecFUN = RecFUN, start = start,
          end = end, min_points = min_points, stepsize = stepsize
        )
      } else {
        RecFUN(
          split_candidates, left = left, mid = new, right = mid, x = x,
          SegmentLossFUN = SegmentLossFUN, RecFUN = RecFUN, start = start,
          end = end, min_points = min_points, stepsize = stepsize
        )
      }
    }
  }
}
