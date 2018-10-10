
#' SectionSearch
#'
#' Implements a variant of golden section search where the stepsize can be
#' choosen freely.
#'
#' Hence its is possible to trade robustness against computational effciency.
#' This function is a closure that implements a cache for already calculated
#' values and returns a function that can be used to find the global minimum
#' numerically.
#'
#' @inheritParams FindBestSplit
#' @param split_candidates A vector of indices where splits of \code{x} can
#'   occur.
#' @param min_points The number of points left between right and left that
#'   triggers a final evaluation of all remaining split candidates.
#' @param stepsize The stepsize for performing section search, should be in (0,
#'   0.5].
#' @param k_sigma Constant part in threshold \eqn{k\sigma \sqrt(log n)} that
#'   loss needs to differ to decied on where to split. If threshold is not met
#'   the loss for the outer segements will be calculated in the algorithmn
#'   proceed on the side with higher loss (which is equal to the variance for
#'   some loss functions).
#'
#' @export
#' @return Returns a function with arguments left, mid, right, RecFUN where
#'   RecFun should always be set to the object name of the function has been
#'   assigned so the function can call itself recursively.
SectionSearch <- function(x, split_candidates, n_obs, SegmentLossFUN, start, end,
                          min_points = 3,
                          stepsize = 0.5,
                          k_sigma = 0) {

  # Implement cache for storing already calculated values
  cache <- NULL

  seg_loss <- SegmentLossFUN(x, start, end)

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

  key_from_inds <- function(start, end){
    paste(start, end, sep = "-")
  }

  # Initialize the cache
  cache_reset()

  # Initialize alternating splits for equal segment length
  cache_set("left", TRUE)

  function(left, mid, right, RecFUN) {

    n_obs_seg <- right - left + 1

    loss_tolerance <- k_sigma * sqrt(log(n_obs) / n_obs)

    # If no mid point is supplied start on cache status
    if (missing(mid)) {
      step <- n_obs_seg * stepsize
      if (cache_get("left")) {
        mid <- ceiling(left + step)
        cache_set("left", FALSE)
      } else {
        mid <- floor(right - step)
        cache_set("left", TRUE)
      }
    }

    # Stopping condition for recursion
    if (n_obs_seg <= min_points) {

      # Check all remaining points for the minimum
      inds <- split_candidates[left:right]
      loss <- sapply(inds, function(y) SplitLoss(x, y,
          SegmentLossFUN = SegmentLossFUN,
          start = start, end = end
        ))
      gain <- seg_loss - loss
      return(list(opt_split = inds[which.max(gain)], gain = max(gain)))
    }

    # check if given index has already been computed before if not compute it and store in cache
    f <- function(ind) {
      key <- as.character(ind)
      if (cache_has_key(key)) {
        cache_get(key)
      } else {
        cache_set(key, seg_loss - SplitLoss(
          x = x, split_point = split_candidates[ind],
          SegmentLossFUN = SegmentLossFUN, start = start, end = end
        ))
        cache_get(key)
      }
    }

    # check if given loss has already been computed before if not compute it and store in cache
    f_loss <- function(x, start, end) {
      key <- key_from_inds(start, end)
      if (cache_has_key(key)) {
        cache_get(key)
      } else {
        cache_set(key, SegmentLossFUN(x, start, end) / (end - start + 1)) # normalize by segment length
        cache_get(key)
      }
    }

    f_mid <- f(mid)

    if (mid - left == right - mid) {
      dir_left <- cache_get("left")
      cache_set("left", !dir_left)
    }
    else if (mid - left < right - mid) {
      dir_left <- FALSE
    } else {
      dir_left <- TRUE
    }

    if (!dir_left) {
      step <- (right - mid) * stepsize
      new <- floor(right - step)
      f_new <- f(new)
      if (f_new <= f_mid - loss_tolerance) {
        RecFUN(left = left, mid = mid, right = new, RecFUN = RecFUN) # go left
      } else if (f_new > f_mid + loss_tolerance) {
        RecFUN(left = mid, mid = new, right = right, RecFUN = RecFUN) # go right
      } else {
        loss_left  <- f_loss(x[left:mid],
                             start = start + split_candidates[left] - 1, # global indices!
                             end = start + split_candidates[mid] - 1) # global indices!
        loss_right <- f_loss(x[new:right],
                             start = start + split_candidates[new] - 1, # global indices!
                             end = start + split_candidates[right] - 1) # global indices!
        if(loss_left >= loss_right ) {
          RecFUN(left = left, mid = mid, right = new, RecFUN = RecFUN) # go left
        } else {
          RecFUN(left = mid, mid = new, right = right, RecFUN = RecFUN) # go right
        }
      }
    }
    else {
      step <- (mid - left) * stepsize
      new <- ceiling(left + step)
      f_new <- f(new)
      if (f_new <= f_mid - loss_tolerance) {
        RecFUN(left = new, mid = mid, right = right, RecFUN = RecFUN) # go right
      } else if (f_new > f_mid + loss_tolerance) {
        RecFUN(left = left, mid = new, right = mid, RecFUN = RecFUN) # go left
      } else {
        loss_left  <- f_loss(x[left:new],
                             start = start + split_candidates[left] - 1, # global indices!
                             end = start + split_candidates[new] - 1) # global indices!
        loss_right <- f_loss(x[mid:right],
                             start = start + split_candidates[mid] - 1, # global indices!
                             end = start + split_candidates[right] - 1) # global indices!
        if(loss_left >= loss_right ) {
          RecFUN(left = left, mid = new, right = mid, RecFUN = RecFUN) # go left
        } else {
          RecFUN(left = new, mid = mid, right = right, RecFUN = RecFUN) # go right
        }
      }
    }
  }
}
