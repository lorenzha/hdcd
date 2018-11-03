
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
SectionSearch <- function(split_candidates, n_obs, SegmentLossFUN, start, end,
                          min_points = 3,
                          stepsize = 0.5,
                          k_sigma = 0) {

  seg_loss <- SegmentLossFUN(start, end)

  tol <- k_sigma * sqrt(log(n_obs) / n_obs)

  select_via_variance <- function(start1, end1, start2, end2){
    if (SegmentLossFUN(start1, end1) / (end1-start1 + 1) >= SegmentLossFUN(start2, end2) / (end2 - start2 + 1)){
      SectionSearch_recursive(start1, end1, start2)
    } else {
      SectionSearch_recursive(end1, start2, end2)
    }
  }

  loss <-  rep(NA, n_obs)

  SectionSearch_recursive <- function(cur_left, cur_middle, cur_right){
    if (cur_right - cur_left <= min_points){
      loss[cur_left : cur_right] <- sapply(cur_left : cur_right, function(y) SplitLoss(y, SegmentLossFUN, start, end))
      return(list(gain = seg_loss - loss, opt_split = which.min(loss)))
    }

    if (cur_right - cur_middle > cur_middle - cur_left){
      w <- cur_right - ceiling((cur_right - cur_middle) * stepsize)
      loss[cur_middle] <<- SplitLoss(cur_middle, SegmentLossFUN, start, end)
      loss[w] <<- SplitLoss(w, SegmentLossFUN, start, end)

      if ( loss[w] + tol <= loss[cur_middle] ){
        SectionSearch_recursive(cur_middle, w, cur_right)
      } else if (loss[cur_middle] + tol <= loss[w]) {
        SectionSearch_recursive(cur_left, cur_middle, w)
      } else {
        select_via_variance(cur_left, cur_middle, w, cur_right)
      }
    } else {
      w <- cur_left + ceiling((cur_middle - cur_left) * stepsize)
      loss[w] <<- SplitLoss(w, SegmentLossFUN, start, end)
      loss[cur_middle] <<- SplitLoss(cur_middle, SegmentLossFUN, start, end)

      if ( loss[w] + tol <= loss[cur_middle]){
        SectionSearch_recursive(cur_left, w, cur_middle)
      } else if (loss[cur_middle] + tol <= loss[w]) {
        SectionSearch_recursive(w, cur_middle, cur_right)
      } else {
        select_via_variance(cur_left, w, cur_middle, cur_right)
      }
    }
  }

  left <- split_candidates[1]
  right <- split_candidates[length(split_candidates)]
  SectionSearch_recursive(left, ceiling( (start + stepsize * end)/(1 + stepsize)), right) #generates symmetrical setup in next step
}


  # # Implement cache for storing already calculated values
  # cache <- NULL
  #
  # seg_loss <- SegmentLossFUN(start, end)
  #
  # cache_reset <- function() {
  #   cache <<- new.env(TRUE, emptyenv())
  # }
  #
  # cache_set <- function(key, value) {
  #   assign(key, value, envir = cache)
  # }
  #
  # cache_get <- function(key) {
  #   get(key, envir = cache, inherits = FALSE)
  # }
  #
  # cache_has_key <- function(key) {
  #   exists(key, envir = cache, inherits = FALSE)
  # }
  #
  # key_from_inds <- function(start, end){
  #   paste(start, end, sep = "-")
  # }
  #
  # # Initialize the cache
  # cache_reset()
  #
  # # Initialize alternating splits for equal segment length
  # cache_set("left", TRUE)
  #
  # function(left, mid, right, RecFUN) {
  #
  #   n_obs_seg <- right - left + 1
  #
  #   loss_tolerance <- k_sigma * sqrt(log(n_obs) / n_obs)
  #
  #   # If no mid point is supplied start on cache status
  #   if (missing(mid)) {
  #     step <- n_obs_seg * stepsize
  #     if (cache_get("left")) {
  #       mid <- ceiling(left + step)
  #       cache_set("left", FALSE)
  #     } else {
  #       mid <- floor(right - step)
  #       cache_set("left", TRUE)
  #     }
  #   }
  #
  #   # Stopping condition for recursion
  #   if (n_obs_seg <= min_points) {
  #
  #     # Check all remaining points for the minimum
  #     loss <- sapply(left : right, function(y) SplitLoss(y,
  #         SegmentLossFUN = SegmentLossFUN,
  #         start = start, end = end
  #       ))
  #     gain <- seg_loss - loss
  #     return(list(opt_split = inds[which.max(gain)], gain = max(gain)))
  #   }
  #
  #   # check if given index has already been computed before if not compute it and store in cache
  #   f <- function(ind) {
  #     key <- as.character(ind)
  #     if (cache_has_key(key)) {
  #       cache_get(key)
  #     } else {
  #       cache_set(key, seg_loss - SplitLoss(
  #         split_point = split_candidates[ind],
  #         SegmentLossFUN = SegmentLossFUN, start = start, end = end
  #       ))
  #       cache_get(key)
  #     }
  #   }
  #
  #   # check if given loss has already been computed before if not compute it and store in cache
  #   f_loss <- function(start, end) {
  #     key <- key_from_inds(start, end)
  #     if (cache_has_key(key)) {
  #       cache_get(key)
  #     } else {
  #       cache_set(key, SegmentLossFUN(start, end) / (end - start + 1)) # normalize by segment length
  #       cache_get(key)
  #     }
  #   }
  #
  #   f_mid <- f(mid)
  #
  #   if (mid - left == right - mid) {
  #     dir_left <- cache_get("left")
  #     cache_set("left", !dir_left)
  #   }
  #   else if (mid - left < right - mid) {
  #     dir_left <- FALSE
  #   } else {
  #     dir_left <- TRUE
  #   }
  #
  #   if (!dir_left) {
  #     step <- (right - mid) * stepsize
  #     new <- floor(right - step)
  #     f_new <- f(new)
  #     if (f_new <= f_mid - loss_tolerance) {
  #       RecFUN(left = left, mid = mid, right = new, RecFUN = RecFUN) # go left
  #     } else if (f_new > f_mid + loss_tolerance) {
  #       RecFUN(left = mid, mid = new, right = right, RecFUN = RecFUN) # go right
  #     } else {
  #       loss_left  <- f_loss(x[left:mid],
  #                            start = start + split_candidates[left] - 1, # global indices!
  #                            end = start + split_candidates[mid] - 1) # global indices!
  #       loss_right <- f_loss(x[new:right],
  #                            start = start + split_candidates[new] - 1, # global indices!
  #                            end = start + split_candidates[right] - 1) # global indices!
  #       if(loss_left >= loss_right ) {
  #         RecFUN(left = left, mid = mid, right = new, RecFUN = RecFUN) # go left
  #       } else {
  #         RecFUN(left = mid, mid = new, right = right, RecFUN = RecFUN) # go right
  #       }
  #     }
  #   }
  #   else {
  #     step <- (mid - left) * stepsize
  #     new <- ceiling(left + step)
  #     f_new <- f(new)
  #     if (f_new <= f_mid - loss_tolerance) {
  #       RecFUN(left = new, mid = mid, right = right, RecFUN = RecFUN) # go right
  #     } else if (f_new > f_mid + loss_tolerance) {
  #       RecFUN(left = left, mid = new, right = mid, RecFUN = RecFUN) # go left
  #     } else {
  #       loss_left  <- f_loss(x[left:new],
  #                            start = start + split_candidates[left] - 1, # global indices!
  #                            end = start + split_candidates[new] - 1) # global indices!
  #       loss_right <- f_loss(x[mid:right],
  #                            start = start + split_candidates[mid] - 1, # global indices!
  #                            end = start + split_candidates[right] - 1) # global indices!
  #       if(loss_left >= loss_right ) {
  #         RecFUN(left = left, mid = new, right = mid, RecFUN = RecFUN) # go left
  #       } else {
  #         RecFUN(left = new, mid = mid, right = right, RecFUN = RecFUN) # go right
  #       }
  #     }
  #   }
  # }
# }
