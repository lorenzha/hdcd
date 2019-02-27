
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
SectionSearch <- function(split_candidates, n_obs, split_fun,
                          min_points = 4,
                          stepsize = 0.5,
                          k_sigma = 0) {

  #seg_loss <- SegmentLossFUN(start, end, lambda) #loss over whole segment
  loss <-  rep(NA, n_obs)

  # if(is.na(seg_loss)){
  #   return(list(gain = loss, opt_split = NA))
  # }

  tol <- k_sigma * sqrt(log(n_obs) / n_obs) # for splitting via variance

  select_via_variance <- function(start1, end1, start2, end2){
    stop('select via variance currently not supported')
    if (split_fun(start1, start2, no_split = T) / (end1 - start1) >= split_fun(start2, end2, no_split = T) / (end2 - start2)){
      SectionSearch_recursive_1(start1, end1, start2)
    } else {
      SectionSearch_recursive_1(end1, start2, end2)
    }
  }

  SectionSearch_recursive_2 <- function(cur_left, w_left, w_right, cur_right){
    # If segment [cur_left, w_left) doesn't have enough non-missing observations, loss[w_left] might be NA.
    # In that case, increase w_left by one and try again
    if (is.na(loss[w_left])){
      w_left <- w_left + 1
      loss[w_left] <<- split_fun(w_left) #SplitLoss(w_left, SegmentLossFUN, start, end, lambda)
      SectionSearch_recursive_2(cur_left, w_left, w_right, cur_right)
    # Same for segment [w_right, cur_right)
    } else if (is.na(loss[w_right])){
      w_right <- w_right - 1
      loss[w_right] <<- split_fun(w_right) # SplitLoss(w_right, SegmentLossFUN, start, end, lambda)
      SectionSearch_recursive_2(cur_left, w_left, w_right, cur_right)
    } else if (w_right - w_left >= min_points){ # Check that Segment is long enough
      # if loss[w_left] <= loss[w_right], the inverse is true for gain & thus discard [w_right, cur_right)
      # TODO: This in not symmetrical for tol = 0.
      if(loss[w_left] + tol <= loss[w_right]){
        SectionSearch_recursive_1(cur_left, w_left, w_right)
      } else if (loss[w_right] +tol <= loss[w_left]){
        SectionSearch_recursive_1(w_left, w_right, cur_right)
      } else {
        select_via_variance(cur_left, w_left, w_right, cur_right)
      }
    } else {
      #If Segment is sufficiently small, choose splitpoint via line_search
      loss[cur_left : cur_right] <- sapply(cur_left : cur_right, split_fun) #function(y) SplitLoss(y, SegmentLossFUN, start, end, lambda))
      return(list(gain = seg_loss - loss, opt_split = catch(which.min(loss))))
    }
  }


  SectionSearch_recursive_1 <- function(cur_left, cur_middle, cur_right){

    stopifnot(!is.na(loss[cur_middle]))

    # # if segment has less than min_points points, do line search
    # if (cur_right - cur_left + 1 <= min_points){ #this can be optimized, as loss[cur_middle was already calculated]
    #   loss[cur_left : cur_right] <<- sapply(cur_left : cur_right, function(y) SplitLoss(y, SegmentLossFUN, start, end))
    #   return(list(gain = seg_loss - loss, opt_split = which.min(loss)))
    # }

    # select new midpoint in smaller segment
    if (cur_right - cur_middle > cur_middle - cur_left){
      w <- cur_right - ceiling((cur_right - cur_middle) * stepsize)
      loss[w] <<- split_fun(w) #SplitLoss(w, SegmentLossFUN, start, end, lambda)
      SectionSearch_recursive_2(cur_left, cur_middle, w, cur_right)
    } else {
      w <- cur_left + ceiling((cur_middle - cur_left) * stepsize)
      loss[w] <<- split_fun(w) #SplitLoss(w, SegmentLossFUN, start, end, lambda)
      SectionSearch_recursive_2(cur_left, w, cur_middle, cur_right)
    }
  }

  cur_left <- split_candidates[1]
  cur_right <- split_candidates[length(split_candidates)]
  w_left <- ceiling((cur_left + stepsize * cur_right) / (1 + stepsize))
  w_right <- floor((stepsize * cur_left + cur_right) / (1 + stepsize))
  loss[w_left] <- split_fun(w_left) #SplitLoss(w_left, SegmentLossFUN, start, end, lambda)
  loss[w_right] <- split_fun(w_left) #SplitLoss(w_right, SegmentLossFUN, start, end, lambda)

  SectionSearch_recursive_2(cur_left, w_left, w_right, cur_right)
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
