#' PostProcess
#'
#' Adjust the position of each changepoint by refitting the model in the respective segments.
#'
#' @inheritParams BinarySegmentation
#' @param cpts A numeric vector containing the found changepoints. Can be of length zero if no changepoints have been found.
#'
#' @return A numeric vector containing the final segment boundaries
#' @export
#'
PostProcess <- function(x, cpts, delta, lambda, method) {
  n_obs <- nrow(x)

  seg_boundaries <- c(1, cpts, n_obs)

  if (length(seg_boundaries) < 3) {
    return(seg_boundaries)
  } else {
    for (i in seq(2, length(seg_boundaries) - 1)) {
      seg_inds <- seg_boundaries[i - 1]:(seg_boundaries[i + 1] - 1)

      opt_split_ind <- FindBestSplit(
        x[seg_inds, ],
        delta = delta, n_obs = n_obs, optimizer = "line_search", control = NULL,
        SegmentLossFUN = SegmentLoss(n_obs = n_obs, lambda = lambda, method = method)
      )$opt_split

      seg_boundaries[i] <- if (is.na(opt_split_ind)) NA else seg_inds[opt_split_ind]
    }
  }
  seg_boundaries[which(!is.na(seg_boundaries))]
}
