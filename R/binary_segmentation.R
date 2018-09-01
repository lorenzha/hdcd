#' BinarySegmentation
#'
#' Applies the binary segmentation algorithmn by recursively calling
#' \code{\link{FindBestSplit}} in order to build a binary tree. The tree can then
#' be pruned using \code{\link{PruneTreeGamma}} in order to obtain a changepoint
#' estimate. Typically this function is not used directly but the interface
#' \code{\link{hdcd}}.
#'
#' @param x A n times p matrix or data frame.
#' @param delta Numeric value between 0 and 0.5. This tuning parameter determines
#'  the minimal segment size proportional to the size of the dataset and hence
#'  an upper bound for the number of changepoints (roughly \eqn{1/\delta}).
#' @param lambda Positive numeric value. This is the regularization parameter in
#'  the single Lasso fits. This value is ignored if FUN is not NULL.
#' @param method Which estimator should be used? Possible choices are \itemize{
#'  \item \strong{nodewise_regression}: Nodewise regression is based on a single
#'  node that needs to be specified with an additional parameter \code{node}
#'  pointing to the column index of the node of interest. Uses
#'  \code{\link[glmnet]{glmnet}} internally. See Kovács (2016) for details.
#'  \item \strong{summed_regression}: Summed nodewise regression sums up the
#'  residual variances of nodewise regression over all nodes. Uses
#'  \code{\link[glasso]{glasso}} internally. See Kovács (2016) for details.
#'  \item \strong{ratio_regression}: Likelihood ratio based regression sums the
#'  pseudo-profile-likelihood over all nodes. Uses \code{\link[glasso]{glasso}}
#'  internally. See Kovács (2016) for details. \item \strong{glasso}: The
#'  graphical Lasso uses the approach of Friedman et al (2007). In contrast to
#'  the other approaches the exact likelihood the whole graphical model is
#'  computed and used as loss. } This value is ignored if \code{FUN} is not
#'  \code{NULL}.
#' @param penalize_diagonal Boolean, should the diagonal elements of the
#'  precision matrix be penalized by \eqn{\lambda}? This value is ignored if FUN
#'  is not NULL.
#' @param optimizer Which search technique should be used for performing
#'  individual splits in the binary segmentation alogrithm? Possible choices are
#'  \itemize{ \item \strong{line_search}: Exhaustive linear search. All possivle
#'  split candidates are evaluated and the index with maximal loss reduction is
#'  returned. \item \strong{section_search}: Iteratively cuts the search space
#'  according by a flexible ratio as determined by parameter \code{stepsize} in
#'  \code{control} parameter list and approximately finds an index at a local
#'  maximum. See Haubner (2018) for details. }
#' @param control A list with parameters that is accessed by the selected
#'  optimizer: \itemize{ \item \strong{stepsize}: Numeric value between 0 and
#'  0.5. Used by section search. \item \strong{min_points}: Integer value larger
#'  than 3. Used by section search.}
#' @param standardize Boolean. If TRUE the penalty parameter \eqn{\lambda} will
#'  be adjusted for every dimension in the single Lasso fits according to the
#'  standard deviation in the data.
#' @param threshold The threshold for halting the iteration in
#'  \code{\link[glasso]{glasso}} or \code{\link[glmnet]{glmnet}}. In the former
#'  it controls the absolute change of single parameters in the latter it
#'  controls the total objective value. This value is ignored if FUN is not
#'  NULL.
#' @param verbose Boolean. If TRUE additional information will be printed.
#' @param FUN A loss function with formal arguments, \code{x}, \code{n_obs} and
#'  \code{standardize} which returns a scalar representing the loss for the
#'  segment the function is applied to.
#' @param ... Supply additional arguments for a specific method (e.g. \code{node}
#'  for \strong{nodewise_regression}) or own loss function \code{FUN}
#'
#' @section References:
#'
#'  Friedman, J., Hastie, T. & Tibshirani, R. Sparse inverse covariance
#'  estimation with the graphical lasso. Biostatistics 9, 432–441 (2008).
#'
#'  Haubner, L. Optimistic binary segmentation: A scalable approach to
#'  changepoint detection in high-dimensional graphical models. (Seminar for
#'  Statistics, ETH Zurich, 2018).
#'
#'  Kovács, S. Changepoint detection for high-dimensional covariance matrix
#'  estimation. (Seminar for Statistics, ETH Zurich, 2016).
#'
#' @return An object of class \strong{bs_tree} and \strong{Node} (as defined in
#'  \code{\link[data.tree]{Node}}).
#' @export
#'
#' @examples
#' # Use summed regression loss function and ChainNetwork
#' dat <- SimulateFromModel(CreateModel(n_segments = 2,n = 50,p = 30, ChainNetwork))
#' res <- BinarySegmentation(dat, delta = 0.1, lambda = 0.01, method = "summed_regression")
#' print(res)
#'
#' # Define your own loss function and pass it to the BinarySegmentation algorithm
#'
#' InitNaiveSquaredLoss <- function(x){
#'
#'   n_obs <- NROW(x)
#'
#'   function(x, start, end){
#'
#'     stopifnot(end >= start && end <= n_obs && start >= 1)
#'
#'     sum((x - mean(x))^2) / n_obs
#'
#'   }
#' }
#'
#' p <- 5
#' n <- 20
#' mean_vecs <- list(rep(0, p), c(rep(1, p-2), rep(5, 2)), rep(-0.5, p))
#'
#' model <- CreateModel(3, n, p, DiagMatrix, equispaced = TRUE, mean_vecs = mean_vecs)
#'
#' x <- SimulateFromModel(model)
#'
#' res <- BinarySegmentation(x, delta = 0.1, lambda = 0.01, FUN = InitNaiveSquaredLoss,
#' optimizer = "section_search")
#' print(res)
#'
#' res <- BinarySegmentation(x, delta = 0.1, lambda = 0.01, FUN = InitNaiveSquaredLoss,
#' optimizer = "line_search")
#' print(res)
#'
BinarySegmentation <- function(x, delta, lambda,
                               gamma = 0,
                               method = c("nodewise_regression", "summed_regression", "ratio_regression"),
                               penalize_diagonal = F,
                               optimizer = c("line_search", "section_search"),
                               control = NULL,
                               standardize = T,
                               threshold = 1e-7,
                               verbose = FALSE,
                               FUN = NULL,
                               ...) {
  n_obs <- NROW(x)

  if (is.null(FUN)) {
    SegmentLossFUN <- SegmentLoss(
      n_obs = n_obs, lambda = lambda, penalize_diagonal = penalize_diagonal,
      method = method, standardize = standardize, threshold = threshold, ...
    )
  } else {
    stopifnot(c("x") %in% methods::formalArgs(FUN))
    SegmentLossFUN <- FUN(x)
    stopifnot(c("x", "start", "end") %in% methods::formalArgs(SegmentLossFUN))
  }


  tree <- data.tree::Node$new("bs_tree", start = 1, end = NROW(x))
  class(tree) <- c("bs_tree", class(tree))


  Rec <- function(x, n_obs, delta, SegmentLossFUN, node, optimizer) {
    n_selected_obs <- NROW(x)

    if (verbose) print(tree)

    if (n_selected_obs / n_obs >= 2 * delta) { # check whether segment is still long enough

      res <- FindBestSplit(
        x, node$start, node$end, delta, n_obs, control,
        SegmentLossFUN, optimizer, gamma
      )

      node$min_loss <- min(res[["loss"]])
      node$loss <- res[["loss"]]
      node$segment_loss <- res[["segment_loss"]]
      split_point <- res[["opt_split"]]

      if (is.na(split_point)) {
        return(NA)
      } else {
        start <- node$start

        child_left <- node$AddChild(
          as.character(start),
          start = start,
          end = start + split_point - 2
        )
        Rec(
          x[1:(split_point - 1), , drop = F], n_obs, delta,
          SegmentLossFUN, child_left, optimizer
        )

        child_right <- node$AddChild(
          as.character(start + split_point - 1),
          start = start + split_point - 1,
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
    x = x, n_obs = n_obs, delta = delta, SegmentLossFUN = SegmentLossFUN,
    node = tree, optimizer = optimizer
  )
  tree
}

#' FindBestSplit
#'
#' Takes a segment of the data and dispatches the choosen \code{method} to the different optimizers.
#'
#' @inheritParams BinarySegmentation
#' @param n_obs The number of observations in the data set.
#' @param SegmentLossFUN A loss function as created by closure \code{\link{SegmentLoss}}.
#' @param start The start index of the given segment \code{x}.
#' @param end The end index of the given segment \code{x}.
FindBestSplit <- function(x, start, end, delta, n_obs, control, SegmentLossFUN,
                          optimizer = c("line_search", "section_search"),
                          gamma = 0) {
  opt <- match.arg(optimizer)

  obs_count <- NROW(x)
  min_seg_length <- ceiling(delta * n_obs)

  if (obs_count < 2 * min_seg_length || obs_count < 4) {
    return(list(opt_split = NA, loss = NA))
  }

  split_candidates <- seq(
    max(3, min_seg_length + 1),
    min(obs_count - 1, obs_count - min_seg_length + 1), 1
  )

  segment_loss <- SegmentLossFUN(x, start, end)
  switch(opt,
    "line_search" = {
      loss <- sapply(
        split_candidates,
        function(y) SplitLoss(x, y,
            SegmentLossFUN = SegmentLossFUN,
            start = start,
            end = end
          )
      )
      result <- list(opt_split = split_candidates[which.min(loss)], loss = loss)
    },
    "section_search" = {
      rec <- SectionSearch()
      min_points <- control[["min_points"]]
      stepsize <- control[["stepsize"]]
      if (is.null(stepsize) || stepsize <= 0) {
        stepsize <- 0.5
      } # set default value if necessary
      if (is.null(min_points) || min_points < 3) {
        min_points <- 3
      } # set default value if necessary
      result <- rec(
        split_candidates,
        left = 1, right = length(split_candidates), x = x,
        SegmentLossFUN = SegmentLossFUN, RecFUN = rec, stepsize = stepsize,
        min_points = min_points, start = start, end = end
      )
    }
  )

  if (round(min(result[["loss"]]), 15) + gamma >= round(segment_loss, 15)) {
    list(opt_split = NA, loss = result[["loss"]], segment_loss = segment_loss)
  } else {
    list(opt_split = result[["opt_split"]], loss = result[["loss"]], segment_loss = segment_loss)
  }
}
