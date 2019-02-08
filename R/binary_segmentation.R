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
#' @param FUN A loss function with formal arguments, \code{x} and \code{lambda},
#' which returns a function handle taking arguments \code{start} and \code{end}, representing the loss for the
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
BinarySegmentation <- function(x, x_test = NULL, test_inds = NULL, lambda = NULL, gamma = NULL, delta = NULL,
                               method = c("nodewise_regression", "summed_regression", "ratio_regression"),
                               NA_method = 'complete_observations',
                               optimizer = c("line_search", "section_search"),
                               control = NULL,
                               FUN = NULL,
                               ...) {

  max_depth <- control_get(control, "max_depth", Inf)
  verbose <- control_get(control, "verbose", TRUE)
  cv_outer <- control_get(control, "cv_outer", FALSE)
  cv_inner <- control_get(control, "cv_inner", FALSE)
  stop_early <- control_get(control, "stop_early", FALSE)
  n_obs <- NROW(x)

  if (is.null(FUN)) {
    SegmentLossFUN <- SegmentLoss(
      x, lambda = lambda, method = method, NA_method = NA_method, ...
    )
  } else {
    stopifnot(c("x") %in% methods::formalArgs(FUN)) ###TODO adapt such that this works with regression
    if ("lambda" %in% methods::formalArgs(FUN)){
      SegmentLossFUN <- FUN(x, lambda = lambda)
    } else {
      SegmentLossFUN <- FUN(x)
    }
  }
  stopifnot(c("start", "end") %in% methods::formalArgs(SegmentLossFUN))

  # initiate tree
  tree <- data.tree::Node$new("bs_tree", start = 1, end = NROW(x))
  class(tree) <- c("bs_tree", class(tree))

  if(cv_inner){
    # calculates cv_loss for whole training data
    cv_loss_global_output <- cv_loss(x, n_obs = n_obs, x_test = x_test,  method = method, NA_method = NA_method, control = control)
    tree$cv_train_loss <- cv_loss_global_output$train_loss
    tree$cv_loss <-  cv_loss_global_output$loss
    tree$cv_test_loss <- cv_loss_global_output$test_loss
    tree$cv_lambda_opt <- cv_loss_global_output$lambda_opt
  }

  BinarySegmentation_recursive <- function(delta, n_obs, SegmentLossFUN, node, optimizer) {

    n_selected_obs <- node$end - node$start + 1

    # check whether segment is still long enough & tree not already to deep
    if (n_selected_obs / n_obs < 2 * delta | length(node$path) - 1 >= max_depth){
      return(NA)
    } else {
      res <- FindBestSplit(
        node$start, node$end, delta, n_obs,
        SegmentLossFUN, control, optimizer,
        gamma
      )

      max_gain <- max(res[["gain"]], na.rm = T)
      node$max_gain <- ifelse(is.finite(max_gain), max_gain, NA)
      node$gain <- res[["gain"]]
      split_point <- res[["opt_split"]]

      if (is.na(split_point)) { #FindBestSplit will return NA if max(gain) <= gamma
        return(NA)
      } else {

        start <- node$start
        end <- node$end

        child_left <- node$AddChild(
          as.character(start),
          start = start,
          end = split_point - 1
        )

        child_right <- node$AddChild(
          as.character(split_point),
          start = split_point,
          end = end
        )

        if(cv_inner){
          if(cv_outer){
            x_test_left <- x_test[test_inds[start] : test_inds[split_point - 1], ]
            x_test_right <- x_test[test_inds[split_point] : test_inds[end], ]
          } else {
            x_test_left <- x_test_right <- NULL
          }
          cv_loss_left <- cv_loss(x[start : (split_point - 1), ],
                                               n_obs = n_obs,
                                               x_test = x_test_left,
                                               method = method, NA_method = NA_method, control = control)
          cv_loss_right<- cv_loss(x[split_point : end, ],
                                    n_obs = n_obs,
                                    x_test =  x_test_right,
                                    method = method, NA_method = NA_method, control = control)
          child_left$cv_train_loss <-  cv_loss_left$train_loss
          child_left$cv_loss <- cv_loss_left$loss
          child_left$cv_test_loss <-  cv_loss_left$test_loss
          child_left$cv_lambda_opt <- cv_loss_left$lambda_opt
          child_right$cv_train_loss <-  cv_loss_right$train_loss
          child_right$cv_loss <- cv_loss_right$loss
          child_right$cv_test_loss <-  cv_loss_right$test_loss
          child_right$cv_lambda_opt <- cv_loss_right$lambda_opt
          node$cv_train_improvement <- (node$cv_train_loss - cv_loss_left$train_loss - cv_loss_right$train_loss) / n_selected_obs
          node$cv_train_improvement_biased <- max(node$cv_loss - cv_loss_left$loss - cv_loss_right$loss) / n_selected_obs

          node$cv_test_improvement <- (node$cv_test_loss - cv_loss_left$test_loss - cv_loss_right$test_loss) / (test_inds[end] - test_inds[start])
        }

        if(!stop_early || !cv_inner || node$cv_train_improvement > 0){
          BinarySegmentation_recursive(delta, n_obs = n_obs, SegmentLossFUN, child_left, optimizer)
          BinarySegmentation_recursive(delta, n_obs = n_obs, SegmentLossFUN, child_right, optimizer)
        }
      }
    }
  }
  BinarySegmentation_recursive(delta = delta, n_obs = n_obs, SegmentLossFUN = SegmentLossFUN,
    node = tree, optimizer = optimizer)
  tree
}

#' FindBestSplit
#'
#' Takes a segment of the data and dispatches the choosen \code{method} to the different optimizers.
#'
#'
#' @param start The start index of the given segment saved in \code{SegmentLossFUN}.
#' @param end The end index of the given segment saved in \code{SegmentLossFUN}.
#' @param n_obs The number of observations in the data set.
#' @param SegmentLossFUN A loss function as created by closure \code{\link{SegmentLoss}}.
#' @param control A list with parameters that is accessed by the selected
#'  optimizer: \itemize{ \item \strong{stepsize}: Numeric value between 0 and
#'  0.5. Used by section search. \item \strong{min_points}: Integer value larger
#'  than 3. Used by section search.}
#'  @param optimizer Which search technique should be used for performing
#'  individual splits in the binary segmentation alogrithm? Possible choices are
#'  \itemize{ \item \strong{line_search}: Exhaustive linear search. All possivle
#'  split candidates are evaluated and the index with maximal loss reduction is
#'  returned. \item \strong{section_search}: Iteratively cuts the search space
#'  according by a flexible ratio as determined by parameter \code{stepsize} in
#'  \code{control} parameter list and approximately finds an index at a local
#'  maximum. See Haubner (2018) for details. }
#'  @param gamma Tuning parameter that controls sensitivity to finding splits.
FindBestSplit <- function(start, end, delta, n_obs, SegmentLossFUN,
                          control = NULL,
                          optimizer = c("line_search", "section_search"),
                          gamma = 0
                          ) {

  opt <- match.arg(optimizer)

  obs_count <- end - start + 1
  min_seg_length <- max(2, ceiling(delta * n_obs))

  if (obs_count < 2 * min_seg_length) {
    return(list(opt_split = NA, gain = NA))
  }

  split_candidates <- seq(
    start + min_seg_length, # need at least two points to estimate loss
    end - min_seg_length + 1, 1
  )

  split_loss <- rep(NA, n_obs)

  switch(opt,
    "line_search" = {
      split_loss[split_candidates] <- sapply(
        split_candidates,
        function(y) SplitLoss(y,
            SegmentLossFUN = SegmentLossFUN,
            start = start,
            end = end
          )
      )
      gain <- SegmentLossFUN(start, end) - split_loss
      result <- list(opt_split = catch(which.max(gain)), gain = gain) #this will return NA if all(is.na(gain))
    },
    "section_search" = {
      min_points <- control[["min_points"]]
      stepsize <- control[["stepsize"]]
      k_sigma <- control[["k_sigma"]]

      if (is.null(stepsize) || stepsize <= 0) {
        stepsize <- 0.5
      } # set default value if necessary
      if (is.null(min_points) || min_points < 3) {
        min_points <- 3
      } # set default value if necessary
      if (is.null(k_sigma) || k_sigma < 0) {
        k_sigma <- 0
      } # set default value if necessary

      result <- SectionSearch(split_candidates = split_candidates, n_obs = n_obs,
                              SegmentLossFUN = SegmentLossFUN, start = start, end = end,
                              min_points = min_points, stepsize = stepsize, k_sigma = k_sigma)
    }
  )


  if (is.na(result[["opt_split"]]) || result[["gain"]][result[["opt_split"]]] - ifelse(is.null(gamma), 0, gamma) <= 0){
    list(opt_split = NA, gain = result[["gain"]])
  } else {
    list(opt_split = result[["opt_split"]], gain = result[["gain"]])
  }
}
