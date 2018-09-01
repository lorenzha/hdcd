#' SplitLoss
#'
#' Caluclates the sum of the loss that results from splitting the segment at the given point
#'
#' @inheritParams FindBestSplit
#' @param split_point Index on which to split the segment
#'
#' @return Sum of the loss for both new segments after split.
SplitLoss <- function(x, split_point, SegmentLossFUN, start, end) {
  SegmentLossFUN(x[1:(split_point - 1), , drop = F], start = start, end = start + (split_point - 2)) +
    SegmentLossFUN(x[split_point:NROW(x), , drop = F], start = start + (split_point - 1), end = end)
}

#' SegmentLoss
#'
#' This closure returns a function which calculates the loss for the given segment of the data.
#'
#' Depending on the desired method and the tuning parameters a different loss function will be parametrized and returned.
#' If method \emph{nodewise_regression} is selected, the additional argument \emph{node} must be supplied to determine on which node (dimension)
#' it should be performed.
#'
#' @inheritParams BinarySegmentation
#' @param n_obs Total number of observations
#' @param ... Further arguments supplied to the select method.
#'
#' @importFrom stats var deviance
#'
#' @return A parametrized loss function
SegmentLoss <- function(n_obs,
                        lambda,
                        penalize_diagonal = FALSE,
                        standardize = TRUE,
                        threshold = 1e-07,
                        method = c("nodewise_regression", "summed_regression", "ratio_regression", "glasso"),
                        ...) {
  args <- list(...)
  mth <- match.arg(method)

  if (mth == "nodewise_regression") {
    p <- args[["node"]]
    stopifnot(length(p) == 1 && is.numeric(p))
  }

  if (mth == "glasso") {
    function(x, ...) {
      obs_count <- NROW(x)
      obs_share <- obs_count / n_obs

      # We need more than one observation to caclculate the covariance matrix
      stopifnot(obs_count > 1)

      n_p <- NCOL(x)

      cov_mat <- (obs_count - 1) / obs_count * cov(x)

      glasso_output <- glasso::glasso(
        cov_mat,
        rho = lambda / sqrt(obs_share) * diag(cov_mat),
        penalize.diagonal = penalize_diagonal,
        thr = threshold
      )

      if (!penalize_diagonal) {
        diag(glasso_output$wi) <- 0
      }

      ((glasso_output$loglik / (-n_p / 2) # Needed to undo transformation of likelihood in glasso package
      - lambda / sqrt(obs_share) * sum(abs(glasso_output$wi))) * obs_share) # Remove regularizer added in glasso package
    }
  } else if (mth == "nodewise_regression") {
    function(x, ...) {
      obs_count <- NROW(x)
      obs_share <- obs_count / n_obs

      # We need more than one observation
      stopifnot(obs_count > 1)

      if (standardize) {
        s_dev_y <- sqrt((obs_count - 1) / obs_count * var(x[, p, drop = F]))
      } else {
        s_dev_y <- 1
      }

      fit <- glmnet::glmnet(
        x[, -p, drop = F], x[, p, drop = F],
        alpha = 1, lambda = lambda / sqrt(obs_share) * s_dev_y,
        thresh = threshold,
        standardize = FALSE
      )
      deviance(fit) / n_obs
    }
  } else if (mth == "ratio_regression") {
    function(x, ...) {
      obs_count <- NROW(x)
      obs_share <- obs_count / n_obs

      # We need more than one observation to caclculate the covariance matrix
      stopifnot(obs_count > 1)

      cov_mat <- (obs_count - 1) / obs_count * cov(x)

      if (standardize) {
        var_x <- diag(cov_mat)
      } else {
        var_x <- 1
      }

      withCallingHandlers({
        glasso_output <- glasso::glasso(
          cov_mat,
          rho = lambda / sqrt(obs_share) * var_x,
          approx = T,
          thr = threshold
        )$wi
      }, warning = HandleGlassoNaN)

      mean_vec <- colMeans(x)
      intercepts <- mean_vec - colSums(glasso_output * mean_vec)
      ss <- t(t(x - x %*% glasso_output) - intercepts)^2
      loss <- obs_count * log(colSums(ss) / obs_count)

      mean(loss / n_obs)
    }
  } else if (mth == "summed_regression") {
    function(x, ...) {
      obs_count <- NROW(x)
      obs_share <- obs_count / n_obs

      # We need more than one observation to caclculate the covariance matrix
      stopifnot(obs_count > 1)

      cov_mat <- (obs_count - 1) / obs_count * cov(x)

      if (standardize) {
        var_x <- diag(cov_mat)
      } else {
        var_x <- 1
      }

      withCallingHandlers({
        glasso_output <- glasso::glasso(
          cov_mat,
          rho = lambda / sqrt(obs_share) * var_x,
          approx = T,
          thr = threshold
        )$wi
      }, warning = HandleGlassoNaN)

      mean_vec <- colMeans(x)
      intercepts <- mean_vec - colSums(glasso_output * mean_vec)
      ss <- t(t(x - x %*% glasso_output) - intercepts)^2
      loss <- colSums(ss)

      mean(loss / n_obs)
    }
  }
}

#### Custom loss functions used with the FUN argument ####

#' Square loss with updates
#'
#' @inheritParams FindBestSplit
#'
#' @export
InitSquaredLoss <- function(x) {
  csum <- cumsum(x)
  csum_2 <- cumsum(x^2)

  n_obs <- NROW(x)

  function(x, start, end) {
    seg_length <- (end - start + 1)

    stopifnot(end >= start && end <= n_obs && start >= 1)

    csum_2_start <- ifelse(start > 1, csum_2[start - 1], 0)
    csum_start <- ifelse(start > 1, csum[start - 1], 0)

    ((csum_2[end] - csum_2_start) -
      (csum[end] - csum_start)^2 / seg_length) / n_obs
  }
}


#' Square loss with naive calculation
#'
#' @inheritParams FindBestSplit
#'
#' @export
InitNaiveSquaredLoss <- function(x) {
  n_obs <- NROW(x)

  function(x, start, end) {
    stopifnot(end >= start && end <= n_obs && start >= 1)

    sum((x - mean(x))^2) / n_obs
  }
}
