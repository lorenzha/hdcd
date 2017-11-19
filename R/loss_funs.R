#' SplitLoss
#'
#' Caluclates the sum of the loss that results from splitting the segment at the given point
#'
#' @param opt_split Index on which to split the segment
#'
#' @return Sum of the loss for both new segments after split.
SplitLoss <- function(x, split_point, SegmentLossFUN) {
  SegmentLossFUN(x[1:(split_point-1),, drop = F]) +
    SegmentLossFUN(x[split_point:nrow(x),, drop = F])
}

#' SegmentLoss
#'
#' This closure returns a function which calculates the loss for the given segment of the data.
#'
#' Depending on the desired method and the tuning paramtersa different loss function will be parametrized and returned.
#' If method \emph{nodewise_regression} is selected, additional argument \emph{p} can be supplied to determine on which node (dimension)
#' it should be performed.
#'
#' @inheritParams BinarySegmentation
#'
#' @return A parametrized loss function
#'
#' @examples
#' dat <- SimulateFromModel(CreateModel(n_segments = 1,n = 100,p = 30, ChainNetwork))
#' lossFUN <- SegmentLoss(100, 0.1, F, method = "summed_regression")
#' lossFUN
#' lossFUN(dat)
SegmentLoss <- function(n_obs, lambda, penalize_diagonal, threshold = 1e-07,
                        method = c("glasso", "nodewise_regression", "summed_regression", "ratio_regression"),
                        ...) {
  args <- list(...)
  meth <- match.arg(method)

  if (meth == "nodewise_regression") {
    p <- args[["p"]]
    stopifnot(length(p) == 1 && is.numeric(p))
  }

  if (meth == "glasso") {
    function(x) {
      obs_count <- nrow(x)
      obs_share <- obs_count / n_obs

      # We need more than one observation to caclculate the covariance matrix
      stopifnot(obs_count > 1)

      n_p <- ncol(x)

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
  } else if (meth == "nodewise_regression") {
    function(x) {
      obs_count <- nrow(x)
      obs_share <- obs_count / n_obs

      # We need more than one observation
      stopifnot(obs_count > 1)

      s_dev_y <- sqrt((obs_count - 1) / obs_count * var(x[, p, drop = F]))
      fit <- glmnet::glmnet(
        x[, -p, drop = F], x[, p, drop = F],
        alpha = 1, lambda = lambda / sqrt(obs_share) * s_dev_y,
        thresh = threshold
      )
      deviance(fit) / n_obs
    }
  } else if (meth == "ratio_regression") {
    function(x) {
      obs_count <- nrow(x)
      obs_share <- obs_count / n_obs

      # We need more than one observation to caclculate the covariance matrix
      stopifnot(obs_count > 1)

      cov_mat <- (obs_count - 1) / obs_count * cov(x)

      withCallingHandlers({
        glasso_output <- glasso::glasso(
        cov_mat,
        rho = lambda / sqrt(obs_share) * diag(cov_mat),
        approx = T,
        thr = threshold
      )$wi
      }, warning = HandleGlassoNaN)

      mean_vec <- colMeans(x)
      intercepts <- mean_vec - colSums(glasso_output * mean_vec)
      ss <- t(t(x - x %*% glasso_output) - intercepts) ^ 2
      loss <- obs_count * log(colSums(ss) / obs_count)

      mean(loss / n_obs)
    }
  } else if (meth == "summed_regression") {
    function(x) {
      obs_count <- nrow(x)
      obs_share <- obs_count / n_obs

      # We need more than one observation to caclculate the covariance matrix
      stopifnot(obs_count > 1)

      cov_mat <- (obs_count - 1) / obs_count * cov(x)

     withCallingHandlers({
       glasso_output <- glasso::glasso(
        cov_mat,
        rho = lambda / sqrt(obs_share) * diag(cov_mat),
        approx = T,
        thr = threshold
      )$wi
     }, warning = HandleGlassoNaN)

      mean_vec <- colMeans(x)
      intercepts <- mean_vec - colSums(glasso_output * mean_vec)
      ss <- t(t(x - x %*% glasso_output) - intercepts) ^ 2
      loss <- colSums(ss)

      mean(loss / n_obs)
    }
  }
}
