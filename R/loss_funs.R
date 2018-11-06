#' SplitLoss
#'
#' Caluclates the sum of the loss that results from splitting the segment at the given point
#'
#' @inheritParams FindBestSplit
#' @param split_point Index on which to split the segment
#'
#' @return Sum of the loss for both new segments after split.
SplitLoss <- function(split_point, SegmentLossFUN, start, end) {
  SegmentLossFUN(start = start, end = split_point - 1) + SegmentLossFUN(start = split_point, end = end)
}



#' av_SegmentLoss
av_SegmentLoss <- function(x,
                           lambda,
                           penalize_diagonal = FALSE,
                           standardize = TRUE,
                           threshold = 1e-07){
  n_obs <- nrow(x)
  n_p <- ncol(x)

  function(start, end){

    for(i in 1:n_p){
      x[start : end, i][is.na(x[start : end, i]==T)] <- mean(x[start : end, i], na.rm=T)
    }

    obs_count <- end - start + 1
    obs_share <- obs_count / n_obs

    # We need more than one observation to caclculate the covariance matrix
    stopifnot(obs_count > 1)


    cov_mat <- (obs_count - 1) / obs_count * cov(x[start : end, ], use = 'pairwise')
    cov_mat[is.na(cov_mat)] <- 0 # adapt value
    cov_mat <- as.matrix(nearPD(cov_mat)$mat)



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
}

#' EM_SegmentLoss
EM_SegmentLoss <- function(x,
                           lambda,
                           method,
                           penalize_diagonal = FALSE,
                           standardize = TRUE,
                           threshold = 1e-7,
                           ...){

  n_obs <- nrow(x)
  n_p <- ncol(x)
  mis <- is.na(x)

  function(start, end){

    mis <- mis[start : end, ]
    x_est <- x[start : end, ]

    obs_count <- end - start + 1
    obs_share <- obs_count / n_obs

    # INIT: Do 1 cov & mean estimation using complete data
    cov_mat <- (obs_count - 1) / obs_count * cov(x_est, use = 'pairwise') #use = 'pairwise' works with missing data
    cov_mat[is.na(cov_mat)] <- 0 # adapt value
    cov_mat <- as.matrix(nearPD(cov_mat)$mat)

    mu <- apply(x_est, 2, function(y) mean(y, na.rm = T))

    #### Do a first glasso fit
    K <- glasso::glasso(cov_mat,
                        rho = lambda / sqrt(obs_share) * diag(cov_mat),
                        penalize.diagonal = penalize_diagonal,
                        thr = threshold
    )$wi

    #### Do covariance imputation on missing values
    # TODO: can this be done more elegantly / faster?

    for (i in 1:nrow(x_est)){

      mis_i <- mis[i,]

      if(all(mis_i)){
        x_est[i, ] <- mu
      } else if(sum(!mis_i) == 1){
        x_est_i <- x_est[i, ]
        x_est_i[mis_i] <- mu[mis_i] - solve(K[mis_i, mis_i], K[mis_i, which(!mis_i)] * (x_est_i[!mis_i] - mu[!mis_i]))
        x_est[i, ] <- x_est_i
      } else {
        x_est_i <- x_est[i, ]
        x_est_i[mis_i] <- mu[mis_i] - solve(K[mis_i, mis_i], K[mis_i, !mis_i] %*% (x_est_i[!mis_i] - mu[!mis_i]))
        x_est[i, ] <- x_est_i
      }
    }

    ##### Do an other glasso fit on imputed data
    cov_mat <- cov(x_est)
    mu <- apply(x_est, 2, function(y) mean(y, na.rm = T))

    glasso_output <- glasso::glasso(
      cov_mat,
      rho = lambda / sqrt(obs_share) * diag(cov_mat),
      penalize.diagonal = penalize_diagonal,
      thr = threshold
    )

    # prepare output
    if (!penalize_diagonal) {
      diag(glasso_output$wi) <- 0
    }

    ((glasso_output$loglik / (-n_p / 2) # Needed to undo transformation of likelihood in glasso package
      - lambda / sqrt(obs_share) * sum(abs(glasso_output$wi))) * obs_share) # Remove regularizer added in glasso package
  }
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
SegmentLoss <- function(x,
                        y = NULL,
                        lambda = 0,
                        alpha = 1,
                        penalize_diagonal = FALSE,
                        standardize = TRUE,
                        threshold = 1e-07,
                        method = c("nodewise_regression", "summed_regression", "ratio_regression", "glasso", "elastic_net"),
                        ...) {
  args <- list(...)
  mth <- match.arg(method)
  n_obs <- NROW(x)

  if (mth == "nodewise_regression") {
    p <- args[["node"]]
    stopifnot(length(p) == 1 && is.numeric(p))
  }

  if (mth == 'elastic_net'){
    stopifnot(length(y) == nrow(x))
  }

  if (mth == "glasso") {
    function(start, end, ...) {
      obs_count <- end - start + 1
      obs_share <- obs_count / n_obs

      # We need more than one observation to calculate the covariance matrix
      stopifnot(obs_count > 1)

      n_p <- NCOL(x)

      cov_mat <- (obs_count - 1) / obs_count * cov(x[start : end, ])

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
    function(start, end, ...) {
      obs_count <- end - start + 1
      obs_share <- obs_count / n_obs

      # We need more than one observation
      stopifnot(obs_count > 1)

      if (standardize) {
        s_dev_y <- sqrt((obs_count - 1) / obs_count * var(x[start : end, p, drop = F]))
      } else {
        s_dev_y <- 1
      }

      fit <- glmnet::glmnet(
        x[start : end, -p, drop = F], x[start : end, p, drop = F],
        alpha = 1, lambda = lambda / sqrt(obs_share) * s_dev_y,
        thresh = threshold,
        standardize = FALSE
      )
      deviance(fit) / n_obs
    }
  } else if (mth == "ratio_regression") {
    function(start, end, ...) {
      obs_count <- end - start + 1
      obs_share <- obs_count / n_obs

      # We need more than one observation to caclculate the covariance matrix
      stopifnot(obs_count > 1)

      cov_mat <- (obs_count - 1) / obs_count * cov(x[start : end, ])

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

      mean_vec <- colMeans(x[start : end, ])
      intercepts <- mean_vec - colSums(glasso_output * mean_vec)
      ss <- t(t(x[start : end, ] - x[start : end, ] %*% glasso_output) - intercepts)^2
      loss <- obs_count * log(colSums(ss) / obs_count)

      mean(loss / n_obs)
    }
  } else if (mth == "summed_regression") {
    function(start, end, ...) {
      obs_count <- end - start + 1
      obs_share <- obs_count / n_obs

      # We need more than one observation to caclculate the covariance matrix
      stopifnot(obs_count > 1)

      cov_mat <- (obs_count - 1) / obs_count * cov(x[start : end, ])

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

      mean_vec <- colMeans(x[start : end, ])
      intercepts <- mean_vec - colSums(glasso_output * mean_vec)
      ss <- t(t(x[start : end, ] - x[start : end, ] %*% glasso_output) - intercepts)^2
      loss <- colSums(ss)

      mean(loss / n_obs)
    }
  } else if (mth == "elastic_net") {
    if (is.null(args$family)){
      family <- 'gaussian'
    } else {
      family <-  args$family
    }

    function(start, end, ...){

      obs_count <- end - start + 1
      obs_share <- obs_count / n_obs
      fit <- glmnet::glmnet(x[start : end, ], y[start : end], alpha = alpha, lambda = sqrt(obs_share) * lambda, standardize = standardize, family = family, thres = threshold)
      deviance(fit)

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

  function(start, end) {
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

  function(start, end) {
    stopifnot(end >= start && end <= n_obs && start >= 1)

    sum((x[start : end] - mean(x[start : end]))^2) / n_obs
  }
}
