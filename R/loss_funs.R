#' SplitLoss
#'
#' Caluclates the sum of the loss that results from splitting the segment at the given point
#'
#' @inheritParams FindBestSplit
#' @param split_point Index on which to split the segment
#'
#' @return Sum of the loss for both new segments after split.
SplitLoss <- function(split_point, SegmentLossFUN, start, end, lambda) {
  SegmentLossFUN(start = start, end = split_point - 1, lambda = lambda) + SegmentLossFUN(start = split_point, end = end, lambda = lambda)
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
                        lambda = NULL,
                        NA_method =c('complete_observations', 'pairwise_covariance',
                                     'loh_wainwright_bias_correction', 'average_imputation',
                                     'expectation_maximisation'),
                        method = c("nodewise_regression", "summed_regression", "ratio_regression", "glasso", "elastic_net"),
                        control = NULL,
                        ...) {
  args <- list(...)
  NA_mth <- match.arg(NA_method)
  mth <- match.arg(method)
  n_obs <- nrow(x)
  p <- NCOL(x)
  lambda_global <- lambda

  if (mth == "glasso") {
    # load parameters for glasso from control
    penalize_diagonal <-  control$glasso_penalize_diagonal
    standardize <-  control$glasso_standardize
    threshold <- control$glasso_threshold
    min_points <- control$segment_loss_min_points

    function(start, end, lambda, ...) {

      if(!is.null(lambda_global)){
        lambda <- lambda_global
      }

      cov_mat_output <- get_cov_mat(x[start  : end , , drop = F], NA_method, min_points = min_points)
      cov_mat <- cov_mat_output$mat
      inds <- cov_mat_output$inds
      p_cur <- sum(inds) # which predictors have sufficient non-missing values

      if (p_cur == 0){
        return(NA)
      }

      obs_share <- cov_mat_output$n_eff_obs / n_obs # for rescaling of lambda

      if (standardize){
        glasso_output <- glasso::glasso(
          cov_mat,
          rho = lambda / sqrt(obs_share) * diag(cov_mat),
          penalize.diagonal = penalize_diagonal,
          thr = threshold
        )
      } else {
        glasso_output <- glasso::glasso(
          cov_mat,
          rho = lambda / sqrt(obs_share),
          penalize.diagonal = penalize_diagonal,
          thr = threshold
        )
      }


      if(T | any(is.na(x[start : end, inds])) & NA_mth != 'complete_observations'){ # TODO ADAPT WHY DOES THIS WORK
        loglikelihood(x[start : end, cov_mat_output$inds], colMeans(x[start : end, cov_mat_output$inds], na.rm = T), glasso_output$w, glasso_output$wi) / n_obs
      } else {
        # The following is to compute the test error using the glasso output
        if (!penalize_diagonal){
          diag(glasso_output$wi) <- 0
        }
        # Returns the loglikelihood divided by n_obs. Needed to undo transformation done by glasso package
        if (!standardize){
          ((-2 / p_cur) * glasso_output$l - sum(abs(lambda / sqrt(obs_share) * glasso_output$wi))) * obs_share
        } else {
          ((-2 / p_cur) * glasso_output$l - sum(abs(lambda / sqrt(obs_share) * sqrt((diag(cov_mat) %*% t(diag(cov_mat)))) * glasso_output$wi))) * obs_share
        }
      }
    }

  } else if (mth == "nodewise_regression") {

    node = control$nodewise_regression_node
    stopifnot(length(node) == 1 && is.numeric(node))

    function(start, end, ...) {
      n_cur <- end - start + 1
      obs_share <- n_cur / n_obs

      # We need more than one observation
      stopifnot(n_cur > 1)

      if (standardize) {
        s_dev_y <- sqrt((n_cur - 1) / n_cur * var(x[start : end, node, drop = F]))
      } else {
        s_dev_y <- 1
      }
      fit <- glmnet::glmnet(
        x[start : end, -node, drop = F], x[start : end, node, drop = F],
        alpha = 1, lambda = lambda / sqrt(obs_share) * s_dev_y,
        thresh = threshold,
        standardize = FALSE
      )
      deviance(fit) / n_obs
    }

  } else if (mth == "ratio_regression") {
    function(start, end, ...) {
      n_cur <- end - start + 1

      stopifnot(n_cur > 1)  # We need more than one observation to calculate the covariance matrix

      cov_mat_output <- get_cov_mat(x[start : end, , drop = F], NA_mth)

      obs_share <- cov_mat_output$n_eff_obs / n_obs

      if (standardize) {
        var_x <- diag(cov_mat_output$mat)
      } else {
        var_x <- 1
      }

      withCallingHandlers({
        glasso_output <- glasso::glasso(
          cov_mat_output$mat,
          rho = lambda / sqrt(obs_share) * var_x,
          approx = T,
          thr = threshold
        )$wi
      }, warning = HandleGlassoNaN)

      inds <- cov_mat_output$inds

      # Verstehe ich nicht
      mean_vec <- colMeans(x[start : end, inds, drop = F])
      intercepts <- mean_vec - colSums(glasso_output * mean_vec)
      ss <- t(t(x[start : end, inds, drop = F] - x[start : end, inds, drop = F] %*% glasso_output) - intercepts)^2
      loss <- n_cur * log(colSums(ss) / n_cur)

      mean(loss / n_obs)
    }

  } else if (mth == "summed_regression") {
    function(start, end, ...) {
      n_cur <- end - start + 1

      stopifnot(n_cur > 1) # We need more than one observation to calculate the covariance matrix

      cov_mat_output <- get_cov_mat(x[start : end, , drop = F], NA_method)
      obs_share = cov_mat_output$n_eff_obs / n_obs

      if (standardize) {
        var_x <- diag(cov_mat_output$mat)
      } else {
        var_x <- 1
      }

      withCallingHandlers({
        glasso_output <- glasso::glasso(
          cov_mat_output$mat,
          rho = lambda / sqrt(obs_share) * var_x,
          approx = T,
          thr = threshold
        )$wi
      }, warning = HandleGlassoNaN)

      inds <- cov_mat_output$inds

      # Verstehe ich nicht. Ist dies mit inds legitim?
      mean_vec <- colMeans(x[start : end, inds, drop = F])
      intercepts <- mean_vec - colSums(glasso_output * mean_vec)
      ss <- t(t(x[start : end, inds, drop = F] - x[start : end, inds, drop = F ] %*% glasso_output) - intercepts)^2
      loss <- colSums(ss)

      mean(loss / n_obs)
    }

  } else if (mth == "elastic_net") {

    alpha = control$elastic_net_alpha
    stopifnot(is.numeric(alpha) && length(alpha) == 1 && 0 <= alpha && alpha <= 1)
    family <- control$elastic_net_family

    function(start, end, ...){
      n_cur <- end - start + 1
      obs_share <- n_cur / n_obs
      fit <- glmnet::glmnet(x[start : end, -1], x[start : end, 1], alpha = alpha, lambda = sqrt(obs_share) * lambda, standardize = standardize, family = family, thres = threshold)
      deviance(fit) / n_obs
    }
  }
}


#########
# Arguments
# x: Matrix
# NA_method: string
# evaluate: boolean
#
# Returns:
# estimate of covariance matrix
get_cov_mat <- function(x, NA_method = c('complete_observations', 'pairwise_covariance',
                                        'loh_wainwright_bias_correction', 'average_imputation',
                                        'expectation_maximisation'), min_points = 2){
  NA_mth <- match.arg(NA_method)
  if(NA_mth == "complete_observations"){
    n_obs <- sum(complete.cases(x))
    list(mat = (n_obs - 1) / n_obs * cov(x, use = 'na.or.complete'), inds = rep(T, ncol(x)), n_eff_obs = n_obs)
  } else {
    available_obs <- apply(!is.na(x), 2, sum)
    inds <- (available_obs >= max(2, min_points))
    cov_mat <- calc_cov_mat(x[, inds], NA_mth)
    list(mat = cov_mat, inds = inds, n_eff_obs = sum(available_obs[inds]) / sum(inds))
  }
}

# calculates the covariance method of x
calc_cov_mat <- function(x, NA_mth){
 if (NA_mth == "pairwise_covariance"){ # calculate pairwise covariances
    cov_mat <- cov(x, use = 'pairwise')
    cov_mat[is.na(cov_mat)] <- 0 # for some pairs of features there might not be two or more available observations
    as.matrix(Matrix::nearPD(cov_mat)$mat) # result might not be posd, which is required for e.g. glasso
  } else if (NA_mth == "loh_wainwright_bias_correction"){
    z <- scale(x, T, F) # center
    z[is.na(z)] <- 0 # impute mean
    miss_frac <- apply(is.na(x), 2, mean) # estimate missingness per predictor
    stopifnot(length(miss_frac) == ncol(z))
    M <- matrix(rep((1 - miss_frac), ncol(z)), ncol = ncol(z))
    M <- t(M) * M
    diag(M) <- 1 - miss_frac
    as.matrix(Matrix::nearPD(cov(z) * (nrow(z) - 1) / nrow(z) / M)$mat)
    #z <- z %*% diag(1 / (1 - miss_frac), nrow = length(miss_frac)) # first step in the loh-wainwright bias correction
    #cov_mat <- cov(z)
    #cov_mat <- cov_mat - diag(miss_frac * diag(cov_mat), nrow = length(miss_frac)) # second step loh-wainwright correction
    #as.matrix(Matrix::nearPD(cov_mat)$mat)
  } else if (NA_mth == "average_imputation"){
    z <- scale(x, T, F) # center
    z[is.na(z)] <- 0 #impute mean
    cov(z) * (nrow(z) - 1) / nrow(z)
  }
}

# Takes an input with missing values and imputes values using the max likelihood estimator
# for a gaussian with mean mu and covariance cov_mat


# n <- length(x)
# k <- sum(is.na(x))
# A <- diag(n)[, is.na(x)]
# U <- cbind(A, cov_mat[, is.na(x)])
# U[, (k + 1) : (2*k)][as.logical(A)] <- 0
# V <- rbind(cov_mat[is.na(x), ], t(A))
# V[1 : k, ][as.logical(t(A))] <- 0
# V[1 : k, ][as.logical(t(A[, k : 1]))] <- 0
# cov_mat_inv_red <- cov_mat_inv - cov_mat_inv %*% U %*% solve(-diag(2 * k) + V %*% cov_mat_inv %*% U) %*% V %*% cov_mat_inv


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
