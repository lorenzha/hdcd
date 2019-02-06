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

#' cv_loss
#'
#' This function returns the loss after fitting a glasso or elastic-net fit on \emph{x_train} and evaluating on \emph{x_test}
#'
#' If data with missing values is supplied then a first estimate of the covariance matrix is estimated with the \emph{NA_method}.
#' Then the fit is improved by doing \emph{nrep_EM} Expectation-Maximisation fits.
#'
#' @inheritParams BinarySegmentation
#' @param x_train Training observations
#' @param x_test Test observations
#' @n_obs_train Total amount of training observations in the whole model
#'
#' @return list of train error, best lambda and possibly test error
cv_loss <- function(x_train, n_obs,
                    x_test = NULL,
                    method = c("nodewise_regression", "summed_regression", "ratio_regression", "glasso", "elastic_net"),
                    NA_method = c("complete_observations", "average_imputation", "pairwise_covariance_estimation", 'loh_wainwright_bias_correction'),
                    control = NULL){
  # TODO: after implementing control lambda should be passed as the current lambda value
  # such that for n_fold_inner the outer lambda is chosen.
  n_cur_train <- nrow(x_train)
  n_folds_inner <- control_get(control, "n_folds_inner", 4)
  randomize_inner <- control_get(control, "randomize_inner", FALSE)
  lambda_inner <- control_get(control, "lambda_inner", 1) #ADAPT!!!

  folds_inner <- sample_folds(n_cur_train, n_folds_inner, randomize_inner)

  # to store information about cv
  loss <- list(cv = numeric(length(lambda_inner)))

  for (i in seq_along(lambda_inner)){
    for (inner in 1:n_folds_inner){

      x_train_train <- x_train[folds_inner != inner, , drop = F]
      x_train_test <- x_train[folds_inner == inner, , drop = F]

      # should we divide by something? if folds are not of equal size
      loss[["cv"]][i] <- loss[["cv"]][i] + cv_fit(x_train_train, x_train_test, lambda = lambda_inner[i], method = method, NA_method = NA_method, n_obs_train = n_obs, control = control)

    }
  }

  loss[["lambda_opt"]] <- lambda_inner[which.min(loss[["cv"]])]
  loss[["train_loss"]] <- min(loss[["cv"]])

  if (!is.null(x_test)){
    loss[["test_loss"]] <- cv_fit(x_train, x_test, lambda = loss[["lambda_opt"]], method = method, NA_method = NA_method, n_obs_train = n_obs_train, control = control)
  }
  loss
}


cv_fit <- function(x_train, x_test, lambda, method, NA_method, n_obs_train, control = list()){

  if(nrow(x_test) == 0){
    return(0)
  }

  mth <- match.arg(method, choices = c("nodewise_regression", "summed_regression", "ratio_regression", 'glasso', 'elastic_net'))
  NA_mth <- match.arg(NA_method, choices = c('complete_observations', 'pairwise_covariance_estimation',
                                             'loh_wainwright_bias_correction', 'average_imputation',
                                             'expectation_maximisation'))

  penalize_diagonal <- control_get(control, "penalize_diagonal", FALSE)
  standardize <- control_get(control, "standardize", TRUE)
  threshold <- control_get(control, "threshold", 1e-7)
  nrep_em <- control_get(control, "nrep_em", 5)

  n_cur_train <- nrow(x_train)
  obs_share_train <- n_cur_train / n_obs_train

  if (mth %in% c("nodewise_regression", "summed_regression", "glasso", "ratio_regression")){

    cov_mat_output <- get_cov_mat(x_train, NA_method = NA_method)
    cov_mat <- cov_mat_output$mat
    inds <- cov_mat_output$inds

    if(any(is.na(cov_mat))){
      warning('Not enough observations in training set to fit, loss 0 returned')
      return(0)
    }

    mu <- colMeans(x_train[, inds, drop = F], na.rm = T)

    glasso_output <- glasso::glasso(
      cov_mat,
      rho = lambda / sqrt(obs_share_train) * diag(cov_mat),
      penalize.diagonal = penalize_diagonal,
      thr = threshold
    )

    cov_mat <- glasso_output$w
    cov_mat_inv <- glasso_output$wi

    if (NA_mth %in% c('average_imputation', 'pairwise_covariance_estimation', 'loh_wainwright_bias_correction') & any(is.na(x_train))){

      for (i in 1:nrep_em) {
        x_impute <- apply(x_train[, inds, drop = F], 1, impute_em, mu = mu, cov_mat = cov_mat)
        cov_mat <- cov(x_impute)
        glasso_output <- glasso::glasso(
          cov_mat,
          rho = lambda / sqrt(obs_share_train) * diag(cov_mat),
          penalize.diagonal = penalize_diagonal,
          thr = threshold
        )
        cov_mat <- glasso_output$w
        cov_mat_inv <- glasso_output$wi
      }
      rm(x_impute)
      rm(glasso_output)
    }

    loglikelihood(x_test, mu, cov_mat, cov_mat_inv)

  } else if (mth == 'elastic_net'){

    family <-  args[['family']]
    if(is.null(family)) family <- 'gaussian'

    alpha <- args[['alpha']]
    if (is.null(args[['alpha']])) alpha <- 0

    glmnet_output <- glmnet(x = x_train[, -1], y = x_train[, 1],
                            lambda = lambda / sqrt(obs_share_train), alpha = alpha, family = family)

    y_pred <- predict(glmnet_output, x_test[, -1])

    mean( (x_test[, 1] - y_pred)^2 )
  }
}




# calculates loglikelihood of x (with missing values) given mu, cov_mat. cov_mat_inv can be used to speed up calculaitons
# of the inverse of the submatrix via the Woodbury matrix identity as a low rank update from cov_mat_inv
# Maybe improve speed??
loglikelihood <- function(x, mu, cov_mat, cov_mat_inv){
  p <- ncol(x)
  inds <- is.na(x) # missingness structure of x
  na_order <- do.call(order, lapply(1:ncol(inds), function(i) inds[, i])) # order x by missingness structure
  cov_mat_inv_cur <- cov_mat_inv
  loss <- 0

  inds_old <- rep(F, ncol(x))

  for(i in 1:nrow(x)){
    inds_cur <- inds[na_order[i], ]
    if(any(inds_cur != inds_old)){
      # fast calculation of cov_mat_inv_cur as update from cov_mat_inv via Woodbury matrix identity
      k <- sum(inds[na_order[i], ]) # amount of missing values
      A <- diag(p)[, inds_cur, drop = F] # helper matrix
      U <- cbind(A, cov_mat[, inds_cur])
      U[, (k + 1) : (2*k)][as.logical(A)] <- 0
      V <- rbind(cov_mat[inds_cur, ], t(A))
      V[1 : k, ][as.logical(t(A))] <- 0
      V[1 : k, ][as.logical(t(A[, k : 1]))] <- 0
      cov_mat_inv_cur <- (cov_mat_inv - cov_mat_inv %*% U %*% solve(-diag(2 * k) + V %*% cov_mat_inv %*% U) %*% V %*% cov_mat_inv)[!inds_cur, !inds_cur]
    }
    v <- (x[na_order[i], ] - mu)[!inds_cur]
    loss <- loss + t(v) %*% cov_mat_inv_cur %*% v - log(abs(det(cov_mat_inv_cur)))
  }
  loss
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
                        lambda = 0,
                        NA_method = "complete_observations",
                        method = c("nodewise_regression", "summed_regression", "ratio_regression", "glasso", "elastic_net"),
                        control = NULL,
                        ...) {
  args <- list(...)
  mth <- match.arg(method)
  n_obs <- nrow(x)
  p <- NCOL(x)

  if (mth == "glasso") {
    # load parameters for glasso from control
    penalize_diagonal <-  control_get(control, "penalize_diagonal", FALSE)
    standardize <-  control_get(control, "standardize", TRUE)
    threshold <-  control_get(control, "threshold", 1e-07)

    function(start, end, ...) {

      cov_mat_output <- get_cov_mat(x[start  : end , , drop = F], NA_method)
      cov_mat <- cov_mat_output$mat
      p_cur <- sum(cov_mat_output$inds) # which predictors have sufficient non-missing values

      if (p_cur == 0){
        return(NA)
      }

      obs_share <- cov_mat_output$n_eff_obs / n_obs # for rescaling of lambda

      glasso_output <- glasso::glasso( # correction with 1 / sqrt(obs_share) from asymptotic theory
        cov_mat,
        rho = lambda / sqrt(obs_share) * diag(cov_mat),
        penalize.diagonal = penalize_diagonal,
        thr = threshold
      )

      # The following is to compute the test error using the glasso output
      if (!penalize_diagonal){
        diag(glasso_output$wi) <- 0
      }
      # Needed to undo transformation of likelihood in glasso package
      # equivalent to ' - log(det(glasso_output$wi)) + psych::tr(cov_mat %*% glasso_output$wi)' before setting diag() <- 0
      #- 2 * glasso_output$loglik / n_p
      #      - lambda  / sqrt(obs_share) * as.numeric(sqrt(diag(cov_mat)) %*% abs(glasso_output$wi) %*% sqrt(diag(cov_mat) ))
      # MESS!!! TODO: cleanup!
      p / p_cur * (((glasso_output$loglik / (-p_cur / 2) # Needed to undo transformation of likelihood in glasso package
        - lambda / sqrt(obs_share) * sum(abs(glasso_output$wi))) * obs_share)) # Remove regularizer added in glasso package
      #n_p * 2 * -glasso_output$loglik - (n_p / cur_p) * lambda * log(cur_p) / log(n_p) * sum(abs(glasso_output$wi)) * sqrt(obs_share)
    }

  } else if (mth == "nodewise_regression") {

    node = control_get(control, "node", 1)
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

    alpha = control_get(control, "alpha", 1)
    stopifnot(is.numeric(alpha) && length(alpha) == 1 && 0 <= alpha && alpha <= 1)
    family <- control_get(control, "family", "gaussian")

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
get_cov_mat <- function(x, NA_method = c('complete_observations', 'pairwise_covariance_estimation',
                                        'loh_wainwright_bias_correction', 'average_imputation',
                                        'expectation_maximisation')){
  n_obs <- nrow(x)
  NA_mth <- match.arg(NA_method)
  if(NA_mth == "complete_observations"){
    list(mat = cov(x, use = 'na.or.complete'), inds = rep(T, ncol(x)), n_eff_obs = n_obs)
  } else {
    available_obs <- n_obs - apply(is.na(x), 2, sum)
    inds <- (available_obs >= 2)
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
    z <- z %*% diag(1 / (1 - miss_frac), nrow = length(miss_frac)) # first step in the loh-wainwright bias correction
    cov_mat <- cov(z)
    cov_mat <- cov_mat - diag(miss_frac * diag(cov_mat), nrow = length(miss_frac)) # second step loh-wainwright correction
    as.matrix(Matrix::nearPD(cov_mat)$mat)
  } else if (NA_mth == "average_imputation"){
    z <- scale(x, T, F) # center
    z[is.na(z)] <- 0 #impute mean
    cov(z)
  }
}

# Takes an input with missing values and imputes values using the max likelihood estimator
# for a gaussian with mean mu and covariance cov_mat
impute_em <- function(x, mu, cov_mat){
  if (!any(is.na(x))) {
    x
  } else if (all(is.na(x))) {
    mu
  } else {
  x[is.na(x)] <- mu[is.na(x)] - solve(cov_mat[is.na(x), is.na(x), drop = F],
                                        cov_mat[is.na(x), !is.na(x), drop = F] %*% (x[!is.na(x), drop = F] - mu[!is.na(x), drop = F])
    )
  x
  }
}


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
