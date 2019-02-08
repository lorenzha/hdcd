
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

  stopifnot(length(n_obs) == 1)

  n_cur_train <- nrow(x_train)
  n_folds_inner <- control_get(control, "n_folds_inner", 4)
  randomize_inner <- control_get(control, "randomize_inner", FALSE)
  lambda_inner <- control_get(control, "lambda_inner", 1) ### ADAPT THIS

  folds_inner <- sample_folds(n_cur_train, n_folds_inner, randomize_inner)

  # to store information about cv
  loss <- numeric(length(lambda_inner))

  for (inner in 1:n_folds_inner){
    x_train_train <- x_train[folds_inner != inner, , drop = F]
    x_train_test <- x_train[folds_inner == inner, , drop = F]

    loss <- loss + cv_fit(x_train_train, x_train_test, lambda = lambda_inner, method = method, NA_method = NA_method, n_obs = n_obs, control = control)
  }

  lambda_opt <- lambda_inner[which.min(loss)]
  train_loss <- min(loss)

  if (!is.null(x_test)){
    test_loss <- cv_fit(x_train, x_test, lambda = lambda_opt, method = method, NA_method = NA_method, n_obs = n_obs, control = control)
    list(loss = loss, lambda_opt = lambda_opt, train_loss = train_loss, test_loss = test_loss)
  } else {
    list(loss = loss, lambda_opt = lambda_opt, train_loss = train_loss)
  }
}


cv_fit <- function(x_train, x_test, lambda, method, NA_method, n_obs, control = list()){

  if(nrow(x_test) == 0){
    return(0)
  }

  mth <- match.arg(method, choices = c("nodewise_regression", "summed_regression", "ratio_regression", 'glasso', 'elastic_net'))
  NA_mth <- match.arg(NA_method, choices = c('complete_observations', 'pairwise_covariance_estimation',
                                             'loh_wainwright_bias_correction', 'average_imputation',
                                             'expectation_maximisation'))

  penalize_diagonal <- control_get(control, "penalize_diagonal", FALSE)
  standardize <- control_get(control, "standardize", TRUE)
  threshold <- control_get(control, "threshold", 1e-4)
  nrep_em <- control_get(control, "nrep_em", 5)

  n_cur_train <- nrow(x_train)
  obs_share_train <- n_cur_train / n_obs

  if (mth %in% c("nodewise_regression", "summed_regression", "glasso", "ratio_regression")){

    cov_mat_output <- get_cov_mat(x_train, NA_method = NA_method)
    inds <- cov_mat_output$inds
    cov_mat <- cov_mat_output$mat

    if(any(is.na(cov_mat))){
      warning('Not enough observations in training set to fit, loss 0 returned')
      return(0)
    }

    mu <- colMeans(x_train[, inds, drop = F], na.rm = T)

    if (standardize == TRUE){ #apparently glassopath doesn't do standardisation
      cov_mat <- cov_mat_inv <- array(NA, dim = c(dim(cov_mat_output$mat), length(lambda)))
      for(j in seq_along(lambda)){
        glasso_output <- glasso::glasso(
          cov_mat_output$mat,
          rho = lambda[j] / sqrt(obs_share_train) * diag(cov_mat_output$mat),
          penalize.diagonal = penalize_diagonal,
          thr = threshold,
          trace = FALSE
        )
        cov_mat[, , j] <- glasso_output$w
        cov_mat_inv[, , j] <- glasso_output$wi
      }
    } else {
      glasso_output <- glasso::glassopath(
        cov_mat,
        rho = lambda / sqrt(obs_share_train) * mean(diag(cov_mat)),
        penalize.diagonal = penalize_diagonal,
        thr = threshold,
        trace = FALSE
      )
      cov_mat <- glasso_output$w
      cov_mat_inv <- glasso_output$wi
    }


  if (NA_mth %in% c('average_imputation', 'pairwise_covariance_estimation', 'loh_wainwright_bias_correction') & any(is.na(x_train))){

    for (j in seq_along(lambda)){
      for (i in seq_len(nrep_em)) {
        x_impute <- t(apply(x_train[, inds, drop = F], 1, impute_em, mu = mu, cov_mat = cov_mat[, , j]))
        cov_mat[, , j] <- cov(x_impute)
        glasso_output <- glasso::glasso(
          cov_mat[, , j],
          rho = lambda[j] / sqrt(obs_share_train) * mean(diag(cov_mat[, , j])),
          penalize.diagonal = penalize_diagonal,
          thr = threshold
        )
        cov_mat[, , j] <- glasso_output$w
        cov_mat_inv[, , j] <- glasso_output$wi
      }
    }
  }

  loss <- numeric(length(lambda))

  for (i in seq_along(lambda)){
    loss[i] <- loglikelihood(x_test, mu, cov_mat[,,i], cov_mat_inv[,,i])
  }

  loss

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
  log_det_cur <- log(det(cov_mat_inv_cur))
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
      log_det_cur <- log(det(cov_mat_inv_cur))
    }
    v <- (x[na_order[i], ] - mu)[!inds_cur]
    loss <- loss + t(v) %*% cov_mat_inv_cur %*% v - log_det_cur
  }
  loss
}

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
