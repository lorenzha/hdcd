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
                      folds_inner,
                      x_test = NULL,
                      lambda = NULL,
                      lambda0 = NULL,
                    method = NULL, #POssibly adapt so this works with elastic net
                      NA_method = c("complete_observations", "average_imputation", "pairwise_covariance_estimation", 'loh_wainwright_bias_correction'),
                      control = NULL){

  n_folds_inner <- control$cv_inner_n_folds
  standardize <- control$glasso_standardize
  search_lambda_inner <- control$cv_inner_search_lambda
  lambda_inner_step <- control$cv_inner_lambda_step
  threshold <- control$glasso_threshold
  penalize_diagonal <- control$glasso_penalize_diagonal
  lambda <- control$cv_inner_lambda

  if(lambda_inner_step < 1){
    lambda_inner_step <- 1 / lambda_inner_step
  }

  stopifnot(length(n_obs) == 1)

  x_train_list <- lapply(1 : n_folds_inner, function(y) x_train[folds_inner != y, ])
  cov_mat_output_list <- lapply(x_train_list, get_cov_mat, NA_method)
  inds_list <- lapply(cov_mat_output_list, '[[', 'inds')
  means_list <- mapply(function(x, j) list(colMeans(x[, j], na.rm = T)), x_train_list, inds_list)

  x_test_list <- mapply( function(i, j) list(x_train[folds_inner == i, j]), as.list(1 : n_folds_inner), inds_list)

  f <- function(lambda){
    if(standardize){
      rho <- mapply(function(cov_mat_output, lambda){
        list(lambda * cov_mat_output$n_eff_obs / n_obs * diag(cov_mat_output$mat))
        },
        cov_mat_output_list, as.list(c(sapply(lambda, rep, n_folds_inner))))
    } else {
      rho <- mapply(function(cov_mat_output, lambda){
        list(lambda * cov_mat_output$n_eff_obs / n_obs)
      },
      cov_mat_output_list, as.list(repc(sapply(lambda, rep, n_folds_inner))))
    }

    cov_mat_list <- lapply(cov_mat_output_list, '[[', 'mat')

    glasso_output_list <- apply(mapply(glasso::glasso,
                                 cov_mat_list,
                                 rho,
                                 MoreArgs = list(thr = threshold, penalize.diagonal = penalize_diagonal)
    ), 2, as.list)

    matrix(mapply(loglikelihood, x_test_list, means_list,
                  lapply(glasso_output_list, '[[', 'w'),
                  lapply(glasso_output_list, '[[', 'wi')),
           nrow = n_folds_inner)
  }


  if(search_lambda_inner){

    lambda <- c(1 / lambda_inner_step, 1, lambda_inner_step) * lambda0

    loss <- f(lambda)
    loss_sum <- apply(loss, 2, sum)
    loss_sd <- apply(loss, 2, sd)
    i <- which.min(loss_sum)
    if(i == 2){
      list(loss_array = loss[, i], lambda_opt = lambda0, train_loss = loss_sum[i])
    } else if (i == 1){
      while(TRUE){
        lambda <- c(lambda[1] / lambda_inner_step, lambda)
        loss_cur <- f(lambda[1])
        loss_sum <- c(sum(loss_cur), loss_sum)
        loss_sd <- c(sd(loss_cur), loss_sd)
        loss <- cbind(loss_cur, loss)
        if(loss_sum[1] >= loss_sum[2]){
          return(list(loss_array = loss[, 2], lambda_opt = lambda[2], train_loss = loss_sum[2]))
        }
      }
    } else {
      while (TRUE){
        k <- length(lambda)
        lambda <- c(lambda, lambda[k] * lambda_inner_step)
        loss_cur <- f(lambda[k + 1])
        loss_sum <- c(loss_sum, sum(loss_cur))
        loss_sd <- c(loss_sd, sd(loss_cur))
        loss <- cbind(loss, loss_cur)
        if(loss_sum[k + 1] >= loss_sum[k]){
          return(list(loss_array = loss[, k], lambda_opt = lambda[k], train_loss = loss_sum[k]))
        }
      }
    }
  } else {
    loss <- f(lambda)
    loss_sum <- apply(loss, 2, sum)
    loss_sd <- apply(loss, 2, sd)
    i <- which.min(loss_sum)
    return(list(loss_array = loss[, i], lambda_opt = lambda[i], train_loss = loss_sum[i]))
  }
}
#
# cv_loss <- function(x_train, n_obs,
#                     folds_inner,
#                     x_test = NULL,
#                     method = c("nodewise_regression", "summed_regression", "ratio_regression", "glasso", "elastic_net"),
#                     NA_method = c("complete_observations", "average_imputation", "pairwise_covariance_estimation", 'loh_wainwright_bias_correction'),
#                     control = NULL){
#   # TODO: after implementing control lambda should be passed as the current lambda value
#   # such that for n_fold_inner the outer lambda is chosen.
#
#   n_folds_inner <- control_get(control, "n_folds_inner", 4)
#   lambda_inner <- control_get(control, "lambda_inner", 1)
#   nrep_em <- control_get(control, "nrep_em", 5)
#
#   stopifnot(length(n_obs) == 1)
#
#   n_cur_train <- nrow(x_train)
#
#   loss <- sapply(1 : n_folds_inner, function(y){
#     cv_fit(x_train[folds_inner != y, , drop = F],
#            x_train[folds_inner == y, , drop = F],
#            lambda = lambda_inner, method = method, NA_method = NA_method, n_obs = n_obs, control = control)
#   })
#
#
#
#   # to store information about cv
#   #loss <- array(NA, dim = c(length(lambda_inner), n_folds_inner))
#
#   # for (inner in 1 : n_folds_inner){
#   #   x_train_train <- x_train[folds_inner != inner, , drop = F]
#   #   x_train_test <- x_train[folds_inner == inner, , drop = F]
#   #
#   #   loss[, inner] <- cv_fit(x_train_train, x_train_test, lambda = lambda_inner, method = method, NA_method = NA_method, n_obs = n_obs, control = control)
#   # }
#
#   loss_sum <- apply(loss, 1, sum)
#   i <- which.min(loss_sum)
#   lambda_opt <- lambda_inner[i]
#   train_loss <- loss_sum[i]
#
#   if (!is.null(x_test)){
#     test_loss <- cv_fit(x_train, x_test, lambda = lambda_opt, method = method, NA_method = NA_method, n_obs = n_obs, control = control)
#     list(loss_array = loss[i, ], lambda_opt = lambda_opt, train_loss = train_loss, test_loss = test_loss)
#   } else {
#     list(loss_array = loss[i, ], lambda_opt = lambda_opt, train_loss = train_loss)
#   }
# }
#
#
# cv_fit <- function(x_train, x_test, lambda, method, NA_method, n_obs, control = list()){
#
#   if(nrow(x_test) == 0){
#     return(0)
#   }
#
#   mth <- match.arg(method, choices = c("nodewise_regression", "summed_regression", "ratio_regression", 'glasso', 'elastic_net'))
#   NA_mth <- match.arg(NA_method, choices = c('complete_observations', 'pairwise_covariance_estimation',
#                                              'loh_wainwright_bias_correction', 'average_imputation',
#                                              'expectation_maximisation'))
#
#   penalize_diagonal <- control_get(control, "penalize_diagonal", FALSE)
#   standardize <- control_get(control, "standardize", TRUE)
#   threshold <- control_get(control, "threshold", 1e-4)
#   nrep_em <- control_get(control, "nrep_em", 5)
#
#   n_cur_train <- nrow(x_train)
#   obs_share_train <- n_cur_train / n_obs
#
#   loss <- numeric(length(lambda))
#
#   if (mth %in% c("nodewise_regression", "summed_regression", "glasso", "ratio_regression")){
#
#     cov_mat_output <- get_cov_mat(x_train, NA_method = NA_method)
#     inds <- cov_mat_output$inds
#
#     if(standardize){
#       standardisation <- diag(cov_mat_output$mat)
#     } else {
#       standardisation <- 1
#     }
#
#     if(any(is.na(cov_mat_output$mat))){
#       warning('Not enough observations in training set to fit, loss 0 returned')
#       return(0)
#     }
#
#     mu <- colMeans(x_train[, inds, drop = F], na.rm = T)
#
#     for(j in seq_along(lambda)){
#       # do a glasso fit
#       glasso_output <- glasso::glasso(
#         cov_mat_output$mat,
#         rho = lambda[j] / sqrt(obs_share_train) * standardisation,
#         penalize.diagonal = penalize_diagonal,
#         thr = threshold
#       )
#
#       if (NA_mth %in% c('average_imputation', 'pairwise_covariance_estimation', 'loh_wainwright_bias_correction') & any(is.na(x_train))){
#         for (i in seq_len(nrep_em)) {
#           x_impute <- t(apply(x_train[, inds, drop = F], 1, impute_em, mu = mu, cov_mat = glasso_output$w))
#           glasso_output$w <- crossprod(x_impute) / nrow(x_impute)
#           if(standardize) standardisation <- diag(glasso_output$w)
#           glasso_output <- glasso::glasso(
#             glasso_output$w,
#             rho = lambda[j] / sqrt(obs_share_train) * standardisation,
#             penalize.diagonal = penalize_diagonal,
#             thr = threshold
#           )
#         }
#       }
#       loss[j] <- loglikelihood(x_test[, inds, drop = F], mu, glasso_output$w, glasso_output$wi)
#     }
#     loss
#
#   } else if (mth == 'elastic_net'){
#
#     stop('cv_inner with elastic net is not implemented yet')
#
#     family <-  args[['family']]
#     if(is.null(family)) family <- 'gaussian'
#
#     alpha <- args[['alpha']]
#     if (is.null(args[['alpha']])) alpha <- 0
#
#     glmnet_output <- glmnet(x = x_train[, -1], y = x_train[, 1],
#                             lambda = lambda / sqrt(obs_share_train), alpha = alpha, family = family)
#
#     y_pred <- predict(glmnet_output, x_test[, -1])
#
#     mean( (x_test[, 1] - y_pred)^2 )
#   }
# }





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
