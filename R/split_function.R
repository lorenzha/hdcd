gain <- function(x, lambda = NULL,
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
  lambda_global <- lambda

  if (mth != "glasso") {
    stop('not supported for anything except glasso yet')
  } else {

    function(start, end, lambda = NULL) {

      if(!is.null(lambda_global)){
        lambda <- lambda_global
      }

      glasso_output_full <- make_glasso_fit(x[start  : end , , drop = F], lambda, n_obs, NA_method, control)

      function(split_point){

        if(is.null(split_point)){
          loglikelihood(x[start : end, glasso_output_full$inds, drop = F],
                        glasso_output_full$mu[glasso_output_full$inds, drop = F],
                        glasso_output_full$w[glasso_output_full$inds, glasso_output_full$inds, drop = F],
                        glasso_output_full$wi[glasso_output_full$inds, glasso_output_full$inds, drop = F])
        } else {

          stopifnot(split_point >= start + 2 & end >= split_point + 1)

          glasso_output_left <- make_glasso_fit(x[start  : (split_point - 1) , , drop = F], lambda, n_obs, NA_method, control)
          glasso_output_right <- make_glasso_fit(x[split_point  : end, , drop = F], lambda, n_obs, NA_method, control)

          (loglikelihood(x[start : (split_point - 1), glasso_output_left$inds, drop = F],
                        glasso_output_full$mu[glasso_output_left$inds, drop = F],
                        glasso_output_full$w[glasso_output_left$inds, glasso_output_left$inds, drop = F],
                        glasso_output_full$wi[glasso_output_left$inds, glasso_output_left$inds, drop = F]) +
            loglikelihood(x[split_point : end, glasso_output_right$inds, drop = F],
                          glasso_output_full$mu[glasso_output_right$inds, drop = F],
                          glasso_output_full$w[glasso_output_right$inds, glasso_output_right$inds, drop = F],
                          glasso_output_full$wi[glasso_output_right$inds, glasso_output_right$inds, drop = F]) -
            loglikelihood(x[start : (split_point - 1), glasso_output_left$inds, drop = F],
                          glasso_output_left$mu[glasso_output_left$inds, drop = F],
                          glasso_output_left$w[glasso_output_left$inds, glasso_output_left$inds, drop = F],
                          glasso_output_left$wi[glasso_output_left$inds, glasso_output_left$inds, drop = F]) -
            loglikelihood(x[split_point : end, glasso_output_right$inds, drop = F],
                          glasso_output_right$mu[glasso_output_right$inds, drop = F],
                          glasso_output_right$w[glasso_output_right$inds, glasso_output_right$inds, drop = F],
                          glasso_output_right$wi[glasso_output_right$inds, glasso_output_right$inds, drop = F])
            )/n_obs
        }
      }
    }
  }
}


make_glasso_fit <- function(x, lambda, n_obs, NA_method, control){

  penalize_diagonal <-  control$glasso_penalize_diagonal
  standardize <-  control$glasso_standardize
  threshold <- control$glasso_threshold
  min_frac <- control$segment_loss_min_frac

  p <- ncol(x)
  cov_mat_output <- get_cov_mat(x, NA_method, min_frac = min_frac)

  stopifnot(any(cov_mat_output$inds))

  obs_share <- cov_mat_output$n_eff_obs / n_obs # for rescaling of lambda

  out <- list(w = matrix(NA, nrow = p, ncol = p), wi = matrix(NA, nrow = p, ncol = p))

  if (standardize){
    glasso_output <- glasso::glasso(
      cov_mat_output$mat,
      rho = lambda / sqrt(obs_share) * diag(cov_mat_output$mat),
      penalize.diagonal = penalize_diagonal,
      thr = threshold
    )
  } else {
    glasso_output <- glasso::glasso(
      cov_mat_output$mat,
      rho = lambda / sqrt(obs_share),
      penalize.diagonal = penalize_diagonal,
      thr = threshold
    )
  }

  out$w[cov_mat_output$inds, cov_mat_output$inds] <- glasso_output$w
  out$wi[cov_mat_output$inds, cov_mat_output$inds] <- glasso_output$wi
  out$inds <- cov_mat_output$inds
  out$mu <- colMeans(x, na.rm = T)

  out
}


# calculates loglikelihood of x (with missing values) given mu, cov_mat. cov_mat_inv can be used to speed up calculaitons
# of the inverse of the submatrix via the Woodbury matrix identity as a low rank update from cov_mat_inv
# Maybe improve speed??
# loglikelihood <- function(x, mu, cov_mat, cov_mat_inv, standardize_loglik = F){
#   if(!is.matrix(x)){
#     x <- as.matrix(x)
#     warning('non matrix input into loglikelihood')
#   }
#   p <- ncol(x)
#   inds <- is.na(x) # missingness structure of x
#   na_order <- do.call(order, lapply(1:ncol(inds), function(i) inds[, i])) # order x by missingness structure
#   cov_mat_inv_cur <- cov_mat_inv
#   log_det_cur <- determinant(cov_mat, logarithm = TRUE)$modulus
#   attributes(log_det_cur) <- NULL
#   loss <- 0
#
#   inds_old <- rep(F, ncol(x))
#   k <- 0
#   for(i in 1:nrow(x)){
#     inds_cur <- inds[na_order[i], ]
#     if(any(inds_cur != inds_old)){
#       inds_old <- inds_cur
#       if (sum(inds_cur) < length(inds_cur) / 3){
#         # fast calculation of cov_mat_inv_cur as update from cov_mat_inv via Woodbury matrix identity
#         k <- sum(inds[na_order[i], ]) # amount of missing values
#         A <- diag(p)[, inds_cur, drop = F] # helper matrix
#         U <- cbind(A, cov_mat[, inds_cur])
#         U[, (k + 1) : (2*k)][as.logical(A)] <- 0
#         V <- rbind(cov_mat[inds_cur, ], t(A))
#         V[1 : k, ][as.logical(t(A))] <- 0
#         V[1 : k, ][as.logical(t(A[, k : 1]))] <- 0
#         S <- -diag(2 * k) + V %*% cov_mat_inv %*% U
#         cov_mat_inv_cur <- tryCatch({
#           #solve(cov_mat[!inds_cur, !inds_cur])
#           (cov_mat_inv - cov_mat_inv %*% U %*% solve(S) %*% V %*% cov_mat_inv)[!inds_cur, !inds_cur]
#         }, condition = function(e){
#           #warning('loglikelihood was not able to apply the Woodbury matrix identity')
#           solve(cov_mat[!inds_cur, !inds_cur])
#         })
#       } else {
#         cov_mat_inv_cur <- solve(cov_mat[!inds_cur, !inds_cur])
#       }
#       log_det_cur <- determinant(cov_mat[!inds_cur, !inds_cur], logarithm = TRUE)$modulus
#       attributes(log_det_cur) <- NULL
#     }
#     v <- (x[na_order[i], ] - mu)[!inds_cur]
#     if (standardize_loglik){
#       loss <- loss + (t(v) %*% cov_mat_inv_cur %*% v - p) / sqrt(p / (p - k))
#     } else {
#       loss <- loss + t(v) %*% cov_mat_inv_cur %*% v + log_det_cur #+ (p - k) * log(2*pi)  #This p+k is essential. WHY??
#       #cat(i, ', addition to loss = ',  t(v) %*% cov_mat_inv_cur %*% v, ' + ', log_det_cur, '\n')
#     }
#   }
#   loss / 2
# }


loglikelihood <- function(x, mu, cov_mat, cov_mat_inv, standardize_loglik = F){
  if(!is.matrix(x)){
    x <- as.matrix(x)
    warning('non matrix input into loglikelihood')
  }
  n <- nrow(x)
  p <- ncol(x)
  loss <- numeric(n)

  for (i in 1 : nrow(x)){
    inds <- !is.na(x[i, ])
    if(any(!inds)){
      log_det <- determinant(cov_mat[inds, inds], logarithm = TRUE)$modulus
      v <- x[i, inds] - mu[inds]
      distance <- t(v) %*% solve(cov_mat[inds, inds], v) ##vectorize, b = matrix, for rows of same missingness structure
      loss[i] <- distance + log_det + sum(inds) * log(2*pi)
    } else {
      v <- x - mu
      distance <- t(v) %*% cov_mat_inv %*% v
      loss[i] <- distance + log_det_full + p * log(2*pi)
    }
  }


  sum(loss) / 2
}
