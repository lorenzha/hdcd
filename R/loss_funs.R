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

    x_est <- x[start : end, ]

    mis_completely <- apply(x_est, 2, function(y) (sum(!is.na(y)) < 2))

    for(i in 1:n_p){
      if (mis_completely[i]){
        x_est[is.na(x_est[, i]), i] <- rnorm(sum(is.na(x_est[, i])))
      } else {
        x_est[is.na(x_est[, i]), i] <- mean(x_est[, i], na.rm=T)
      }
    }

#     if (any(!is.finite(x_est))){
#       x_est[!is.finite(x_est)] <- runif(sum(!is.finite(x_est)))
#       ('warning in av_Segloss: whole pred missing NA')
#     }

    obs_count <- end - start + 1
    obs_share <- obs_count / n_obs

    # We need more than one observation to caclculate the covariance matrix
    stopifnot(obs_count > 1)


    cov_mat <- (obs_count - 1) / obs_count * cov(x_est)
    #cov_mat <- as.matrix(Matrix::nearPD(cov_mat)$mat)


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
    cov_mat[is.na(cov_mat)] <- runif(sum(is.na(cov_mat))) # adapt value
    cov_mat <- as.matrix(Matrix::nearPD(cov_mat, posd.tol = 1e-3, eig.tol = 1e-3)$mat)

    mu <- apply(x_est, 2, function(y) mean(y, na.rm = T))
    mu[!is.finite(mu)] <- 0
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

      if (all(mis_i)) {
        x_est[i, ] <- mu
      } else if(sum(!mis_i) == 1) {
        x_est_i <- x_est[i, ]
        x_est_i[mis_i] <- mu[mis_i] - solve(K[mis_i, mis_i], K[mis_i, !mis_i] * (x_est_i[!mis_i] - mu[!mis_i]))
        x_est[i, ] <- x_est_i
      } else if (any(mis_i)) {
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


#' cv.Loss
#
#

cv.Loss <- function(x_train, x_test,
                    lambda = 0,
                    penalize_diagonal = FALSE,
                    standardize = TRUE,
                    threshold = 1e-07,
                    method = c("nodewise_regression", "summed_regression", "ratio_regression", 'glasso', 'elastic_net'),
                    NA_method = c('complete_observations', 'pairwise_covariance_estimation', 'loh_wainwright_bias_correction'),
                    ...){

  mth <- match.arg(method)
  args <- list(...)

  stopifnot(ncol(x_train) == ncol(x_test))
  n_obs_train <- nrow(x_train)
  n_obs_test <- nrow(x_test)

  SegmentLossFUN <- SegmentLoss(x_train, lambda = lambda, penalize_diagonal = penalize_diagonal,
                                threshold = threshold, method = method, NA_method = NA_method, cv = TRUE, ...)

  if (mth %in% c("nodewise_regression", "summed_regression", 'glasso', "ratio_regression")){
    function(start_train, end_train, start_test, end_test){

      stopifnot(end_train - start_train >= 1 & end_test - start_test >= 0)

      glasso_output <- SegmentLossFUN(start_train, end_train)

      #R <- chol(glasso_output$w)

      mu <- colMeans(x[start_train : end_train, glasso_output$inds, drop = F], na.rm = T)

      n_eff_obs_test <- sum(!is.na(x_test[start_test : end_test, ])) / ncol(x_test)

      # optimize this!
      lossfun <- function(y){
        if(all(is.na(y))){
          0
        } else {
        Sigma <- glasso_output$w[!is.na(y), !is.na(y), drop = F]
       # R_cur <- R
        v <- y - mu
        v <- v[!is.na(y)]
       # R_cur[, is.na(y)] <- 0
        #diag(R_cur)[is.na(y)] <- 1
        t(v) %*% solve(Sigma, v) + log(abs(det(Sigma)))
        }
      }

      sum(apply(x_test[start_test : end_test, glasso_output$inds, drop = F], 1, lossfun))
      # (- log(det(abs(glasso_output$wi))) + sum(diag( cov_mat_test %*% glasso_output$wi )))
    }

  } else if (mth == 'elastic_net'){

    family <-  args[['family']]
    if(is.null(family)) family <- 'gaussian'

    alpha <- args[['alpha']]
    if (is.null(args[['alpha']])) alpha <- 0

    function(start_train, end_train, start_test, end_test){

      obs_share_train <- (end_train - start_train + 1) / n_obs_train

      glmnet_output <- glmnet(x = x_train[start_train : end_train, -1], y = x_train[start_train : end_train, 1],
                              lambda = lambda / sqrt(obs_share_train), alpha = alpha, family = family)
      y_pred <- predict(glmnet_output, x_test[start_test : end_test, -1])

      sum( (x_test[start_test : end_test, 1] - y_pred)^2 ) / (end_test - start_test + 1)
    }
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
                        lambda = 0,
                        penalize_diagonal = FALSE,
                        standardize = TRUE,
                        threshold = 1e-07,
                        use_cache = TRUE,
                        cv = FALSE,
                        NA_method = 'complete_observations',
                        method = c("nodewise_regression", "summed_regression", "ratio_regression", "glasso", "elastic_net"),
                        ...) {


  args <- list(...)
  mth <- match.arg(method)
  n_obs <- NROW(x)
  n_p <- NCOL(x)

  if (mth == "glasso") {
    get_cov_mat <- get_covFUN(x, NA_method)

    get_loss <- function(start, end, ...) {

      obs_count <- end - start + 1

      stopifnot(obs_count > 1) # We need more than one observation to calculate the covariance matrix

      cov_mat_output <- get_cov_mat(start, end)
      inds <- cov_mat_output$inds # which predictors have sufficient non-missing values
      cur_p <- sum(inds)
      cov_mat <- cov_mat_output$mat

      obs_share <- cov_mat_output$n_eff_obs / n_obs

      if (cur_p == 0){
        return(NA)
      }

      glasso_output <- glasso::glasso( # correction with log(cur_p) / log(n_p) / sqrt(obs_share) from asymptotic theory
        cov_mat,
        rho = lambda / sqrt(obs_share) * log(cur_p) / log(n_p) * diag(cov_mat),
        penalize.diagonal = penalize_diagonal,
        thr = threshold
      )

      if (cv){
        glasso_output[['inds']] <- inds
        glasso_output
      } else {
        if (!penalize_diagonal) {
          diag(glasso_output$wi) <- 0
        }
        # Needed to undo transformation of likelihood in glasso package
        # equivalent to ' - log(det(glasso_output$wi)) + psych::tr(cov_mat %*% glasso_output$wi)' before setting diag() <- 0
        #- 2 * glasso_output$loglik / n_p
        #      - lambda  / sqrt(obs_share) * as.numeric(sqrt(diag(cov_mat)) %*% abs(glasso_output$wi) %*% sqrt(diag(cov_mat) ))
        # MESS!!! TODO: cleanup!
        n_p / cur_p * (((glasso_output$loglik / (-cur_p / 2) # Needed to undo transformation of likelihood in glasso package
          - lambda * (log(cur_p) / log(n_p)) / sqrt(obs_share) * sum(abs(glasso_output$wi))) * obs_share)) # Remove regularizer added in glasso package

        #n_p * 2 * -glasso_output$loglik - (n_p / cur_p) * lambda * log(cur_p) / log(n_p) * sum(abs(glasso_output$wi)) * sqrt(obs_share)
      }
    }

  } else if (mth == "nodewise_regression") {

    p <- args[["node"]]
    stopifnot(length(p) == 1 && is.numeric(p))

    get_loss <- function(start, end, ...) {
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
    get_loss <- function(start, end, ...) {
      obs_count <- end - start + 1
      obs_share <- obs_count / n_obs

      stopifnot(obs_count > 1)  # We need more than one observation to calculate the covariance matrix

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
    get_loss <- function(start, end, ...) {
      obs_count <- end - start + 1
      obs_share <- obs_count / n_obs

      stopifnot(obs_count > 1) # We need more than one observation to calculate the covariance matrix

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

    alpha <- args[['alpha']]
    stopifnot(is.numeric(alpha) && 0 <= alpha && alpha <= 1)

    family <-  args[['family']]
    if (is.null(family)){
      family <- 'gaussian'
    }

    get_loss <- function(start, end, ...){

      obs_count <- end - start + 1
      obs_share <- obs_count / n_obs
      fit <- glmnet::glmnet(x[start : end, -1], x[start : end, 1], alpha = alpha, lambda = sqrt(obs_share) * lambda, standardize = standardize, family = family, thres = threshold)
      deviance(fit) / n_obs
    }
  }


  # If
  if (!use_cache){
    get_loss
  } else {
    cache <- new.env(TRUE, emptyenv())
    key <- function(start, end){
        paste(start, end, sep = '-')
    }
    function(start, end, ...){
      if(!exists(key(start, end), envir = cache, inherits = FALSE)){
        assign(key(start, end),
               get_loss(start, end, ...),
               envir = cache)
      }
    get(key(start, end), envir = cache, inherits = FALSE)
    }
  }
}


get_covFUN <- function(x, NA_method = c('complete_observations', 'pairwise_covariance_estimation',
                                        'loh_wainwright_bias_correction')){
  NA_mth <- match.arg(NA_method)
  if (NA_mth == 'complete_observations'){
    function(start, end){
      cov_mat <- cov(x[start : end, , drop = F], use = 'na.or.complete')
      list(mat = cov_mat, inds = rep(T, ncol(x)), n_eff_obs = sum(complete.cases(x[start : end, , drop = F])))
    }

  } else if (NA_mth  == 'pairwise_covariance_estimation'){
    av_obs <- rbind(rep(0, ncol(x)), apply(x, 2, function(y) cumsum(!is.na(y))))
    function(start, end){
      cur_av <- av_obs[end + 1, ] - av_obs[start, ]
      obs_count <- end - start + 1
      if(all(cur_av == obs_count)){
        list(mat = cov(x[start : end, , drop = F]), inds = rep(T, ncol(x)), n_eff_obs = obs_count)
      } else {
        inds <- (cur_av >= 2) # for which predictors can we even estimate covariance?
        cov_mat <- cov(x[start : end, inds, drop = F], use = 'pairwise')
        cov_mat[is.na(cov_mat)] <- 0
        cov_mat <- (obs_count - 1) / obs_count * as.matrix(Matrix::nearPD(cov_mat)$mat) #TODO: this de-debiasing should be differen
        list(mat = cov_mat, inds = inds, n_eff_obs = sum(cur_av) / sum(inds))
      }
    }

  } else if (NA_mth == 'loh_wainwright_bias_correction') {
    av_obs <- rbind(rep(0, ncol(x)), apply(x, 2, function(y) cumsum(!is.na(y))))
    function(start, end){
      cur_av <- av_obs[end + 1, ] - av_obs[start, ]
      obs_count <- end - start + 1
      if(all(cur_av == obs_count)){
        list(mat = cov(x[start : end, , drop = F]), inds = rep(T, ncol(x)), n_eff_obs = obs_count)
      } else {
        obs_count <-  end - start + 1
        inds <- (cur_av >= 2) # need at least two observations to estimate variance
        stopifnot(any(inds))
        mis_frac = 1 - cur_av[inds] / obs_count  # estimate prob of misssingness per column
        z <- scale(x[start : end, inds, drop = F], T, F) # center
        z[is.na(z)] <- 0 # impute with mean (since centered)
        z <- z %*% diag(1 - mis_frac) #first step in the loh-wainwright correction
        cov_mat <- cov(z)
        cov_mat <- cov_mat - diag(mis_frac * diag(cov_mat)) # second step loh-wainwright correction
        n_eff_obs <- sum(cur_av) / sum(inds)
        cov_mat <- (obs_count - 1) / n_eff_obs * as.matrix(Matrix::nearPD(cov_mat)$mat) #TODO: divide bz n_eff_obs???
        list(mat = cov_mat, inds = inds, n_eff_obs = n_eff_obs)
      }
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
