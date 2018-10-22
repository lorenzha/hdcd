NA_SegmentLoss <- function(NA_method){

  NA_mth <- match.arg(NA_method, c('interpolation', 'averaging', 'pairwise_covariance'))

  if (NA_mth == 'interpolation'){

    function(lambda,
             n_obs,
             penalize_diagonal = FALSE,
             standardize = TRUE,
             threshold = 1e-07){

      n_p <- ncol(x)

      function(x, start, end){

        x_interp <- na.interpolation(x)

        obs_count <- NROW(x)
        obs_share <- obs_count / n_obs

        # We need more than one observation to calculate the covariance matrix
        stopifnot(obs_count > 1)

        n_p <- NCOL(x_interp)

        cov_mat <- (obs_count - 1) / obs_count * cov(x_interp)

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
  } else if (NA_mth == 'averaging'){

    function(n_obs,
             lambda,
             penalize_diagonal = FALSE,
             standardize = TRUE,
             threshold = 1e-07){

      function(x, start, end){

        n_p <- NCOL(x)
        x_av <- x

        for(i in 1:n_p){
          x_av[,i][is.na(x[,i]==T)] <- mean(x[,i], na.rm=T)
        }

        obs_count <- NROW(x)
        obs_share <- obs_count / n_obs

        # We need more than one observation to caclculate the covariance matrix
        stopifnot(obs_count > 1)



        cov_mat <- (obs_count - 1) / obs_count * cov(x_av)

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
  } else if (NA_mth == 'pairwise_covariance'){

    function(n_obs,
             lambda,
             penalize_diagonal = FALSE,
             standardize = TRUE,
             threshold = 1e-07){

      function(x, start, end){

        obs_count <- NROW(x)
        obs_share <- obs_count / n_obs

        # We need more than one observation to caclculate the covariance matrix
        stopifnot(obs_count > 1)

        n_p <- NCOL(x)

        cov_mat <- (obs_count - 1) / obs_count * cov(x, use='pairwise')

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
  } else {
    print('error, loss funciton trivial')
  }
}
