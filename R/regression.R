#' FullRegression
#'
#' Neighbourhood selection for a given set of changepoints and a data sequence
#'
#' @param x A n times p data matrix
#' @param cpts A set of changepoints correspoding to row indices.
#' @param lambda Sparsity penalty parameter
#'
#' @return \describe{
#'   \item{est_coeffs}{A p times p matrix of estimated coefficients.}
#'   \item{est_mean}{A p vector of estimated means.}
#'   \item{est_intercepts}{A p vector of estimated intercepts.}
#' }
FullRegression <- function(x, cpts, lambda = 0.1) {
  est_mean <- est_wi <- est_intercepts <- list()
  n_obs <- nrow(x)

  cpts <- c(1, cpts, (nrow(x) + 1))

  for (i in 1:(length(cpts) - 1)) {
    start <- cpts[i]
    end <- cpts[i + 1] - 1

    obs_count <- end - start + 1
    obs_share <- obs_count / n_obs

    cov_mat <- (obs_count - 1) / obs_count * cov(x[start:end, ])
    est_mean[[i]] <- colMeans(x[start:end, ])
    withCallingHandlers({
    est_coeffs[[i]] <- glasso::glasso(cov_mat, rho = lambda / sqrt(obs_share) * diag(cov_mat), approx = T)$wi
    }, warning = HandleGlassoNaN)
    est_intercepts[[i]] <- est_mean[[i]] - colSums(est_coeffs[[i]] * est_mean[[i]])
  }
  list(est_coeffs = est_coeffs, est_mean = est_mean, est_intercepts = est_intercepts)
}