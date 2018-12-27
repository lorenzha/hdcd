#' HandleGlassoNaN
#'
#' Suppress specific warning in glasso package
#'
#' To be used with withCallingHandlers. Glasso produces non-pd covariance matrices
#' when using the Meinshausen-Bühlmann approximation. It then tries to calculate the
#' likelihood an hence log(det(cov_mat_i)) which will throw a warning if cov_mat is
#' not pd.
#'
#' @param w Warning to be handled
HandleGlassoNaN <- function(w) {
  if (any(grepl("NaNs produced", w))) {
    invokeRestart("muffleWarning")
  }
}

catch <- function(x){
  if (length(x) > 0){
    x
  } else {
    NA
  }
}

sample_folds <- function(n, k, randomize = TRUE){
  if (k == 1){
    as.factor(rep(0, n))
  } else  if (randomize){
    random_draw <- runif(n)
    k_quantiles <- quantile(random_draw, 0:k/k)
    cut(random_draw, k_quantiles, labels = 1:k, include.lowest = TRUE)
  } else {
    as.factor(rep(1:k, ceiling(n/k))[1:n])
  }
}

hdcd_control <- function(){
  control <- list()
}
