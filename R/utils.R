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

#returns, when specified, control$name, else returns default
control_get <- function(control, name, default){
  value <- control[[name]]
  if(!is.null(value)){
    value
  } else {
    default
  }
}

# turns numeric(0) into NA
catch <- function(x){
  if (length(x) > 0){
    x
  } else {
    NA
  }
}

# returns folds sampled randomly / equispaced folds. If k == 1 returns list of zeros
# such that x_train = x[inds != 1, ] = x.
sample_folds <- function(n, k, randomize = FALSE){
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

# adapt this
hdcd_control <- function(){
  control <- list()
}
