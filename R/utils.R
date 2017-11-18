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
HandleGlassoNaN <- function(w){
  if(any(grepl("NaNs produced", w)))
    invokeRestart("muffleWarning")
}
