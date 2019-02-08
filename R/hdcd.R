#' hdcd
#'
#' High Dimensional Changepoint Detection
#'
#' @inheritParams BinarySegmentation
#' @inheritParams CrossValidation
#'
#' @return For a single fit a list with elements
#'
#' \describe{
#'   \item{changepoints}{A numeric list with the indices of the changepoints}
#'   \item{tree}{The fully grown binary tree}
#' }
#'
#' For cross-validation a list with elements
#'
#' #' \describe{
#'   \item{changepoints}{A numeric list with the indices of the changepoints}
#'   \item{cv_results}{A multi-dimensional array with the cross-validation results}
#'   \item{cv_gamma}{Best gamma value}
#'   \item{cv_lambda}{Best lambda value}
#'   \item{cv_delta}{Best delta value}
#' }
#'
#' If only a single fit was performed a list with the found changepoints as well as the fully grown binary tree are returned. For cross-validation the a list with the found changepopints, the optimal parameter values and the full results is returned.
#' @export
#'
#' @examples
#' dat <- SimulateFromModel(CreateModel(n_segments = 2,n = 100,p = 30, ChainNetwork))
#' \dontrun{
#' hdcd(dat, 0.1, 0.1, 0.05, method = "summed_regression", verbose = T)
#' }
hdcd <- function(x,
                 lambda,
                 y= NULL,
                 delta = 0.1,
                 gamma = 0,
                 method = c('glasso', "nodewise_regression", "summed_regression", "ratio_regression", 'elastic_net'),
                 NA_method = 'complete_observations',
                 optimizer = c("line_search", "section_search"),
                 FUN = NULL,
                 control = list(),
                 ...) {

  verbose <- control_get(control, "verbose", F)

  if(!is.matrix(x)){
    x <- as.matrix(x)
    warning("Input data x has been coerced to matrix by hdcd.")
  }
  stopifnot(nrow(x) > 1)

  # If y is supplied, it is bound to x
  if(!is.null(y)){
    stopifnot(nrow(x) == length(y))
    x <- cbind(y, x)
  }

  cv_inner <- control_get(control, "cv_inner", FALSE)
  if(is.null(control[["lambda_inner"]]) & cv_inner){
    lambda_inner_min_ratio = control_get(control, "lambda_inner_min_ratio", 0.01)
    lambda_inner_grid_size = control_get(control, "lambda_inner_grid_size", 4)
    # choose lambda as grid around the asymptotic value
    cov_mat <- get_cov_mat(x, NA_method)$mat
    lambda_max <- max(abs(cov_mat[upper.tri(cov_mat)]))
    control$lambda_inner <- LogSpace(lambda_inner_min_ratio * lambda_max, lambda_max, length.out = lambda_inner_grid_size)
    if (verbose) cat('lambda inner set by asymptotic theory \n')
  }

  # If a individual loss function is supplied, check that it has the required form
  if( !is.null(FUN) ){
    stopifnot('x' %in% methods::formalArgs(FUN))
    if ( is.null(lambda) && !('lambda' %in% methods::formalArgs(FUN(x))))
      lambda <- 0 #don't do lambda CV if FUN doesn't depend on lambda
  }

  # NA_mth <- match.arg(NA_method)
  # if(NA_mth == 'complete_observations'){
  #   cases <- complete.cases(x)
  #   x <- x[cases, ]
  #   if (any(!cases)){
  #     warning(paste('There are', sum(!cases), 'incomplete cases in x that are discarded. Maybe choose another NA_method', sep = ' '))
  #   }
  #   train_inds <- which(cases)
  # } else {
  #   train_inds <- 1 : nrow(x)
  # }

  tree <- BinarySegmentation(
    x = x, delta = delta, lambda = lambda, method = method, NA_method = NA_method,
    optimizer = optimizer, control = control, FUN = FUN,
    ...
  )

  tree
}





#' Print method for objects of class bs_cv
#'
#' @param x Object of class bs_cv
#' @param ... Only included to be consist with print generic.
#'
#' @export
print.bs_cv <- function(x, ...) {
  stopifnot("bs_cv" %in% class(x))

  cat("\nBest parameters:\n\n")
  cat("-> Lambda: ", x$cv_lambda, "\n")
  cat("-> Delta : ", x$cv_delta, "\n")
  cat("-> Gamma : ", x$cv_gamma, "\n")
  cat("\nFound changepoints: \n\n", x$changepoints, "\n")
}
