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
                 y= NULL,
                 delta = 0.1,
                 lambda = NULL,
                 gamma = 0,
                 method = c('glasso', "nodewise_regression", "summed_regression", "ratio_regression", 'elastic_net'),
                 NA_method = 'complete_observations',
                 optimizer = c("line_search", "section_search"),
                 FUN = NULL,
                 control = list(),
                 ...) {

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

  # If a individual loss function is supplied, check that it has the required form
  if( !is.null(FUN) ){
    stopifnot('x' %in% methods::formalArgs(FUN))
    if ( is.null(lambda) && !('lambda' %in% methods::formalArgs(FUN(x))))
      lambda <- 0 #don't do lambda CV if FUN doesn't depend on lambda
  }

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
