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
                 delta = NULL,
                 lambda = NULL,
                 lambda_min_ratio = 0.01,
                 lambda_grid_size = 10,
                 gamma = NULL,
                 method = c('glasso', "nodewise_regression", "summed_regression", "ratio_regression", 'elastic_net'),
                 NA_method = 'complete_observations',
                 penalize_diagonal = F,
                 alpha = 1,
                 optimizer = c("line_search", "section_search"),
                 control = NULL,
                 standardize = T,
                 threshold = 1e-7,
                 n_folds = 10,
                 verbose = T,
                 parallel = T,
                 FUN = NULL,
                 max_depth = Inf,
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

  cv <- FALSE

  # If a individual loss function is supplied, check that it has the required form
  if( !is.null(FUN) ){
    stopifnot('x' %in% methods::formalArgs(FUN))
    if ( is.null(lambda) && !('lambda' %in% methods::formalArgs(FUN(x))))
      lambda <- 0 #don't do lambda CV if FUN doesn't depend on lambda
  }

  # do cross-validation if any of lambda, gamma or delta is not supplied
  if ((is.null(lambda) || is.null(gamma) || is.null(delta)) | length(c(gamma, lambda, delta)) > 3) {
    cv <- TRUE
    if (verbose) cat("\nPerforming ", n_folds, "- fold cross-validation...\n")

    cv_res <- CrossValidation(
      x = x, delta = delta, method = method, lambda = lambda,
      lambda_min_ratio = lambda_min_ratio, lambda_grid_size = lambda_grid_size,
      NA_method = NA_method,
      gamma = gamma, n_folds = n_folds,
      optimizer = optimizer,
      control = control,
      standardize = standardize,
      penalize_diagonal = penalize_diagonal,
      alpha = alpha,
      verbose = verbose,
      parallel = parallel,
      threshold = threshold,
      FUN = FUN,
      max_depth = max_depth,
      ...
    )
    lambda <- cv_res$opt$lambda
    gamma <- cv_res$opt$gamma
    delta <- cv_res$opt$delta
  }

  tree <- BinarySegmentation(
    x = x, delta = delta, lambda = lambda, method = method, NA_method = NA_method,
    threshold = threshold, penalize_diagonal = penalize_diagonal, alpha = alpha,
    optimizer = optimizer, control = control, standardize = standardize,
    FUN = FUN, max_depth = max_depth,
    ...
  )

  # Prune the tree to the cross validated gamma
  res <- PruneTreeGamma(tree, gamma)

  if (cv) {
    res <- list(
      res = res[['pruned_tree']][[1]], changepoints = res[["cpts"]][[1]], cv_results = cv_res[["cv_results"]],
      cv_gamma = gamma, cv_lambda = lambda, cv_delta = delta
    )
    if (verbose) cat('\nFinal tree for cross-validated gamma = ', gamma,' and lambda = ', lambda,':\n \n')
    class(res) <- "bs_cv"
  } else {
    res <- list(changepoints = res[["cpts"]][[1]], tree = tree)
    class(res) <- "bs"
  }
  print(res)
  res
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
