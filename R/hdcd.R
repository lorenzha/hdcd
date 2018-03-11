#' hdcd
#'
#' High Dimensional Changepoint Detection
#'
#' @inheritParams BinarySegmentation
#' @param gamma Split penalty parameter for pruning the tree. If NULL k-fold cross-validation will be conducted.
#' @param n_folds Number of folds in cross validation. Default is 10.
#' @param verbose Should additional information be printed?
#' @param ... Supply additional arguments for a specific method (e.g. p for nodewise_regression)
#'
#' @return The indices of the found changepoints in the sequence
#' @export
#'
#' @examples
#' dat <- SimulateFromModel(CreateModel(n_segments = 2,n = 100,p = 30, ChainNetwork))
#' hdcd(dat, 0.1, 0.1, 0.05, method = "summed_regression", verbose = T)
hdcd <- function(x,
                 delta = 0.1,
                 lambda = NULL,
                 gamma = NULL,
                 method = c("nodewise_regression", "summed_regression", "ratio_regression"),
                 penalize_diagonal = F,
                 optimizer = c("line_search", "ternary_search", "section_search"),
                 control = NULL,
                 standardize = T,
                 threshold = 1e-7,
                 n_folds = 10,
                 verbose = F,
                 parallel = T,
                 ...) {
  stopifnot(nrow(x) > 1)
  x_mat <- as.matrix(x)
  cv <- FALSE

  if ((is.null(lambda) || is.null(gamma) || is.null(delta)) | length(c(gamma, lambda, delta)) > 3) {
    cv <- TRUE
    if (verbose) cat("\nPerforming ", n_folds, "- fold cross-validation...\n")
    cv_res <- CrossValidation(
      x = x_mat, delta = delta, method = method, lambda = lambda,
      gamma = gamma, n_folds = n_folds,
      optimizer = optimizer,
      control = control,
      standardize = standardize,
      penalize_diagonal = penalize_diagonal,
      verbose = verbose,
      parallel = parallel,
      threshold = threshold,
      ...
    )
    lambda <- cv_res$best_lambda
    gamma <- cv_res$best_gamma
    delta <- cv_res$best_delta
  }

  tree <- BinarySegmentation(
    x = x_mat, delta = delta, lambda = lambda, method = method,
    threshold = threshold, penalize_diagonal = penalize_diagonal,
    optimizer = optimizer, control = control, standardize = standardize,
    ...
  )
  res <- PruneTreeGamma(tree, gamma)
  if (verbose) {
    cat("\nFinal tree for cross-validated gamma and lambda:\n \n")
    print(res[["pruned_tree"]])
  }

  if (cv) {
    res <- list(
      changepoints = res[["cpts"]][[1]], cv_results = cv_res[["cv_results"]],
      cv_gamma = gamma, cv_lambda = lambda, cv_delta = delta
    )
    class(res) <- "bs_cv"
  } else {
    res <- list(changepoints = res[["cpts"]][[1]])
    class(res) <- "bs"
  }
  res
}
