#' hdcd
#'
#' High Dimensional Changepoint Detection
#'
#' @param x A n times p matrix for which to find the best splitting point.
#' @param delta Value between 0 and 0.5. Tuning param which determines the minimal segment size
#' proportional to the size of the dataset and hence an upper bound for the number of changepoints. A good value to start with is 0.1.
#' @param lambda Sparsity penality parameter in single lasso fits. If NULL k-fold cross-validation will be conducted.
#' @param gamma Split penalty parameter for pruning the tree. If NULL k-fold cross-validation will be conducted.
#' @param method Which method should be used for fitting the model? See defaults for possible choices.
#' @param penalize_diagonal Boolean, should the diagonal elements of the precision matrix be penalized?
#' @param threshold The threshold for halting the iteration in glasso or glmnet. In the former it controls the change of single parameters
#' in the latter it controls the total objective value.
#' @param use_ternary_search Use a ternary search algorithm in each level of the recursion to find a local optimum (EXPERIMENTAL)
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
                 use_ternary_search = F,
                 standardize = T,
                 threshold = 1e-7,
                 n_folds = 10,
                 verbose = F,
                 parallel = T,
                 ...) {

  stopifnot(nrow(x) > 1)
  x_mat <- as.matrix(x)
  mth   <- match.arg(method)

  if(is.null(lambda) || is.null(gamma) || is.null(delta)){
    cv <- TRUE
    if (verbose) cat("\nPerforming ",n_folds,"- fold cross-validation...\n")
    cv_res <- CrossValidation(x = x_mat, delta = delta, method = mth, lambda = lambda,
                              gamma = gamma, n_folds = n_folds,
                              use_ternary_search = use_ternary_search,
                              standardize = standardize,
                              penalize_diagonal = penalize_diagonal,
                              verbose = verbose,
                              parallel = parallel,
                              threshold = threshold,
                              ...)
    lambda <- cv_res$best_lambda
    gamma  <- cv_res$best_gamma
    delta  <- cv_res$best_delta
  }

  tree <- BinarySegmentation(x = x_mat, delta = delta, lambda = lambda, method = mth,
                             threshold = threshold, penalize_diagonal = penalize_diagonal,
                             use_ternary_search = use_ternary_search, standardize = standardize, ...)
  res <- PruneTreeGamma(tree, gamma)
  if (verbose){
    cat("\nFinal tree for cross-validated gamma and lambda:\n \n")
    print(res[["pruned_tree"]])
  }

  if (cv){
    res <- list(changepoints = res[["cpts"]][[1]], cv_results = cv_res[["cv_results"]],
                cv_gamma = gamma, cv_lambda = lambda)
    class(res) <- "bs_cv"
  } else {
    res <- list(changepoints = res[["cpts"]][[1]])
    class(res) <- "bs"
  }
  res
}
