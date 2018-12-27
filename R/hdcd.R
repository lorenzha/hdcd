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
                 node = NULL,
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
      node = NULL,
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
    FUN = FUN, max_depth = max_depth, node = NULL,
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


#'
#'
#'
#
cv_hdcd <- function(x, y = NULL, method = "glasso", NA_method = "complete_observations",
                    optimizer = "line_search", delta = NULL, lambda = NULL, gamma = NULL, node = NULL, alpha = NULL, control = NULL){

  mth <- match.arg(method, c("glasso", "nodewise_regression", "summed_regression", "ratio_regression", "elastic_net"))
  NA_mth <- match.arg(NA_method, c("complete_observations", "average_imputation", "loh_wainwright_bias_correction", "pairwise_covariance_estimation"))
  opt <- match.arg(optimizer, c("line_search", "section_search"))

  n <- nrow(x)

  ### Here should be something that reads parameters from control
  n_folds_outer <- control[["n_folds_outer"]]
  randomize_outer_folds <- control[["randomize_outer_folds"]]
  verbose <- control[["verbose"]]

  # This is used as long as no control functionality is enabled
  n_folds_outer <- 10
  randomize_outer_folds <- F
  verbose <- T


  if(!is.matrix(x)){
    x <- as.matrix(x)
    warning("Input data x has been coerced to matrix by hdcd.")
  }

  if(!is.null(y) & method != "elastic_net"){
    warning("Input y is ignored since method is not elastic_net")
  }
  if (!is.null(node) & method != "nodewise_regression") {
    warning("Input node is ignored since method is not nodewise_regression")
  }
  if (!is.null(alpha) & method != "elastic_net"){
    warning("Input alpha is ignored since method is not elastic_net")
  }
  if (NA_method == "complete_observations" & mean(complete.cases(x)) < 0.5){
    warning("Less than 50% of observations are complete. Consider using a different NA_method")
  }

  if (!is.null(y) & method == "elastic_net"){
    x <- cbind(y,x)
  }

  # choose lambda as grid around the asymptotic value
  if (is.null(lambda) && NCOL(x) > 1) {
    cov_mat <- get_cov_mat(x, NA_method)$mat
    lambda_max <- max(abs(cov_mat[upper.tri(cov_mat)]))
    lambda <- LogSpace(0.01 * lambda_max, lambda_max, length.out = 10)
    if (verbose) cat("Values for lambda chosen by asymptotic theory are", round(lambda, 3), "\n", sep = " ")
  }

  # choose three sensible values for delta
  if (is.null(delta)) {
    delta <- c(0.05, 0.1, 0.2)
    if (verbose) cat("Values for delta chosen are", delta, "\n", sep = " ")
  }

  folds_outer <- sample_folds(n, n_folds_outer, randomize_outer_folds)

  # do outer cross validation
  cv_results <- data.table::data.table()
  for (outer in 1:n_folds_outer){
    for (lam in lambda){
      for (del in delta){

        train_inds <- which(folds_outer != outer)
        test_inds <- which(folds_outer == outer)

        tree <- BinarySegmentation(x[train_inds, , drop = F], method = mth,
                                   NA_method = NA_mth, optimizer = opt, delta = del, lambda = lam,
                                   gamma = gamma, node = node, alpha = alpha, control = control)

        cat("tree fit finished. \n")
        gam <- gamma
        if (is.null(gamma)) {
          gam <- c(0, sort(tree$Get("max_gain")))
          gam <- gam[gam >= 0]
        }

        cat("gamma =", gamma, "\n", sep = " ")

        # get changepoints for different gammas
        res <- PruneTreeGamma(tree, gam) #this possibly creates duplcates. Fix this!
        rm(tree)

        for(g in gam){
          test_loss <- train_loss <- 0
          cpts <- train_inds[ res[["cpts"]][[as.character(g)]] ]
          alpha <- c(1, cpts , n + 1)
          for (j in 1:(length(alpha) - 1)){
            cat(g, j, sep = ', ')
            segment <- alpha[j] : (alpha[j+1] - 1)
            train_inds_segment <- intersect(train_inds, segment)
            test_inds_segment <- intersect(test_inds, segment)
            loss_output <- cv_loss(x[train_inds_segment, , drop = F], x[test_inds_segment, , drop = F], length(train_inds))
            test_loss <- test_loss + loss_output$test_loss
            train_loss <- train_loss + loss_output$train_loss
          }
          cv_results <- rbind(cv_results, data.table::data.table(fold = outer, lambda = lam, delta = del, gamma = g, test_loss = test_loss, train_loss = train_loss, cpts = list(cpts)))
        }
        if (verbose){
          temp <- cv_results[delta == del & lambda == lam & fold == outer, .(gamma, train_loss, test_loss, cpts)]
          cpts_min_train <- temp[train_loss == min(train_loss), cpts][[1]]
          cpts_min_test <- temp[test_loss == min(test_loss), cpts][[1]]
          cat("Fit finished for lambda = ", lam, " and delta = ", del,". Changepoints corresponding to minimal train / test error are ",
              cpts_min_train, " and ", cpts_min_test, " respectively. \n", sep = "")
        }
      }
    }
  }
cv_results
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
