#' CrossValidation
#'
#' EXPERIMENTAL. Crossvalidation for the desired method and parameter combinations.
#' Evaluating different lambda values will lead to refitting the entire model whereas gamma values can be evaluted cheaply using a
#' single fit.
#'
#' Will make use of a registered parallel backend over the folds and lambda values. Therefore the whole cross-validation will make
#' use even of a high number of compute nodes (up to number of folds times number of lambda values) efficiently.
#'
#' @inheritParams BinarySegmentation
#' @param n_folds Number of folds. Test data will be selected equi-spaced, i.e. each n_fold-th observation.
#' @param gamma_max If NULL the range of gamma will be determined by doing an extra fit on the full model and taking the
#' difference between the full segment loss and the loss of the first split.
#' @param verbose If TRUE additional information will be printed.
#'
#' @return A nested list with the cv results and the full fitted models for each gamm, lambda combination.
#' @export
#'
#'@examples
#' dat <- SimulateFromModel(CreateModel(n_segments = 2,n = 100,p = 30, ChainNetwork))
#' CrossValidation(dat, delta = 0.1, method = "summed_regression")
CrossValidation <- function(x, delta,
                            lambda = NULL,
                            method = c("nodewise_regression", "summed_regression", "ratio_regression"),
                            threshold = 1e-4,
                            penalize_diagonal = F,
                            use_ternary_search = F,
                            n_folds = 10,
                            gamma_max = NULL,
                            verbose = T,
                            ...) {
  n_obs <- nrow(x)
  n_p <- ncol(x)
  mth <- match.arg(method)

  # necessary because parser won't allow 'foreach' directly after a foreach object
  if (foreach::getDoParWorkers() == 1) {
    cat("\n No parallel backend registered. Cross-validation will be performed using a single node and might take very long. \n See for instance https://cran.r-project.org/web/packages/doParallel/index.html \n")
    `%hdcd_do%` <- foreach::`%do%`
  } else {
    `%hdcd_do%` <- foreach::`%dopar%`
  }
  `%:%` <- foreach::`%:%`

  # determine a range for gamma automatically
  if(is.null(gamma_max)){
    tree      <- BinarySegmentation(x = x, delta = delta, lambda = min(lambda), method = mth,
                                    penalize_diagonal = penalize_diagonal, use_ternary_search = use_ternary_search, ...)
    gamma_max <- tree$Get(function(x) abs(x$min_loss - x$segment_loss), filterFun = data.tree::isRoot)
    rm(tree)
  }

  # determine a range for lambda, not sure how to pick it in a smarter way yet...
  if(is.null(lambda)){
    lambda <- seq(0.01, 0.15, 0.01)
  }
  n_lambdas <- length(lambda)

  folds <- seq_len(n_folds)

  if (verbose) cat("\n")

  cv_results <- foreach::foreach(fold = folds, .inorder = F, .packages = "hdcd", .verbose = F) %:%
    foreach::foreach(lam = lambda, .inorder = T) %hdcd_do% {

      test_inds  <- seq(fold, n_obs, n_folds)
      train_inds <- setdiff(1:n_obs, test_inds)
      n_g <- length(test_inds)

      tree <- BinarySegmentation(x = x[train_inds, ], delta = delta, lambda = lam,
                                 method = mth, penalize_diagonal = penalize_diagonal,
                                 use_ternary_search = use_ternary_search, ...)
      res  <- PruneTreeGamma(tree, gamma_max)
      rm(tree)
      rss_gamma <- numeric(length(res$gamma))
      cpts      <- list()
      for (gam in seq_along(res$gamma)){
        fit <- FullRegression(x[train_inds, ], cpts = res$cpts[[gam]], lambda = lam)

        segment_bounds  <- c(1, train_inds[res$cpts[[gam]]], n_obs) # transform cpts back to original indices

        rss <- 0
        for (seg in seq_along(fit[[1]])){

          wi <- fit$est_coefs[[seg]]
          intercepts <- fit$est_intercepts[[seg]]

          seg_test_inds <- test_inds[which(test_inds >= segment_bounds[seg] & test_inds < segment_bounds[seg + 1])]

          if(length(seg_test_inds) == 0){warning("Segment had no test data. Consider reducing the number of folds."); next}

          rss <- rss +
            sum(sapply(1:n_p, function(z) RSS(x[seg_test_inds, -z], x[seg_test_inds, z, drop = F], wi[-z, z, drop = F], intercepts[z]))) / n_obs
        }
        rss_gamma[gam] <- rss / n_g
        cpts[[gam]]    <- segment_bounds[-c(1, length(segment_bounds))]
      }
      g <- res$gamma
      rm(res)
      if (verbose) cat(paste(Sys.time(),"  FINISHED fit for fold ", fold, " and lambda value ", lam, " \n"))
      list(rss = rss_gamma, cpts = cpts, gamma = g, lambda = lam, fold = fold)
    }
  n_gammas <- length(cv_results[[1]][[1]][["gamma"]])
  gamma_names  <- round(cv_results[[1]][[1]][["gamma"]], 3)

  res <- array(data = NA, dim = c(2, n_folds, n_lambdas, n_gammas),
               dimnames = list(type = c("rss", "n_cpts"), fold = seq_len(n_folds), lambda = lambda,
                               gamma = gamma_names))

  for (fold in seq_len(n_folds)){
    for (lam in seq_len(n_lambdas)){
      res["rss", fold, lam, ]    <- cv_results[[fold]][[lam]][["rss"]]
      res["n_cpts", fold, lam, ] <- unlist(lapply(X = cv_results[[fold]][[lam]][["cpts"]], FUN = length))
    }
  }

  # Crude rule to determine final model, should use something like 1SE rule instead!
  inds <- arrayInd(which.min(apply(res["rss", , ,], 2:3, mean)), .dim = c(n_lambdas, n_gammas))

  list(cv_results = res, best_gamma = gamma_names[inds[2]], best_lambda = lambda[inds[1]])
}

RSS <- function(x, y, beta, intercepts){
  sum((y - x %*% beta - intercepts)^2)
}

#' plot.bs_cv
#'
#' S3 method for plotting the results of cross-validation.
#'
#' @param cv_results An object of class \strong{bs_cv}
#'
#' @export
plot.bs_cv <- function(cv_results){

  col.l <- heat.colors(100)

  #lab <- dimnames(cv_results["rss", 1, ,])

  cpts <- lattice::contourplot(apply(cv_results["n_cpts", , ,], 2:3, mean),
                               aspect = "xy",
                               xlab = "Lambda",
                               ylab = "Gamma",
                               main = "Average number of changepoints",
                               col.regions = col.l,
                               cuts = 20,
                               region = T)
  rss <- lattice::contourplot(apply(cv_results["rss", , ,], 2:3, mean),
                              aspect = "xy",
                              xlab = "Lambda",
                              ylab = "Gamma",
                              main = "Average RSS",
                              col.regions = col.l,
                              cuts = 20,
                              region = T)

  print(cpts, split = c(1, 1, 2, 1), more=TRUE)
  print(rss, split = c(2, 1, 2, 1))
}
