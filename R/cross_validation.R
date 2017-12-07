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
CrossValidation <- function(x,
                            delta = NULL,
                            lambda = NULL,
                            gamma = NULL,
                            n_folds = 10,
                            grid_size = 20,
                            method = c("nodewise_regression", "summed_regression", "ratio_regression"),
                            threshold = 1e-4,
                            penalize_diagonal = F,
                            use_ternary_search = F,
                            verbose = T,
                            parallel = T,
                            ...) {
  n_obs <- nrow(x)
  n_p <- ncol(x)
  mth <- match.arg(method)

  # necessary because parser won't allow 'foreach' directly after a foreach object
  if (foreach::getDoParWorkers() == 1 && parallel) {
    cat("\n No parallel backend registered. \n\n Cross-validation will be performed using a single node and might take very long. See for instance https://cran.r-project.org/web/packages/doParallel/index.html to install a parallel backend.\n")
    `%hdcd_do%` <- foreach::`%do%`
  } else if (!parallel){
    `%hdcd_do%` <- foreach::`%do%`
  } else {
    `%hdcd_do%` <- foreach::`%dopar%`
  }
  `%:%` <- foreach::`%:%`

  # choose lambda as grid around the asymptotic value
  if(is.null(lambda)){
    lambda <- seq(0.05, 1.5, length.out = grid_size) * sqrt(log(n_p)/n_obs)
  }
  n_lambdas <- length(lambda)

  if (verbose) cat("\n")
  cv_results <- foreach::foreach(fold = seq_len(n_folds), .inorder = F, .packages = "hdcd", .verbose = F) %:%
    foreach::foreach(lam = lambda, .inorder = T) %hdcd_do% {

      test_inds  <- seq(fold, n_obs, n_folds)
      train_inds <- setdiff(1:n_obs, test_inds)
      n_g <- length(test_inds)

      tree <- BinarySegmentation(x = x[train_inds, ], delta = delta, lambda = lam,
                                 method = mth, penalize_diagonal = penalize_diagonal,
                                 use_ternary_search = use_ternary_search, ...)
      if (is.null(gamma)){
        gamma_diff <- abs(tree$Get("segment_loss") - tree$Get("min_loss"))
        gamma  <- seq(min(gamma_diff, na.rm = T), max(gamma_diff, na.rm = T), length.out = grid_size)
      }

      res  <- PruneTreeGamma(tree, gamma)
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
      if (verbose) cat(paste(Sys.time(),"  FINISHED fit for fold ", fold, " and lambda value ", round(lam, 2), " \n"))
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

  # Crude rule to determine final model
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
plot.bs_cv <- function(results){

  res <- results[["cv_results"]]

  cpts <- lattice::contourplot(apply(res["n_cpts", , ,], 2:3, mean),
                               aspect = "xy",
                               xlab = "Lambda",
                               ylab = "Gamma",
                               main = "Average number of changepoints",
                               col.regions = rainbow(20),
                               cuts = 20,
                               region = T)
  rss <- lattice::contourplot(apply(res["rss", , ,], 2:3, mean),
                              aspect = "xy",
                              xlab = "Lambda",
                              ylab = "Gamma",
                              main = "Average RSS",
                              col.regions = rev(heat.colors(100)),
                              cuts = 20,
                              region = T)

  print(cpts, split = c(1, 1, 2, 1), more=TRUE)
  print(rss, split = c(2, 1, 2, 1))
}
