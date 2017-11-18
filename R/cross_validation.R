#' CrossValidation
#'
#' EXPERIMENTAL. Crossvalidation for the desired method and parameter combinations.
#' Evalusate lambda values will refit the model whereas gamma values can be evaluted using a
#' single fit.
#'
#' Will make use of a registered parallel backend over the folds. This will be changed to parallelize over lambdas
#' in each fold to make use of all nodes if the number of folds is smaller than the number of nodes.
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
                            lambda = seq(0.01, 0.1, 0.01),
                            method = c("summed_regression", "ratio_regression"),
                            threshold = 1e-4,
                            penalize_diagonal = F,
                            use_ternary_search = F,
                            n_folds = 10,
                            gamma_max = NULL,
                            verbose = T) {
  n_obs <- nrow(x)
  n_p <- ncol(x)
  mth <- match.arg(method)
  `%dopar%` <- foreach::`%dopar%`

  if(is.null(gamma_max)){
    tree      <- BinarySegmentation(x = x, delta = delta, lambda = min(lambda), method = mth,
                                    penalize_diagonal = penalize_diagonal, use_ternary_search = use_ternary_search)
    gamma_max <- tree$Get(function(x) abs(x$min_loss - x$segment_loss), filterFun = data.tree::isRoot)
    rm(tree)
  }

  folds <- seq_len(n_folds)

  cv_results <- foreach::foreach(fold = folds, .inorder = F, .packages = "hdcd", .verbose = verbose) %dopar% {
    if (verbose) cat("\n CV fold ",fold, " of ",n_folds,"\n")

    test_inds  <- seq(fold, n_obs, n_folds)
    train_inds <- setdiff(1:n_obs, test_inds)
    n_g <- length(test_inds)
    fold_results <- list()

    for (lam in seq_along(lambda)){
      tree <- BinarySegmentation(x = x[train_inds, ], delta = delta, lambda = lambda[lam],
                                 method = mth, penalize_diagonal = penalize_diagonal, use_ternary_search = use_ternary_search)
      res  <- PruneTreeGamma(tree, gamma_max)

      rss_gamma <- numeric(length(res$gamma))
      cpts      <- list()
      for (gam in seq_along(res$gamma)){
        fit <- FullRegression(x[train_inds, ], cpts = res$cpts[[gam]], lambda = lambda[lam])

        segment_bounds  <- c(1, train_inds[res$cpts[[gam]]], n_obs) # transform cpts back to original indices

        rss <- 0
        for (seg in seq_along(fit[[1]])){

          wi <- fit$est_coeffs[[seg]]
          intercepts <- fit$est_intercepts[[seg]]

          segment_test_inds <- test_inds[which(test_inds >= segment_bounds[seg] & test_inds < segment_bounds[seg + 1])]

          if(length(segment_test_inds) == 0){warning("Segment had no test data. Consider reducing the number of folds."); next}

          rss <- rss +  sum(sapply(1:n_p, function(z) RSS(x[segment_test_inds, -z], x[segment_test_inds, z, drop = F], wi[-z, z, drop = F], intercepts[z]))) / n_obs
        }
        rss_gamma[gam] <- rss / n_g
        cpts[[gam]]    <- segment_bounds[-c(1, length(segment_bounds))]
      }
      fold_results[[lam]] <- list(rss = rss_gamma, gamma = res$gamma, cpts = cpts, tree = tree, lambda = lambda[lam])
    }
    fold_results
  }
  class(cv_results) <- "bs_cv"
  cv_results
}

RSS <- function(x, y, beta, intercepts){
  sum((y - x %*% beta - intercepts)^2)
}

#' plot.bs_cv
#'
#' S3 method for plotting the results of cross-validation. Ony works for a specific lambda value at the moment.
#'
#'EXPERIMENTAL
#'
#' @param cv_results An object of class \strong{bs_cv}
#'
#' @export
plot.bs_cv <- function(cv_results){

  n_cpts <- sapply(cv_results, function(x) sapply(x[["cpts"]], function(y) length(y)))
  gamma  <- cv_results[[1]][["gamma"]]
  rss <- rowMeans(sapply(cv_results, function(x) x[["rss"]]))
  sd  <- apply(sapply(cv_results, function(x) x[["rss"]]), 1, sd)

  par(mar = c(5, 4, 4, 4) + 0.3)  # Leave space for z axis
  plot (gamma, rss, ylim = c(0, 1.1*max(rss+sd)), col = "red", pch = 18, ylab = "rss / n_obs")
  segments(gamma,rss-sd,gamma,rss+sd)
  epsilon <- 0.0005
  segments(gamma-epsilon,rss-sd,gamma+epsilon,rss-sd)
  segments(gamma-epsilon,rss+sd,gamma+epsilon,rss+sd)
  par(new = TRUE)
  plot(gamma, n_cpts[, 1], type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "", col = rgb(0,0,1,alpha=0.2))
  for (i in 2:ncol(n_cpts)){
    lines(gamma, n_cpts[, i], col = rgb(0,0,1,alpha=0.2))
  }
  axis(side=4, at = pretty(range(n_cpts[, 1])))
  mtext("n_cpts", side=4, line=3)
}
