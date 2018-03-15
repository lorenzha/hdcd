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
#' @param gamma If NULL the range of gamma will be determined by doing an extra fit on the full model and taking the
#' difference between the full segment loss and the loss of the first split.
#' @param verbose If TRUE additional information will be printed.
#'
#' @return A nested list with the cv results and the full fitted models for each gamm, lambda combination.
#' @export
#'
#' @examples
#' dat <- SimulateFromModel(CreateModel(n_segments = 2,n = 100,p = 30, ChainNetwork))
#' CrossValidation(dat, delta = 0.1, method = "summed_regression")
CrossValidation <- function(x,
                            delta = c(0.1, 0.25),
                            lambda = NULL,
                            lambda_gamma = NULL,
                            gamma = NULL,
                            n_folds = 10,
                            grid_size = 20,
                            method = c("nodewise_regression", "summed_regression", "ratio_regression"),
                            penalize_diagonal = F,
                            optimizer = c("line_search", "ternary_search", "section_search"),
                            control = NULL,
                            standardize = T,
                            threshold = 1e-7,
                            parallel = T,
                            verbose = T,
                            ...) {
  n_obs <- nrow(x)
  n_p <- ncol(x)

  # necessary because parser won't allow 'foreach' directly after a foreach object
  if (foreach::getDoParWorkers() == 1 && parallel) {
    cat("\nNo parallel backend registered. \n\nCross-validation will be performed using a single compute node and might take a very long time. See for instance https://cran.r-project.org/web/packages/doParallel/index.html to install a parallel backend.\n")
    `%hdcd_do%` <- foreach::`%do%`
  } else if (!parallel) {
    `%hdcd_do%` <- foreach::`%do%`
  } else {
    `%hdcd_do%` <- foreach::`%dopar%`
  }
  `%:%` <- foreach::`%:%`

  # choose lambda as grid around the asymptotic value
  if (is.null(lambda)) {
    lambda <- seq(0.05, 1.5, length.out = grid_size) * sqrt(log(n_p) / n_obs)
  }
  n_lambdas <- length(lambda)

  # choose three sensible values for delta
  if (is.null(delta)) {
    delta <- c(0.05, 0.1, 0.25)
  }
  n_delta <- length(delta)

  # Take smallest delta and specific or median lambda to have the broadest range of the loss
  if (is.null(gamma)) {
    l <- ifelse(!is.null(lambda_gamma), lambda_gamma, median(lambda))
    tree <- BinarySegmentation(
      x = x, delta = min(delta), lambda = l,
      method = method, penalize_diagonal = penalize_diagonal,
      optimizer = optimizer, control = control,
      threshold = threshold, standardize = standardize, ...
    )
    gamma_diff <- abs(tree$Get("segment_loss") - tree$Get("min_loss"))
    gamma <- seq(0, max(gamma_diff, na.rm = T), length.out = grid_size)
  }
  n_gammas <- length(gamma)

  if (verbose) cat("\n")
  cv_results <- foreach::foreach(fold = seq_len(n_folds), .inorder = F, .packages = "hdcd", .verbose = F) %:%
    foreach::foreach(del = delta, .inorder = T) %:%
    foreach::foreach(lam = lambda, .inorder = T) %hdcd_do% {
      test_inds <- seq(fold, n_obs, n_folds)
      train_inds <- setdiff(1:n_obs, test_inds)
      n_g <- length(test_inds)

      tree <- BinarySegmentation(
        x = x[train_inds, ], delta = del, lambda = lam,
        method = method, penalize_diagonal = penalize_diagonal,
        optimizer = optimizer, control = control, threshold = threshold,
        standardize = standardize, ...
      )

      res <- PruneTreeGamma(tree, gamma)
      rm(tree)
      rss_gamma <- numeric(length(gamma))
      cpts <- list()
      for (gam in seq_along(gamma)) {
        fit <- FullRegression(
          x[train_inds, ], cpts = res$cpts[[gam]], # TODO: Can we somehow cache the fits from before instead of refitting the model? Should be the endpoints of the pruned tree!
          lambda = lam, standardize = standardize,
          threshold = threshold
        )

        segment_bounds <- c(1, train_inds[res$cpts[[gam]]], n_obs) # transform cpts back to original indices

        rss <- 0
        for (seg in seq_along(fit[[1]])) {
          wi <- fit$est_coefs[[seg]]
          intercepts <- fit$est_intercepts[[seg]]

          seg_test_inds <- test_inds[which(test_inds >= segment_bounds[seg] & test_inds < segment_bounds[seg + 1])]

          if (length(seg_test_inds) == 0) {
            warning("Segment had no test data. Consider reducing the number of folds.")
            next
          }
          # TODO: Instead of calculating the RSS, we could take the likelihood ratio here again, can we store the loss for one segment per
          # dimension before and reuse it here?
          rss <- rss +
            sum(sapply(1:n_p, function(z) RSS(x[seg_test_inds, -z], x[seg_test_inds, z, drop = F], wi[-z, z, drop = F], intercepts[z]))) / n_obs
        }
        rss_gamma[gam] <- rss / n_g
        cpts[[gam]] <- segment_bounds[-c(1, length(segment_bounds))]
      }
      rm(res)
      if (verbose) cat(paste(Sys.time(), "  FINISHED fit -  Fold: ", fold, " Lambda: ", round(lam, 3), " Delta: ", round(del, 3), " \n"))
      list(rss = rss_gamma, cpts = cpts)
    }


  res <- array(
    data = NA, dim = c(2, n_folds, n_delta, n_lambdas, n_gammas),
    dimnames = list(
      type = c("rss", "n_cpts"), fold = seq_len(n_folds),
      delta = delta, lambda = lambda, gamma = gamma
    )
  )

  for (fold in seq_len(n_folds)) {
    for (lam in seq_len(n_lambdas)) {
      for (del in seq_len(n_delta)) {
        res["rss", fold, del, lam, ] <- cv_results[[fold]][[del]][[lam]][["rss"]]
        res["n_cpts", fold, del, lam, ] <- unlist(lapply(X = cv_results[[fold]][[del]][[lam]][["cpts"]], FUN = length))
      }
    }
  }

  # Crude rule to determine final model
  inds <- arrayInd(which.min(apply(res["rss", , , , , drop = FALSE], 3:5, mean)), .dim = c(n_delta, n_lambdas, n_gammas))

  list(cv_results = res, best_delta = delta[inds[1]], best_gamma = gamma[inds[3]], best_lambda = lambda[inds[2]])
}

RSS <- function(x, y, beta, intercepts) {
  sum((y - x %*% beta - intercepts) ^ 2)
}

#' plot.bs_cv
#'
#' S3 method for plotting the results of cross-validation.
#'
#' @param cv_results An object of class \strong{bs_cv}
#'
#' @export
plot.bs_cv <- function(results) {
  res <- results[["cv_results"]]

  gam_names <- round(as.numeric(dimnames(res)$gamma),3)
  lam_names <- round(as.numeric(dimnames(res)$lambda),3)
  del_names <- round(as.numeric(dimnames(res)$delta),3)

  n_delta <- length(del_names)

  for (delta in seq_len(n_delta)){

    if(delta == 1) {
      main_cpts <-  "Average number of changepoints"
      main_rss <-  "Average RSS"
    } else {
      main_cpts <-  NULL
      main_rss <- NULL
    }

    cpts <- lattice::contourplot(
      apply(res["n_cpts", , delta, , ], 2:3, mean),
      aspect = "xy",
      xlab = "Lambda",
      xlab.top = paste("Delta:", del_names[delta]),
      ylab = "Gamma",
      main = main_cpts,
      col.regions = rainbow(20),
      cuts = 20,
      region = T,
      scales = list(y = list(label = gam_names),
                    x = list(label = lam_names, rot = 90))
    )
    rss <- lattice::contourplot(
      apply(res["rss", , delta, , ], 2:3, mean),
      aspect = "xy",
      xlab = "Lambda",
      xlab.top = paste("Delta:", del_names[delta]),
      ylab = "Gamma",
      main = main_rss,
      col.regions = rev(heat.colors(100)),
      cuts = 20,
      region = T,
      scales = list(y = list(label = gam_names),
                    x = list(label = lam_names, rot = 90))
    )

    if (delta == n_delta) m <- FALSE else m <- TRUE
    print(cpts, split = c(1 , delta, 2, n_delta), more = TRUE)
    print(rss, split = c(2, delta, 2, n_delta),  more = m)
  }
}
