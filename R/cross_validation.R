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
                            lambda_min = 0.001,
                            lambda_grid_size = 10,
                            gamma = NULL,
                            n_folds = 10,
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
    cov_mat <- cov(x)
    lambda_max <- max(abs(cov_mat[upper.tri(cov_mat)]))
    lambda <- LogSpace(lambda_min, lambda_max, length.out = lambda_grid_size)
  }
  n_lambdas <- length(lambda)

  # choose three sensible values for delta
  if (is.null(delta)) {
    delta <- c(0.05, 0.1, 0.25)
  }
  n_delta <- length(delta)

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

      gamma <- c(0, sort(tree$Get("segment_loss") - tree$Get("min_loss")))
      gamma <- gamma[which(gamma >= 0)]

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
        cpts[[gam]] <- segment_bounds
      }
      rm(res)
      if (verbose) cat(paste(Sys.time(), "  FINISHED fit -  Fold: ", fold, " Lambda: ", round(lam, 3), " Delta: ", round(del, 3), " \n"))
      list(fold = fold, lambda = lam, delta = del, gamma = gamma, rss = rss_gamma, cpts = cpts)
    }

  res <- data.frame()
  for (fold in seq_len(n_folds)) {
    for (lam in seq_len(n_lambdas)) {
      for (del in seq_len(n_delta)) {
        r <- cv_results[[fold]][[del]][[lam]]
        new_res <- data.frame(fold = r[["fold"]],
                              lambda = r[["lambda"]],
                              delta = r[["delta"]],
                              gamma = r[["gamma"]],
                              rss = r[["rss"]])
        # Need to assign those separately since it is a list
        new_res$cpts <- r[["cpts"]]
        res <- rbind(res, new_res)
      }
    }
  }

  single <- split(res, f = list(res$lambda, res$delta), drop = T)

  results <- lapply(single, SolutionPaths)

  opt <- do.call(rbind,lapply(results, GetOpt))

  list(opt = opt[which.min(opt$rss),, drop = TRUE], cv_results = results)
}


#' plot.bs_cv
#'
#' S3 method for plotting the results of cross-validation.
#'
#' @param cv_results An object of class \strong{bs_cv}
#'
#' @export
plot.bs_cv <- function(results, ...) {
  res <- results[["cv_results"]]


  res_long <- do.call(rbind,
                      lapply(res,
                             function(x) if (nrow(x[["cpts"]]) == 1)
                               data.frame(delta = x[["delta"]],
                                          lambda = x[["lambda"]],
                                          gamma = x[["rss"]][,1],
                                          rss = rowMeans(x[["rss"]][,-1]),
                                          n_cpts = 0,
                                          key = paste(x[["delta"]], formatC(x[["lambda"]], format = "e", digits = 2)))
                             else
                              data.frame(delta = x[["delta"]],
                                                    lambda = x[["lambda"]],
                                                    gamma = x[["rss"]][,1],
                                                    rss = rowMeans(x[["rss"]][,-1]),
                                                    n_cpts = rowMeans(apply(x[["cpts"]][,-1, drop = F], 2, function(x) sapply(x, length) - 2)), #substract first and last segment boundary
                                                    key = paste(x[["delta"]], formatC(x[["lambda"]], format = "e", digits = 2)))))
  levs <- levels(res_long$key)
  cols <- rainbow(length(levs))

  op <- par(mfrow = c(2,1), mar = c(4,5,4,2))
  plot(res_long$gamma, res_long$rss, type = "n", axes = F, ylab = "Average loss", xlab ="Gamma", ...)
  axis(side = 2)
  axis(side = 1)
  box()
  for (i in seq_along(levs)){
    plot_dat <- res_long[res_long$key == levs[i], ]
    lines(plot_dat[, "gamma"], plot_dat[, "rss"], col = cols[i], type = "s", xlim = c(min(plot_dat[, "gamma"]), max(plot_dat[, "gamma"])))
  }

  par(mar = c(6,5,2,2))
  plot(res_long$gamma, res_long$n_cpts, type = "n", axes = F, ylab = "Average # of changepoints", xlab ="", ...)
  axis(side = 2)
  axis(side = 1)
  box()
  for (i in seq_along(levs)){
    plot_dat <- res_long[res_long$key == levs[i], ]
    lines(plot_dat[, "gamma"], plot_dat[, "n_cpts"], col = cols[i], type = "s")
  }
  par(fig = c(0, 1, 0, 1), mar = c(0,0,0,0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("bottom", legend = levels(factor(formatC(res_long$lambda, format = "e", digits = 2))), xpd = TRUE, ncol = ceiling(length(levs)/3),
         lty = 1, bty = "n", col = cols, cex = 0.8)
}


SolutionPaths <- function(dat){
  dat <- dat[order(dat$gamma), ]
  list(lambda = dat$lambda[1], delta = dat$delta[1],
       rss = ImputeMatrix(reshape2::dcast(dat, gamma ~ fold, value.var = "rss")),
       cpts = ImputeMatrix(reshape2::dcast(dat, gamma ~ fold, value.var = "cpts", fill = NA)))
}

ImputeMatrix <- function(mat, cols = seq(2, ncol(mat))){
  temp_mat <- mat[,cols]
  for (i in seq(2, nrow(temp_mat))) {
    for (j in seq_len(ncol(temp_mat))){
      if (anyNA(temp_mat[i,j][[1]]) | is.null(temp_mat[i,j][[1]]))
        temp_mat[i,j][[1]] <- ifelse((is.null(temp_mat[i - 1,j]) | length(temp_mat[i - 1,j]) < 1), NA, temp_mat[i - 1,j])

    }
  }
  mat[, cols] <- temp_mat
  mat
}

RSS <- function(x, y, beta, intercepts) {
  sum((y - x %*% beta - intercepts) ^ 2)
}

GetOpt <- function(param_res){

  avg_rss <- rowMeans(param_res$rss[,-1])
  opt <- which.min(avg_rss)

  data.frame(lambda = param_res$lambda,
             delta = param_res$delta,
             gamma = param_res$rss[opt, 1],
             rss = avg_rss[opt])

}
