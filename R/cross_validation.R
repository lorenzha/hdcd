#' CrossValidation
#'
#' Cross-validation for the desired method and parameter combinations.
#'
#' Evaluating different lambda values will lead to refitting the entire model whereas gamma values can be evaluted cheaply using a
#' single fit. Suitable values for lambda as well as gamma are choosen automatically if they are not supplied.
#'
#' Will make use of a registered parallel backend over the folds and lambda values. Therefore the whole cross-validation will make
#' use even of a high number of compute nodes (up to number of folds times number of lambda values).
#'
#' @inheritParams BinarySegmentation
#' @param n_folds Number of folds. Test data will be selected equi-spaced, i.e. each n_fold-th observation.
#' @param lambda_min_ratio Numeric value between 0 and 1. If the \eqn{\lambda_max} is determined internally this will pick \eqn{\lambda_min = lambda_min_ratio * \lambda_max}.
#' @param lambda_grid_size Integer value determining the number of values between \eqn{\lambda_min} and \eqn{\lambda_max} to will be equally spaced on a logarithmic scale.
#' @param gamma Numeric value or vector. If NULL the full solution path for gamma will be caluclated for every combination of \eqn{\lambda} and \eqn{\delta}
#' @param parallel If TRUE and a parallel backend is registered, the cross-validation will be performed in parallel.
#' @param verbose Boolean. If TRUE additional information will be printed.
#'
#' @return A nested list with the cv results and the full fitted models for each combination of \eqn{\delta}, \eqn{lambda} and \eqn{gamma} combination.
#' @export
#'
#' @importFrom stats cov
#'
#' @examples
#' \dontrun{
#' dat <- SimulateFromModel(CreateModel(n_segments = 2,n = 100,p = 30, ChainNetwork))
#' CrossValidation(dat, method = "summed_regression")
#' }
CrossValidation <- function(x,
                            delta = c(0.1, 0.25),
                            lambda = NULL,
                            lambda_min_ratio = 0.01,
                            lambda_grid_size = 10,
                            gamma = NULL,
                            n_folds = 10,
                            method = c("nodewise_regression", "summed_regression", "ratio_regression"),
                            penalize_diagonal = F,
                            optimizer = c("line_search", "section_search"),
                            control = NULL,
                            standardize = T,
                            threshold = 1e-7,
                            parallel = T,
                            verbose = T,
                            FUN = NULL,
                            ...) {
  n_obs <- NROW(x)
  n_p <- NCOL(x)

  if (is.null(FUN)){
    SegmentLossFUN <- SegmentLoss(
      n_obs = NROW(x), lambda = lambda, penalize_diagonal = penalize_diagonal,
      method = method, standardize = standardize, threshold = threshold, ...
    )
  } else {
    stopifnot(c("x", "n_obs", "standardize") %in% formalArgs(FUN))
    SegmentLossFUN <- functional::Curry(FUN, n_obs = NROW(x), standardize = standardize)
  }


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
  if (is.null(lambda) && is.null(FUN) && NCOL(x) > 1) {
    cov_mat <- cov(x)
    lambda_max <- max(abs(cov_mat[upper.tri(cov_mat)]))
    lambda <- LogSpace(lambda_min_ratio * lambda_max, lambda_max, length.out = lambda_grid_size)
  } else if (!is.null(FUN)){
    lambda <- c(1, 2)
  }

  n_lambdas <- length(lambda)

  # choose three sensible values for delta
  if (is.null(delta)) {
    delta <- c(0.05, 0.1, 0.2)
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
        x = x[train_inds, , drop = F], delta = del, lambda = lam,
        method = method, penalize_diagonal = penalize_diagonal,
        optimizer = optimizer, control = control, threshold = threshold,
        standardize = standardize, FUN = FUN, ...
      )

      if (is.null(gamma)) {
        final_gamma <- c(0, sort(tree$Get("segment_loss") - tree$Get("min_loss")))
        final_gamma <- final_gamma[which(final_gamma >= 0)]
      } else {
        final_gamma <- gamma
      }

      res <- PruneTreeGamma(tree, final_gamma)
      rm(tree)
      rss_gamma <- numeric(length(final_gamma))
      n_params_gamma <- numeric(length(final_gamma))
      cpts <- list()
      for (gam in seq_along(final_gamma)) {
        fit <- FullRegression(
          x[train_inds, , drop =F], cpts = res$cpts[[gam]], # TODO: Can we somehow cache the fits from before instead of refitting the model? Should be the endpoints of the pruned tree!
          lambda = lam, standardize = standardize, threshold = threshold
        )

        segment_bounds <- c(1, train_inds[res$cpts[[gam]]], n_obs) # transform cpts back to original indices

        rss <- 0
        n_params <- 0

        for (seg in seq_len(length(segment_bounds) - 1)) {

          seg_test_inds <- test_inds[which(test_inds >= segment_bounds[seg] & test_inds < segment_bounds[seg + 1])]

          if (length(seg_test_inds) == 0) {
            warning("Segment had no test data. Consider reducing the number of folds.")
            next
          }
          # TODO: Instead of calculating the RSS, we could take the likelihood ratio here again, can we store the loss for one segment per
          #  dimension before and reuse it here?

          if(n_p > 1) {
            loss <- loss +
              sum(sapply(1:n_p, function(z) RSS(x[seg_test_inds, -z, drop = F], x[seg_test_inds, z, drop = F], wi[-z, z, drop = F], intercepts[z]))) / n_obs
            n_params <- n_params + length(which(wi[upper.tri(wi)] != 0))
          } else {
            loss <- loss + RSS(x[seg_test_inds, , drop = F], x[seg_test_inds, , drop = F], wi, intercepts) / n_obs
          }
        }
        rss_gamma[gam] <- rss / n_g
        n_params_gamma[gam] <- n_params
        cpts[[gam]] <- segment_bounds
      }
      rm(res)
      if (verbose) cat(paste(Sys.time(), "  FINISHED fit -  Fold: ", fold, " Lambda: ", round(lam, 3), " Delta: ", round(del, 3), " \n"))
      list(fold = fold, lambda = lam, delta = del, gamma = final_gamma, rss = rss_gamma, cpts = cpts, n_params = n_params_gamma)
    }

  res <- data.frame()
  for (fold in seq_len(n_folds)) {
    for (lam in seq_len(n_lambdas)) {
      for (del in seq_len(n_delta)) {
        r <- cv_results[[fold]][[del]][[lam]]
        new_res <- data.frame(
          fold = r[["fold"]],
          lambda = r[["lambda"]],
          delta = r[["delta"]],
          gamma = r[["gamma"]],
          rss = r[["rss"]],
          n_params = r[["n_params"]]
        )
        # Need to assign those separately since it is a list
        new_res$cpts <- r[["cpts"]]
        res <- rbind(res, new_res)
      }
    }
  }

  single <- split(res, f = list(res$lambda, res$delta), drop = T)

  results <- lapply(single, SolutionPaths)

  opt <- do.call(rbind, lapply(results, GetOpt))

  list(opt = opt[which.min(opt$rss), , drop = TRUE], cv_results = results)
}


#' plot.bs_cv
#'
#' S3 method for plotting the results of cross-validation.
#'
#' @param x An object of class \strong{bs_cv}
#' @param ... Only included to be consist with plot generic.
#' @param show_legend If TRUE the legend is shown on the bottom of the plot.
#'
#' @importFrom grDevices rainbow
#' @importFrom reshape2 melt
#' @importFrom latex2exp TeX
#' @import ggplot2
#'
#' @export
plot.bs_cv <- function(x, ..., show_legend = T) {
  res_long <- do.call(
    rbind,
    lapply(
      x[["cv_results"]],
      function(x) if (nrow(x[["cpts"]]) == 1) {
          data.frame(
            delta = factor(x[["delta"]]),
            lambda = formatC(x[["lambda"]], format = "e", digits = 2),
            gamma = x[["rss"]][, 1],
            rss = rowMeans(x[["rss"]][, -1]),
            n_params = rowMeans(x[["n_params"]][, -1]),
            n_cpts = 0
          )
        } else {
          data.frame(
            delta = factor(x[["delta"]]),
            lambda = formatC(x[["lambda"]], format = "e", digits = 2),
            gamma = x[["rss"]][, 1],
            rss = rowMeans(x[["rss"]][, -1]),
            n_params = rowMeans(x[["n_params"]][, -1]),
            n_cpts = rowMeans(apply(x[["cpts"]][, -1, drop = F], 2, function(x) sapply(x, length) - 2))
          )
        }
    )
  ) # substract first and last segment boundary

  res_long <- reshape2::melt(res_long, measure.vars = c("rss", "n_cpts", "n_params"), value.name = "value", variable.name = "metric")

  metrics_names <- c(
    "n_params" = "# of parameters",
    "rss" = "cv-loss",
    "n_cpts" = "# of changepoints"
  )

  if (show_legend) {
    l_pos <- "bottom"
  } else {
    l_pos <- "none"
  }

  ggplot(res_long) +
    geom_line(aes_(x = ~gamma, y = ~value, color = ~lambda, linetype = ~delta)) +
    ylab("") +
    xlab(latex2exp::TeX("$\\gamma$")) +
    theme_minimal() +
    theme(
      legend.position = l_pos,
      axis.text = element_text(size = 14, colour = "black"),
      strip.text = element_text(size = 14, colour = "black"),
      legend.text = element_text(size = 9),
      axis.title = element_text(size = 14),
      strip.background = element_rect(colour = NA, fill = NA),
      strip.placement = "outside"
    ) +
    facet_grid(metric~., scales = "free_y", switch = "y", labeller = labeller(metric = metrics_names)) +
    scale_colour_manual(
      values = grDevices::rainbow(length(levels(res_long$lambda))),
      name = latex2exp::TeX("$\\lambda$"),
      breaks = levels(res_long$lambda),
      labels = levels(res_long$lambda)
    ) +
    scale_linetype_manual(
      values = seq(1:length(levels(res_long$delta))),
      name = latex2exp::TeX("$\\delta$"),
      breaks = levels(res_long$delta),
      labels = levels(res_long$delta)
    )
}


SolutionPaths <- function(dat) {
  dat <- dat[order(dat$gamma), ]
  list(
    lambda = dat$lambda[1], delta = dat$delta[1],
    rss = ImputeMatrix(reshape2::dcast(dat, gamma ~ fold, value.var = "rss")),
    n_params = ImputeMatrix(reshape2::dcast(dat, gamma ~ fold, value.var = "n_params")),
    cpts = ImputeMatrix(reshape2::dcast(dat, gamma ~ fold, value.var = "cpts", fill = NA))
  )
}

ImputeMatrix <- function(mat, cols = seq(2, ncol(mat))) {
  temp_mat <- mat[, cols]
  for (i in seq(2, nrow(temp_mat))) {
    for (j in seq_len(ncol(temp_mat))) {
      if (anyNA(temp_mat[i, j][[1]]) | is.null(temp_mat[i, j][[1]])) {
        temp_mat[i, j][[1]] <- ifelse((is.null(temp_mat[i - 1, j]) | length(temp_mat[i - 1, j]) < 1), NA, temp_mat[i - 1, j])
      }
    }
  }
  mat[, cols] <- temp_mat
  mat
}

RSS <- function(x, y, beta, intercepts) {
  sum((y - x %*% beta - intercepts) ^ 2)
}

GetOpt <- function(param_res) {
  avg_rss <- rowMeans(param_res$rss[, -1])
  opt <- which.min(avg_rss)

  data.frame(
    lambda = param_res$lambda,
    delta = param_res$delta,
    gamma = param_res$rss[opt, 1],
    rss = avg_rss[opt]
  )
}

LogSpace <- function(from, to, length.out) {
  exp(seq(log(from), log(to), length.out = length.out))
}
