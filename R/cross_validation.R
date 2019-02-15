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
                            delta = NULL,
                            lambda = NULL,
                            lambda_min_ratio = 0.01,
                            lambda_grid_size = 10,
                            gamma = NULL,
                            n_folds = 10,
                            method = c("nodewise_regression", "summed_regression", "ratio_regression", 'glasso', 'elastic_net'),
                            NA_method = c('complete_observations', 'pairwise_covariance_estimation', 'loh_wainwright_bias_correction'),
                            penalize_diagonal = F,
                            alpha = 1,
                            optimizer = c("line_search", "section_search"),
                            control = NULL,
                            standardize = T,
                            threshold = 1e-7,
                            parallel = T,
                            verbose = T,
                            FUN = NULL,
                            cv.FUN = NULL,
                            ...) {

  n_obs <- NROW(x)
  n_p <- NCOL(x)

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
  if (is.null(lambda) && NCOL(x) > 1) {
    cov_mat <- get_cov_mat(x, NA_method)$mat
    lambda_max <- max(abs(cov_mat[upper.tri(cov_mat)]))
    lambda <- LogSpace(lambda_min_ratio * lambda_max, lambda_max, length.out = lambda_grid_size)
  }

  n_lambdas <- length(lambda)

  # choose three sensible values for delta
  if (is.null(delta)) {
    delta <- c(0.05, 0.1, 0.2)
  }
  n_delta <- length(delta)

  if (verbose) cat("\n")
  # for (lam in lambda){
  #   for (del in delta){
  #     for (fold in seq_len(n_folds)){
  cv_results <- foreach::foreach(fold = seq_len(n_folds), .inorder = F, .packages = "hdcd", .verbose = F) %:%
    foreach::foreach(del = delta, .inorder = T) %:%
    foreach::foreach(lam = lambda, .inorder = T) %hdcd_do% {
      test_inds <- seq(fold, n_obs, n_folds)
      train_inds <- setdiff(1:n_obs, test_inds)

      tree <- BinarySegmentation(
        x = x[train_inds, , drop = F], delta = del, lambda = lam,
        method = method, NA_method = NA_method, penalize_diagonal = penalize_diagonal,
        optimizer = optimizer, control = control, threshold = threshold,
        standardize = standardize, FUN = FUN, ...
      )

      if (is.null(gamma)) {
        final_gamma <- c(0, sort(tree$Get("max_gain")))
        final_gamma <- final_gamma[which(final_gamma >= 0)]
      } else {
        final_gamma <- gamma
      }

      res <- PruneTreeGamma(tree, final_gamma)
      rm(tree)
      loss_gamma <- numeric(length(final_gamma))
      cpts <- list()

      for (gam in seq_along(final_gamma)) {
        loss <- 0
        alpha <- c(1, train_inds[res$cpts[[gam]]], n_obs + 1) # transform cpts back to original indices
        for (i in 1:(length(alpha) - 1)){

          train_range <- intersect(train_inds, alpha[i] : (alpha[i + 1] - 1))
          test_range <- intersect(test_inds, alpha[i] : (alpha[i + 1] - 1))

          if (length(test_range) == 0){
            warning("Segment has no test data. Consider reducing the number of folds.")
          } else {
            loss <- loss + CrossValidationLoss(x_train = x[train_range, , drop = F], x_test = x[test_range, , drop = F], n_obs_train = length(train_inds),
                                   lambda = lam, penalize_diagonal = penalize_diagonal, threshold = threshold,
                                   method = method, NA_method = NA_method)
          }
        }

        loss_gamma[gam] <- loss / length(test_inds)
        cpts[[gam]] <- train_inds[res$cpts[[gam]]]
        if(verbose) cat(paste('Changepoints corresponding to Gamma = ', final_gamma[gam], ' are ', paste(cpts[[gam]], collapse = ', '), 'corresponding to loss = ', loss, '\n'))
      }
      rm(res)
      if (verbose) cat(paste(Sys.time(), "  FINISHED fit -  Fold: ", fold, " Lambda: ", round(lam, 3), " Delta: ", round(del, 3), " \n"))


      list(fold = fold, lambda = lam, delta = del, gamma = final_gamma, loss = loss_gamma, cpts = cpts)
    #   }
    # }
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
          loss = r[["loss"]]
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

  list(opt = opt[which.min(opt$loss), , drop = TRUE], cv_results = results)
}

#'
#'
#'
#
cv_hdcd <- function(x, y = NULL, method = "glasso", NA_method = "complete_observations",
                    optimizer = "line_search", delta = NULL, lambda = NULL, gamma = NULL, control = NULL){

  # parse values NA_method, method and optimizer
  mth <- match.arg(method, c("glasso", "nodewise_regression", "summed_regression", "ratio_regression", "elastic_net"))
  NA_mth <- match.arg(NA_method, c("complete_observations", "average_imputation", "loh_wainwright_bias_correction", "pairwise_covariance_estimation"))
  opt <- match.arg(optimizer, c("line_search", "section_search"))

  n <- nrow(x)

  # read parameters from control
  n_folds_outer <- control_get(control, "n_folds_outer", 5)
  randomize_outer_folds <- control_get(control, "randomize_outer_folds", FALSE)
  verbose <- control_get(control, "verbose", TRUE)

  # Series of checks to ensure functionality in special cases and return warnings
  if(!is.matrix(x)){
    x <- as.matrix(x)
    warning("Input data x has been coerced to matrix by hdcd.")
  }
  if(!is.null(y) & method != "elastic_net"){
    warning("Input y is ignored since method is not elastic_net")
  }
  # if (!is.null(node) & method != "nodewise_regression") {
  #   warning("Input node is ignored since method is not nodewise_regression")
  # }
  # if (!is.null(alpha) & method != "elastic_net"){
  #   warning("Input alpha is ignored since method is not elastic_net")
  # }
  if (NA_method == "complete_observations" & mean(complete.cases(x)) < 0.5){
    warning("Less than 50% of observations are complete. Consider using a different NA_method")
  }

  if (!is.null(y) & method == "elastic_net"){
    x <- cbind(y,x)
  }

  # choose lambda as grid around the asymptotic value if no specific value is supplied
  if (is.null(lambda) && NCOL(x) > 1) {
    cov_mat <- get_cov_mat(x, NA_method)$mat
    lambda_max <- max(abs(cov_mat[upper.tri(cov_mat)]))
    lambda <- LogSpace(0.01 * lambda_max, lambda_max, length.out = 10)
    if (verbose) cat("Values for lambda chosen by asymptotic theory are", round(lambda, 3), "\n", sep = " ")
  }

  # choose three sensible values for delta if no specific value is supplied
  if (is.null(delta)) {
    delta <- c(0.05, 0.1, 0.2)
    if (verbose) cat("Values for delta chosen are", delta, "\n", sep = " ")
  }

  folds_outer <- sample_folds(n, n_folds_outer, randomize_outer_folds)

  # set up data.table to save results in
  cv_results <- data.table::data.table()

  # do outer cross validation
  for (outer in 1:n_folds_outer){
    for (lam in lambda){
      for (del in delta){

        if (verbose) cat(c(outer, lam), "\n", sep = "") # only for test purposes
        test_inds <- cumsum(folds_outer == outer)
        test_inds <- c(test_inds[folds_outer != outer], test_inds[n])

        tree <- BinarySegmentation(x[folds_outer != outer, , drop = F],
                                   x_test = x[folds_outer == outer, ], test_inds = test_inds,
                                   method = mth,
                                   NA_method = NA_mth, optimizer = opt, delta = del, lambda = lam,
                                   gamma = gamma, control = control)
        if (verbose) cat("fit finished. \n")
        print(tree)
      }}}
  tree
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
            gamma = x[["loss"]][, 1],
            loss = rowMeans(x[["loss"]][, -1]),
            n_cpts = 0
          )
        } else {
          data.frame(
            delta = factor(x[["delta"]]),
            lambda = formatC(x[["lambda"]], format = "e", digits = 2),
            gamma = x[["loss"]][, 1],
            loss = rowMeans(x[["loss"]][, -1]),
            n_cpts = rowMeans(apply(x[["cpts"]][, -1, drop = F], 2, function(x) sapply(x, length) - 2))
          )
        }
    )
  ) # substract first and last segment boundary

  res_long <- reshape2::melt(res_long, measure.vars = c("loss", "n_cpts"), value.name = "value", variable.name = "metric")

  metrics_names <- c(
    "loss" = "cv-loss",
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
    facet_grid(metric ~ ., scales = "free_y", switch = "y", labeller = labeller(metric = metrics_names)) +
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


SolutionPaths <- function(dat, var = "loss") {
  dat <- dat[order(dat$gamma), ]
  list(
    lambda = dat$lambda[1], delta = dat$delta[1],
    loss = ImputeMatrix(reshape2::dcast(dat, gamma ~ fold, value.var = var)),
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
  sum((y - x %*% beta - intercepts)^2)
}

GetOpt <- function(param_res) {
  avg_loss <- rowMeans(param_res$loss[, -1])
  opt <- which.min(avg_loss)

  data.frame(
    lambda = param_res$lambda,
    delta = param_res$delta,
    gamma = param_res$loss[opt, 1],
    loss = avg_loss[opt]
  )
}

LogSpace <- function(from, to, length.out) {
  exp(seq(log(from), log(to), length.out = length.out))
}

sample_folds <- function(n, k, randomize = TRUE){
  if (randomize){
    random_draw <- runif(n)
    k_quantiles <- quantile(random_draw, 0:k/k)
    cut(random_draw, k_quantiles, labels = 1:k, include.lowest = TRUE)
  } else {
    as.factor(rep(1:k, ceiling(n/k))[1:n])
  }
}
