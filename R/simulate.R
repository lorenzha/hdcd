#' SimulateFromModel
#'
#' @param model An object as created by \link{CreateModel}
#'
#' @return A n times p matrix of simulated data.
#' @export
#'
#' @examples
#'
#' # Simulate 100 observations from a 10-dimensional chain network with two changepoints
#' mod <- CreateModel(3, 100, 10, ChainNetwork)
#' SimulateFromModel(mod)
SimulateFromModel <- function(model) {
  seg_lengths <- model$segment_lengths

  data <- matrix(NA, nrow = sum(seg_lengths), ncol = length(model[["segment_means"]][[1]]))

  for (i in seq_along(seg_lengths)) {
    seg_start <- ifelse(i == 1, 1, sum(seg_lengths[(i - 1):1]) + 1)
    seg_end <- seg_start + seg_lengths[i] - 1
    data[seg_start:seg_end, ] <- MASS::mvrnorm(seg_lengths[[i]], model[["segment_means"]][[i]], model[["cov_mats"]][[i]])
  }
  return(data)
}

#' CreateModel
#'
#' Create a model to generate data from for simulating the detection of changepoints
#'
#' @param n_segments Number of segments in the model. There will be one less changepoints than segments.
#' @param n Number of observations
#' @param p Number of dimensions
#' @param equispaced If TRUE, the segments will be of equal length (up to rounding) and hence the changepoints will be equispaced.
#' Otherwise the changepoints will be drawn randomly and the distance between them will differ.
#' @param mean_vecs If NULL, the mean for each segment will be zero.
#' Otherwise mean_vecs should be a list containing a p-dimensional numeric vector with the means for each segment.
#' @param modelFUN A function that spawns covariance matrices of dimension p.
#' @param ... Addtional arguments to be supplied to modelFUN
#'
#' @return An object to be used by \link{SimulateFromModel}
#' @export
CreateModel <- function(n_segments, n, p, modelFUN, changepoints = NULL, equispaced = T, mean_vecs = NULL, ...) {
  model_args <- list(...)

  if (equispaced) {
    segment_lengths <- c(rep(ceiling(n / n_segments), times = n_segments - 1), n - (n_segments - 1) * ceiling(n / n_segments))
    changepoints <- (cumsum(segment_lengths) + 1)[-length(segment_lengths)]
  } else {
    if(is.null(changepoints)){
      changepoints <- sort(sample(2:(n - 1), size = n_segments - 1, replace = F))
    }
    if(length(changepoints) == 0){
      segment_lengths <- n
    } else {
      segment_lengths <- c(changepoints - c(0, changepoints[-length(changepoints)]), n - changepoints[length(changepoints)])
    }
  }

  if (is.null(mean_vecs)) {
    segment_means <- replicate(n_segments, rep(0, p), simplify = F)
  } else {
    segment_means <- mean_vecs
  }

  cov_mats <- replicate(n_segments, do.call(modelFUN, c(list(p = p), model_args)), simplify = F)

  list(
    segment_lengths = segment_lengths,
    segment_means = segment_means,
    cov_mats = cov_mats,
    true_changepoints = changepoints
  )
}


delete_values <- function(x, m, missingness = 'mcar', x_comp = F){
  n <- nrow(x)
  p <- ncol(x)
  if (missingness == 'mcar'){
    inds <- sample(which(!is.na(x)), floor(m * n * p), replace = F)
  } else if (missingness == 'mar'){
    i <- min(ceiling(m * p * 2), p)
    inds <- sample(1 : (n * i), floor(m * n * p), replace = F)
  } else if (missingness == 'nmar'){
    inds <- sample(1 : (n * p), floor(m * n * p), replace = F, prob = abs(c(x)))
  } else if (missingness == 'blockwise'){
    total <- prod(dim(x))
    missing <- sum(is.na(x))
    m <- m + missing / total
    missing_max <- ceiling(m * total)
    while (missing < missing_max){
      l <- min(floor(rexp(1, 8/n)), missing_max - missing)
      j <- sample(p, 1)
      k <- sample((1 - floor(l/2)) : (n + floor(l/2) - 1), 1)
      int <- max((k - floor(l/2)), 1) : min((k + floor(l/2)), n)
      x[int, j] <- NA
      missing <- sum(is.na(x))
    }
    return(x)
  } else if (missingness == 'both'){
    x <- delete_values(x, m/2, 'blockwise')
    return(delete_values(x, m/2))
  } else if (missingness == 'groundwater'){
    total <- prod(dim(x))
    missing <- sum(is.na(x))
    m <- m + missing / total
    missing_max <- ceiling(m * total)
    while (missing < missing_max){
      o <- rpois(1, p / 20)
      l <- min(floor(rexp(1, 8/n)), ceiling((missing_max - missing)/o))
      j <- sample(p, o)###lala
      k <- sample((1 - floor(l/2)) : (n + floor(l/2) - 1), 1)
      int <- max((k - floor(l/2)), 1) : min((k + floor(l/2)), n)
      x[int, j] <- NA
      missing <- sum(is.na(x))
    }
    return(x)
  } else if (missingness == 'none'){
    inds <- numeric(0)
  }

  x_del <- x
  x_del[inds] <- NA

  if (x_comp) {
    x_comp <- x
    x_comp[sample(1 : n, ceiling(m * n), replace = F), ] <- NA
    return(list(x_del = x_del, x_comp = x_comp))
  } else {
    return(x_del)
  }
}


plot_missingness_structure <- function(x){
  dt <- data.table::data.table(is.na(x))
  colnames(dt) <-  as.character(1 : ncol(x))
  dt$index <- 1 : nrow(x)
  dt <- data.table::melt(dt, id.vars = 'index')
  dt$value <-  ifelse(dt$value, 'missing', 'not_missing')
  ggplot2::ggplot(dt, ggplot2::aes(x = index, y = variable, col = value)) + ggplot2::geom_point()
}
