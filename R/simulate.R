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

  for(i in seq_along(seg_lengths)){
    seg_start <- ifelse(i == 1, 1 ,sum(seg_lengths[(i-1):1]) + 1)
    seg_end   <- seg_start + seg_lengths[i] - 1
    data[seg_start:seg_end, ] <- MASS::mvrnorm(seg_lengths[[i]], model[["segment_means"]][[i]], solve(model[["icov_mats"]][[i]]))
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
#' @param mean_zero If TRUE, the mean for each segment will be zero.
#' Otherwise it will be sampled from a standard normal distribution for each segment.
#' @param modelFUN A function that spawns covariance matrices of dimension p.
#' @param ... Addtional arguments to be supplied to modelFUN
#'
#' @return An object to be used by \link{SimulateFromModel}
#' @export
CreateModel <- function(n_segments, n, p, modelFUN, equispaced = T, mean_zero = T, ...){

  model_args <- list(...)

  if (equispaced){
    segment_lengths   <- c(rep(ceiling(n/n_segments), times = n_segments - 1), n - (n_segments - 1)* ceiling(n/n_segments))
    changepoints      <- (cumsum(segment_lengths) + 1)[- length(segment_lengths)]
  } else {
    changepoints      <- sort(sample(2:(n-1), size = n_segments - 1, replace = F))
    segment_lengths   <- c(changepoints - c(0, changepoints[-length(changepoints)]), n - changepoints[length(changepoints)])
  }
  segment_means      <- replicate(n_segments, rep(if(mean_zero) 0 else rnorm(1), p), simplify = F)
  icov_mats          <- replicate(n_segments, do.call(modelFUN, c(list(p = p), model_args)), simplify = F)


  list(segment_lengths = segment_lengths,
       segment_means = segment_means,
       icov_mats = icov_mats,
       true_changepoints = changepoints)
}

