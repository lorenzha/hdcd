#' CompareChangepointsRand
#'
#' Rand type performance indices
#'
#' Calculate rand type performance indices for two sets of changepoints. Typically one
#' of them will be the oracle estimate. See clues package for more details.
#'
#' @param cpts_a A sequence of changepoints.
#' @param cpts_b A sequence of changepoints.
#' @param n Total size of dataset from which both changepoint estimates originate.
#'
#' @return Returns a vector of the index values.
#' @export
#'
#' @examples
#' CompareChangepointsRand(c(20, 50), c(30, 70), 100)
CompareChangepointsRand <- function(cpts_a, cpts_b, n) {

  cpts_a <- sort(cpts_a)[!duplicated(sort(cpts_a))]
  cpts_b <- sort(cpts_b)[!duplicated(sort(cpts_b))]

  MarkGroupings <- function(cpts){
    diffs <- c(cpts, n) - c(0, cpts)
    rep(1:length(diffs), diffs)
  }

  clues::adjustedRand(MarkGroupings(cpts_a), MarkGroupings(cpts_b))
}

