#' ChainNetwork
#'
#' Spawn a chain network covariance matrix
#'
#' @param p Positive integer.The desired number of dimensions.
#' @param n_perm Positive integer. The first n_perm dimensions will be permuted randomly.
#' @param a Positive float between 0 and 1. Scale parameter for the elements of the covariance matrix
#'
#' @return A covariance matrix.
#' @export
#'
#' @examples
#' ChainNetwork(50)
ChainNetwork <- function(p, n_perm = p, a = 0.5) {
  stopifnot(p >= n_perm)
  s_vec <- cumsum(runif(p, 0.5, 1))
  if (!is.null(n_perm) && n_perm >= 0){

    perm_inds <- sample(1:n_perm, n_perm, replace = FALSE)

    if(n_perm < p) perm_inds <- c(perm_inds, (n_perm + 1):p)

    s_vec <- s_vec[perm_inds]
  }

  icov_mat <- matrix(0, nrow = p, ncol = p)

  for(i in seq_len(p)) {
    for(j in seq_len(p)) {
      icov_mat[i, j] <- exp(-a*abs(s_vec[i]- s_vec[j]))
    }
  }
  icov_mat
}

#' HubNetwork
#'
#' Spawn a hub network covariance matrix
#'
#' @inheritParams ChainNetwork
#' @param n_hubs Number of hubs in the network
#'
#' @return A covariance matrix.
#' @export
#'
#' @examples
#' HubNetwork(50)
HubNetwork <- function(p, max_hubs = p/10) {
  n_hubs <- sample(seq_len(max_hubs), 1)
  huge::huge.generator(n = 2, d = p, graph = "hub", g = n_hubs, verbose = FALSE)$sigma
}

#' RandomNetwork
#'
#' Spawn a erdos renyi type random network covariance matrix
#'
#' @inheritParams ChainNetwork
#' @param prob Probabilty that a pair of nodes have a common edge.
#'
#' @return A covariance matrix.
#' @export
#'
#' @examples
#' RandomNetwork(50)
RandomNetwork <- function(p, prob = 3/p) {
  huge::huge.generator(n = 2, d = p, graph = "random", prob = prob, verbose = FALSE)$sigma
}
