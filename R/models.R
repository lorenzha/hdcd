#' ChainNetwork
#'
#' Spawn a chain network covariance matrix
#'
#' @param p Positive integer.The desired number of dimensions.
#' @param n_perm Positive integer. The first n_perm dimensions will be permuted randomly.
#' @param a Positive float between 0 and 1. Scale parameter for the elements of the covariance matrix
#' @param cov_mat Should the precision matrix be returned? If false the covariance matrix will be returned (default).
#'
#' @return A covariance matrix.
#' @export
#'
#' @examples
#' ChainNetwork(50)
ChainNetwork <- function(p, n_perm = p, a = 0.5, prec_mat = F) {
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
  if (prec_mat)
    icov_mat
  else
    solve(icov_mat)
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
#' @param u Constant added to the diagonal elements of the precision matrix for controlling the magnitude of partial correlations.
#' @param v Constant added to the off diagonal of the precision matrix for controlling the magnitude of partial correlations.
#'
#' @return A covariance matrix.
#' @export
#'
#' @examples
#' RandomNetwork(50)
RandomNetwork <- function(p, prob = min(1, 3/p), prec_mat = F, u = 0.1, v = 0.3) {

  theta <- matrix(0, p, p)

  tmp <- matrix(runif(p^2, 0, 0.5), p, p)
  tmp <- tmp + t(tmp)
  theta[tmp < prob] = 1

  omega <- theta * v
  diag(omega) <- abs(min(eigen(omega)$values)) + 0.1 + u

  if (prec_mat)
    omega
  else
    solve(omega)
}


#' MoveEdges
#'
#' Radomly move a share of the edges in a random graph
#'
#' In order to create a slightly different graph with the same level of sparsity
#' the selected share of edges will be randomly moved to positions where no edge existed
#' before. Make sure to choose share_moves low for dense graphs.
#'
#' @param mat A precision matrix
#' @param share_moves Share of edges to be moved.
#'
#' @return
#' @export
#'
#' @examples
MoveEdges <- function(prec_mat, share_moves = 0.1){

  if (share_moves == 0) return(prec_mat)

  edges <- which(upper.tri(prec_mat) & prec_mat != 0)
  not_edges <- which(upper.tri(prec_mat) & prec_mat == 0)

  n_moves <- floor(share_moves * length(edges))

  if(length(not_edges) <= length(n_moves)){
    n_moves <- length(not_edges)
    warning("Cannot move edges because graph is not sparse enough. Will move as many as possible.")
  }

  d <- diag(prec_mat)
  diag(prec_mat) <- 0

  while (any(eigen(prec_mat)$values <= 0)){

    sel_edges <- sample(edges, n_moves)
    prec_mat[sample(not_edges, n_moves)] <- prec_mat[sel_edges]
    prec_mat[sel_edges] <- 0

    prec_mat[lower.tri(prec_mat)] <- t(prec_mat)[lower.tri(prec_mat)]
    diag(prec_mat) <- d
  }
  prec_mat
}
