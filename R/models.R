#' ChainNetwork
#'
#' Spawn a chain network covariance matrix
#'
#' @param p Positive integer.The desired number of dimensions.
#' @param n_perm Positive integer. The first n_perm dimensions will be permuted randomly.
#' @param a Positive float between 0 and 1. Scale parameter for the elements of the covariance matrix
#' @param prec_mat Should the precision matrix be returned? If false the covariance matrix will be returned (default).
#' @param scaled If TRUE the created precision matrix will be inverted and scaled to a correlation matrix.
#'
#' @return A covariance or precision matrix.
#' @export
#'
#' @importFrom stats runif cov2cor
#'
#' @examples
#' ChainNetwork(50)
ChainNetwork <- function(p, n_perm = p, a = 0.5, prec_mat = F, scaled = T) {
  stopifnot(p >= n_perm)
  s_vec <- cumsum(runif(p, 0.5, 1))
  if (!is.null(n_perm) && n_perm >= 0) {
    perm_inds <- sample(1:n_perm, n_perm, replace = FALSE)

    if (n_perm < p) perm_inds <- c(perm_inds, (n_perm + 1):p)

    s_vec <- s_vec[perm_inds]
  }

  omega <- matrix(0, nrow = p, ncol = p)

  for (i in seq_len(p)) {
    for (j in seq_len(p)) {
      omega[i, j] <- exp(-a * abs(s_vec[i] - s_vec[j]))
    }
  }

  if (scaled) {
    sigma <- cov2cor(solve(omega))
    if (prec_mat) {
      solve(sigma)
    } else {
      sigma
    }
  } else {
    if (prec_mat) {
      omega
    } else {
      solve(omega)
    }
  }
}

#' ScaleNetwork
#'
#' Spawn a hub network covariance matrix
#'
#' @inheritParams ChainNetwork
#' @inheritParams RandomNetwork
#' @param preferential_power Power coefficient alpha for weighting of degree number as alpha in prefential attachment mechanism.
#'
#' @return A covariance or precision matrix.
#' @export
#'
#' @examples
#' ScaleNetwork(50)
ScaleNetwork <- function(p, preferential_power = 1, u = 0.1, v = 0.3, prec_mat = T, scaled = F) {
  theta <- matrix(0, p, p)
  probs <- numeric(p)

  theta[1, 2] <- theta[2, 1] <- TRUE

  for (i in seq(3, p)) {
    probs <- colSums(theta) ^ preferential_power
    probs <- probs / sum(probs)

    edges <- sample.int(i - 1, 1, prob = probs[1:(i - 1)])

    theta[edges, i] <- theta[i, edges] <- TRUE
  }

  diag(theta) <- 0

  omega <- theta * v

  diag(omega) <- abs(min(eigen(omega)$values)) + u

  if (scaled) {
    sigma <- cov2cor(solve(omega))
    if (prec_mat) {
      solve(sigma)
    } else {
      sigma
    }
  } else {
    if (prec_mat) {
      omega
    } else {
      solve(omega)
    }
  }
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
RandomNetwork <- function(p, prob = min(1, 5 / p), u = 0.1, v = 0.3, prec_mat = F, scaled = F) {
  theta <- matrix(0, p, p)

  tmp <- matrix(runif(p ^ 2, 0, 0.5), p, p)
  tmp <- tmp + t(tmp)
  theta[tmp < prob] <- 1
  diag(theta) <- 0

  omega <- theta * v
  diag(omega) <- abs(min(eigen(omega)$values)) + u

  if (scaled) {
    sigma <- cov2cor(solve(omega))
    if (prec_mat) {
      solve(sigma)
    } else {
      sigma
    }
  } else {
    if (prec_mat) {
      omega
    } else {
      solve(omega)
    }
  }
}

#' DiagMatrice
#'
#' Spawn a p-dimensional identity matrix.
#'
#' @inheritParams ChainNetwork
#'
#' @return A covariance matrix.
#' @export
#'
#' @examples
#' DiagMatrice(50)
DiagMatrice <- function(p) {
  diag(p)
}


#' MoveEdges
#'
#' Radomly move a share of the edges in a random graph
#'
#' In order to create a slightly different graph with the same level of sparsity
#' the selected share of edges will be randomly moved to positions where no edge existed
#' before. Make sure to choose share_moves low for dense graphs.
#'
#' @param x A precision matrix
#' @param share_moves Share of edges to be moved.
#' @param tol Tolerance for zero values.
#'
#' @export
MoveEdges <- function(x, share_moves = 0.1, tol = 1e-16) {
  if (share_moves == 0) return(x)

  edges <- which(upper.tri(x) & abs(x) >= tol)
  not_edges <- which(upper.tri(x) & x == 0)

  n_moves <- floor(share_moves * length(edges))

  if (length(not_edges) <= n_moves) {
    n_moves <- length(not_edges)
    warning("Cannot move edges because graph is not sparse enough. Will move as many as possible.")
  }

  d <- diag(x)

  # ensure that matrix is not pd to begin with
  diag(x) <- diag(x) - abs(min(eigen(x)$values)) - 0.1

  # sample until matrix is pd again
  while (any(eigen(x)$values <= 0)) {
    sel_edges <- sample(edges, n_moves)
    x[sample(not_edges, n_moves)] <- x[sel_edges]
    x[sel_edges] <- 0

    # make symmetric again
    x[lower.tri(x)] <- t(x)[lower.tri(x)]
    diag(x) <- d
  }
  x
}


#' RegrowNetwork
#'
#' Prune graph and regrow edges according of scale-free network
#'
#' Note that v and preferential_power agruments need to be equal to the ones which
#' initially created omega.
#'
#' @param omega A precision matrix as created by ScaleNetwork
#' @param n_nodes Number of nodes to prune and regrow. Default is 0.1 of all nodes.
#' @inheritParams ScaleNetwork
#' @inheritParams RandomNetwork
#'
#' @export
#'
#' @examples
#' omega <- ScaleNetwork(20, v = 1)
#' omega_hat <- RegrowNetwork(omega, v = 1)
#' omega
#' omega_hat
RegrowNetwork <- function(omega, n_nodes = ncol(omega) * 0.1, preferential_power = 1, v = 0.3) {
  d <- diag(omega)
  diag(omega) <- 0

  n_nodes <- floor(n_nodes)

  p <- ncol(omega)

  # prune
  pruned <- numeric(n_nodes)
  for (i in seq(1, n_nodes)) {
    candidates <- which(colSums(omega != 0) == 1)
    if (length(candidates) == 0) stop("Not enough nodes to prune from graph!")
    pruned[i] <- sample(candidates, 1)
    omega[, pruned[i]] <- omega[pruned[i], ] <- 0
  }

  # regrow
  pruned <- rev(pruned)
  probs <- numeric(p)

  omega_temp <- omega

  # Regrow network and discard permutation if it leads to non pd matrix
  while (any(eigen(omega)$values <= 0)) {
    omega <- omega_temp
    for (i in seq(1, n_nodes)) {
      probs <- colSums(omega != 0) ^ preferential_power
      if (sum(probs) == 0 | any(is.na(probs))) probs <- rep(1, p)
      probs <- probs / sum(probs)

      edges <- sample.int(p, 1, prob = probs)

      omega[edges, pruned[i]] <- omega[pruned[i], edges] <- v
    }
    diag(omega) <- d
  }
  omega
}
