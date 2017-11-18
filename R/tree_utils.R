#' PruneTreeGamma
#'
#' Prune a bs_tree object for different values of tuning parameter gamma
#'
#' This functions adds the gamma penalty to each split in the tree. If the penalized loss is higher
#' than the loss of the base segment the node will be pruned from the tree. The tree will be pruned for the specified
#' path of gamma and for each gamma value a set of changepoints is returned.
#'
#'
#' @param tree An object of class \strong{bs_tree}
#' @param gamma_max Upper limit of gamma. Range will be [0, gamma_max].
#' @param gamma_length Number of equispaced points in the range.
#'
#' @return
#' \describe{
#'   \item{cpts}{A list of lists with changepoints corresponding to a value of gamma each.}
#'   \item{gamma}{Gamma values.}
#' }
#'
#' @export
PruneTreeGamma <- function(tree, gamma_max = 3, gamma_length = 50) {
  stopifnot(is(tree, "bs_tree"))

  cpts <- list()
  if (gamma_length > 1)
    gamma_seq <- seq(0, gamma_max, length.out = gamma_length)
  else
    gamma_seq <- gamma_max

  for (i in seq_along(gamma_seq)) {
    FUN <- PenalizeSplitsFUN(gamma_seq[i])

    clone_tree <- data.tree::Clone(tree, pruneFun = FUN)

    cpts[[i]] <- GetChangePointsFromLeafs(clone_tree)
  }
  list(cpts = cpts, gamma = gamma_seq)
}

#' GetChangePointsFromLeafs
#'
#' Utility function to get the changepoints from the name of the leaf nodes in the tree.
#'
#' @param tree An object of class \strong{bs_tree}
#'
#' @return A vector with the sorted changepoints.
GetChangePointsFromLeafs <- function(tree) {
  stopifnot(is(tree, "bs_tree"))

  filter <- function(x)
    if (x$name != "1" && x$name != "bs_tree") x$name
  else
    NA

  unname(sort(as.numeric(tree$Get(filter, filterFun = data.tree::isLeaf))))
}

#' PenalizeSplitsFUN
#'
#' Utility function to determine whether a node should be pruned from tree for
#' given gamma value. This is a closure and hence returns a function parametrized with gamma
#' to be used as prune function.
#'
#' @param gamma Split penality
#'
#' @return FALSE if the node should be pruned from tree
PenalizeSplitsFUN <- function(gamma) {
  function(node) {
    node <- data.tree::Navigate(node, "..") # We want to prune the parent tree where the split would've occured

    if (is.null(node$min_loss) || is.null(node$segment_loss))
      TRUE
    else
      node$min_loss + gamma < node$segment_loss
  }
}

#' print.bs_tree
#'
#' S3 method for printing a bs_tree object.
#'
#' Decorate the print method of the data.tree package to see more details at each node.
#'
#' @param tree A data.tree node
#' @export
print.bs_tree <- function(tree, ...) {
  NextMethod(generic = NULL, object = NULL, "start", "end", "min_loss", "segment_loss", ...)
  # data.tree::print.Node(node, "start", "end", "min_loss", "segment_loss")
}
