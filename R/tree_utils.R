#' PruneTreeGamma
#'
#' Prune a bs_tree object for different values of tuning parameter gamma
#'
#' This functions adds the gamma penalty to each split in the tree. If the penalized loss is higher
#' than the loss of the base segment the node will be pruned from the tree. The tree will be pruned for the specified
#' path of gamma and for each gamma value a set of changepoints is returned.
#'
#' @param x An object of class \strong{bs_tree}
#' @param gamma A numeric vector with values for gamma for which the tree shall be pruned.
#'
#' @return
#' \describe{
#'   \item{cpts}{A list of lists with changepoints corresponding to a value of gamma each.}
#'   \item{gamma}{Gamma values.}
#' }
#'
#'
#' @importFrom methods is
#'
#' @export
PruneTreeGamma <- function(x, gamma = seq(0, 3, length.out = 50)) {
  stopifnot(is(x, "bs_tree"))

  cpts <- list()
  pruned_tree <- list()
  for (gam in gamma) {
    FUN <- PenalizeSplitsFUN(gam)

    clone_tree <- data.tree::Clone(x, pruneFun = FUN) # TODO: Check if copy can be avoided

    cpts[[as.character(gam)]] <- GetChangePointsFromLeafs(clone_tree)
    pruned_tree[[as.character(gam)]] <- clone_tree
  }
  list(cpts = cpts, gamma = gamma, pruned_tree = pruned_tree)
}

#' GetChangePointsFromLeafs
#'
#' Utility function to get the changepoints from the name of the leaf nodes in the tree.
#'
#' @importFrom methods is
#'
#' @param x An object of class \strong{bs_tree}
#'
#' @export
#' @return A vector with the sorted changepoints.
GetChangePointsFromLeafs <- function(x) {
  stopifnot(is(x, "bs_tree"))

  filter <- function(x)
    if (x$name != "1" && x$name != "bs_tree") {
      x$name
    } else {
      NA
    }

  unname(sort(as.numeric(x$Get(filter, filterFun = data.tree::isLeaf))))
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

    if (is.null(node$max_gain)) {
      TRUE
    } else {
      node$max_gain - gamma > 0
    }
  }
}

#' print.bs_tree
#'
#' S3 method for printing a bs_tree object.
#'
#' Decorate the print method of the data.tree package to see more details at each node.
#'
#' @param x A data.tree node.
#' @param ... Further arguments passed to print generic.
#'
#' @export
print.bs_tree <- function(x, ...) {
  NextMethod(generic = NULL, object = NULL, "start", "end", "max_gain", ...)
}
