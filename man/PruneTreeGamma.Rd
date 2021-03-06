% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree_utils.R
\name{PruneTreeGamma}
\alias{PruneTreeGamma}
\title{PruneTreeGamma}
\usage{
PruneTreeGamma(x, gamma = seq(0, 3, length.out = 50))
}
\arguments{
\item{x}{An object of class \strong{bs_tree}}

\item{gamma}{A numeric vector with values for gamma for which the tree shall be pruned.}
}
\value{
\describe{
  \item{cpts}{A list of lists with changepoints corresponding to a value of gamma each.}
  \item{gamma}{Gamma values.}
}
}
\description{
Prune a bs_tree object for different values of tuning parameter gamma
}
\details{
This functions adds the gamma penalty to each split in the tree. If the penalized loss is higher
than the loss of the base segment the node will be pruned from the tree. The tree will be pruned for the specified
path of gamma and for each gamma value a set of changepoints is returned.
}
