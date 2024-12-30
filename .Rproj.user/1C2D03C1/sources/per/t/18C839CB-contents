#' Compute the proportion needed to be sorted between two variables
#'
#' This function calculates the proportion needed to be sorted for the targeted correlation.
#'
#'
#' @param cor_mat  Specified correlation matrix
#' @param pairbounds Calculated correlation bounds from the function \code{Compute.PairBounds}
#'    based on Demirtas and Hedeker (2011).
#'
#' @return It will return the proportions needed to be sorted.
#'
#' @seealso \code{\link{Compute.PairBounds}}
#'
#' @examples
#'  # Calculate correlation bounds of each pairs using the function Compute.PairBounds.
#'  # Three variables of normal,rnorm(n), Poisson,rpois(n,1), Binomial,
#'  # rbinom(n,10,0.5) are generated.
#'  # Define the 3x3 correlation matrix, cor_mat having the correlation coefficient
#'  # 0.5 for all three pairs.
#'
#'  f1 = function(n){rnorm(n)}
#'  f2 = function(n){rpois(n,1)}
#'  f3 = function(n){rbinom(n,10,0.5)}
#'  pairbound = Compute.PairBounds(list(f1,f2,f3))
#'  cor_mat = matrix(0.5, 3, 3) + diag(0.5,3)
#'  Compute.SortProp(cor_mat, pairbound)
#'
#'
#'
#' @export

Compute.SortProp = function(cor_mat, pairbounds) {
  low_bdd = pairbounds$low_bdd
  up_bdd = pairbounds$up_bdd
  dim = dim(cor_mat)[1]
  prop = matrix(0, dim, dim)
  for (i in 1:(dim - 1)) {
    for (j in (i+1): dim) {
      if (cor_mat[i, j] <= 0) {prop[i,j] = cor_mat[i, j]/low_bdd[i, j]}
      else{prop[i,j] = cor_mat[i, j]/up_bdd[i, j]}
    }
  }
  prop[upper.tri(prop)]
}


