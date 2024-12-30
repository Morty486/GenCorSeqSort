#' Validate the correlation matrix
#'
#' This function validates the correlation matrix whether meets the conditions for this method (1) positive semi-definite;
#'  (2) correlation coefficient is not 0; (3) is within correlation bounds
#'  calculated from the package \code{Compute.PairBounds} by Demirtas and Hedeker(2011).
#'
#'
#' @param cor_mat  Specified correlation matrix
#' @param pairbounds Computed correlation bounds from the package \code{Compute.PairBounds}.
#'
#' @return An error message will be displayed when the matrix doesn't satisfy any of the three
#'  conditions.
#'
#' @seealso \code{\link{Compute.PairBounds}}
#'
#' @examples
#'
#'  f1 = function(n){rnorm(n)}
#'  f2 = function(n){rpois(n,1)}
#'  f3 = function(n){rbinom(n,10,0.5)}
#'  pairbound = Compute.PairBounds(list(f1,f2,f3))
#'  cor_mat = matrix(0.5, 3, 3) + diag(0.5,3)
#'  Compute.SortProp(cor_mat, pairbound)
#'  Validate.Correlation(cor_mat = cor_mat,pairbound)
#'
#' @importFrom matrixcalc is.positive.semi.definite
#'
#'
#'
#' @export


Validate.Correlation = function(cor_mat, pairbounds){
  if (!matrixcalc::is.positive.semi.definite(cor_mat)) {
    stop("The correlation matrix is not positive semi-definite!\n")
  }

  if (all(cor_mat[upper.tri(cor_mat)] == 0)){
    stop("No correaltions between any variables!\n")
  }

  low_bdd = pairbounds$low_bdd
  up_bdd = pairbounds$up_bdd

  for (i in 1:(nrow(cor_mat) - 1)) {
    for (j in (i + 1):nrow(cor_mat)) {
      if (cor_mat[i, j] > up_bdd[i, j] || cor_mat[i, j] < low_bdd[i, j]) {
        err = paste(c("correlation between variables ", as.character(i),
                      " and ", as.character(j)," should be between ",
                      as.character(round(low_bdd[i, j], 4)), " and ", as.character(round(up_bdd[i, j], 4)),
                      " (Demirtas-Hedeker Bounds)"), collapse = "")
        stop(err)
      }
    }
  }
}

