#' Compute lower and upper bounds of correlations for each pair of variables
#'
#' This function computes the lower and upper bounds of pairwise correlations via a sorting approach.
#'  Sorting in the same order either in ascending or descending makes two variables have the positive upper bound correlation,
#'  whereas the opposite sorting order produces the negative lower bound coefficients.
#'  This widely-applicable and practical sorting method was proposed by Demirtas and Hedeker (2011).
#'
#' @param lst A list of functions which generate data under specified marginal distributions
#'
#' @return A list is returned which contains the lower and upper bound of correlation matrices.
#'   \item{\code{low_bdd}}{A matrix showing the lower bound of correlation coefficients of each pair of the functions listed in (\code{lst}).}
#'   \item{\code{up_bdd}}{A matrix showing the upper bound of correlation coefficients of each pair of the functions listed in (\code{lst}).}
#'
#' @references Demirtas and Hedeker (2011), A Practical Way for Computing Approximate Lower and Upper Correlation Bounds,
#'             The American Statistician,65(2):104-109
#'
#' @examples
#'  # Three variables f1,f2,f3 were entered for the trivariate correlations
#'  # having the standard normal distribution, Poisson with lambda=1, and
#'  # the binomial size=10 with p=0.5 using randomly generated n=10^6 samples.
#'  # Then, the upper and lower bound correlations were obtained between three pairs, f1,f2,f3.
#'  f1 = function(n){rnorm(n)}
#'  f2 = function(n){rpois(n,1)}
#'  f3 = function(n){rbinom(n,10,0.5)}
#'  Compute.PairBounds(list(f1,f2,f3))
#'
#' @export

Compute.PairBounds = function(lst){
  n = 10^6
  dim = length(lst)
  low_bdd = matrix(0, dim, dim)
  up_bdd = matrix(0, dim, dim)
  for (i in 1:(dim - 1)) {
    for (j in (i+1):dim) {
      x = lst[[i]](n)
      y = lst[[j]](n)
      low_bdd[i, j] = stats::cor(x[order(x)], rev(y[order(y)]))
      up_bdd[i, j] = stats::cor(x[order(x)], y[order(y)])
    }
  }
  low_bdd = low_bdd + t(low_bdd) + diag(1, dim)
  up_bdd = up_bdd + t(up_bdd) + diag(1, dim)
  list(low_bdd = low_bdd, up_bdd = up_bdd)
}





