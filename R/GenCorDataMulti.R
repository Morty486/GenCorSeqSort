#' Generate multivariate data with specified marginal distributions as well as a specified correlation structure
#'
#' This function generates multivariate mixed simulation data (n by k) given the specified marginal distributions and correlation structure.
#'
#' Generate n random samples from the intended k marginal distributions independently. To introduce a specified correlation structure for these random samples,
#' sort (k-1) sequences of numbers sequentially. In this function, the sequential sorting is achieved by keeping i th sequence of numbers unchanged and sorting
#' (i+1) th to k th sequences of numbers accordingly for i from 1 to (k-1).
#' For example, if we have 4 random variables X, Y, Z, and Q, and the input order is (X, Y, Z, Q), first we keep the sequence for X unchanged and sort Y, Z and Q accordingly.
#' Then we keep Y unchanged and sort Z, Q accordingly. Finally, we keep Z unchanged and sort Q accordingly. The input order will affect the range of feasible correlation structures
#' which can be achieved by this method. For higher dimensions (> 3), it is hard to calculate the conditional bounds for this method.
#' If k <= 3, we recommend to use \code{\link{GenCorData}} function.
#'
#'
#'
#' @param n  Sample size for the simulation data
#' @param lst A list of functions which generate data under specified marginal distributions separately
#' @param cor_mat  Specified correlation matrix
#'
#'
#' @return A list is returned which contains matrices of simulated data, generated correlations and specified correlations.
#'   \item{\code{sim_data}}{n by k matrix. Each column corresponds to a variable, and each row is one random sample.}
#'   \item{\code{gen_cor}}{Correlation matrix calculated from the simulated data}
#'   \item{\code{spec_cor}}{Specified correlation matrix}
#'
#' @references Demirtas (2019), Inducing Any Feasible Level of Correlation to
#'             Bivariate Data With Any Marginals, The American Statistician 73.3 (2019): 273-277
#'
#' @seealso \code{\link{Rank.Sort}}, \code{\link{Validate.Correlation}}, \code{\link{Compute.SortProp}}, \code{\link{Compute.PairBounds}}
#'
#'
#'
#' @examples
#' f1 = function(n){rnorm(n)}
#' cor_mat = matrix(
#' c(1,.1,.2,.3,.1,1,.4,.2,.2,.4,1,.2,.3,.2,.2,1),
#' nrow = 4)
#' GenCorDataK(10^5, list(f1,f1,f1,f1), cor_mat)
#'
#' @importFrom stats cor
#'
#' @export


GenCorDataMulti = function(n, lst, cor_mat){
  if (!nrow(cor_mat) == length(lst)) {
    stop("Dimension of correlation matrix does not match the number of variables! \n")
  }
  sim = sapply(lst, function(a){a(n)})
  n.var = length(lst)
  cor_mat_input = cor_mat
  pairbounds = Compute.PairBounds(lst)
  Validate.Correlation(cor_mat, pairbounds)

  prop_matrix = matrix(0, n.var, n.var)
  prop_matrix[upper.tri(prop_matrix)] = Compute.SortProp(cor_mat, pairbounds)

  low_bdd = pairbounds$low_bdd
  up_bdd = pairbounds$up_bdd

  n_start = 1
  for (i in 1:(n.var - 1)) {
    # Get the rank of ith sequence
    # Sort (n-i) sequences after the ith sequence according to the rank of ith sequence
    for (j in (i + 1):n.var) {
      n_end = floor(prop_matrix[i, j]*n)
      if ((n_start + n_end) > n){
        err = paste(c("After fixing the front correlations, ", "the correlation between ", as.character(i),
                      " and ", as.character(j), " is out of the bounds of this method"), collapse = "")
        stop(err)
      }
      sim[,j] = Rank.Sort(sim[,i], sim[,j], n_start:(n_start + n_end), cor_mat[i,j])
    }
    n_start = n_start + max(floor(n*prop_matrix[i, (i+1):n.var]))
    # Update correlation matrix
    if ((i+1) <= (n.var - 1)){
      for (x in (i+1):(n.var - 1)) {
        for (y in (x+1):n.var) {
          cor_mat[x, y] = cor_mat[x, y] -
            min(prop_matrix[i,x], prop_matrix[i,y])*ifelse((cor_mat[i,x]*cor_mat[i,y] > 0),
                                                           up_bdd[x,y], low_bdd[x,y])
        }
      }
    }
    # Update proportion we need to sort
    prop_matrix[upper.tri(prop_matrix)] = Compute.SortProp(cor_mat, pairbounds)
  }
  l = list(sim, round(cor(sim),4), round(cor_mat_input,4))
  names(l) = c("sim_data", "gen_cor", "spec_cor")
  l
}


