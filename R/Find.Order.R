#' Find the order for sorting
#'
#' This function finds the order for trivariate sorting
#'
#'
#' @param lst A list of functions which generate data under specified marginal distributions separately
#' @param cor_mat  Specified correlation matrix
#'
#'
#' @return A list is returned which contains the sorting order, ordered variables,  ordered lower bounds, ordered upper bounds, names of ordered correlation, ordered proportion vectors, ordered proportion vectors and ordered correlations.
#'   \item{\code{order}}{The order for sorting which achieves the maximum bounds}
#'   \item{\code{ordered list}}{An ordered list which generates the random variables}
#'   \item{\code{ordered lower bounds}}{The lower bounds corresponding to the ordered variables}
#'   \item{\code{ordered upper bounds}}{The upper bounds corresponding to the ordered variables}
#'   \item{\code{names of ordered correlation}}{Names of ordered correlation}
#'   \item{\code{ordered proportion vectors}}{The sorting proportion corresponding to the ordered variables}
#'   \item{\code{ordered correlations}}{Ordered correlations}
#'
#'
#' @seealso \code{\link{Compute.SortProp}}, \code{\link{Compute.PairBounds}}, \code{\link{CorOrder.ToVar}}
#'
#'
#'
#' @examples
#'  f1 = function(n){rnorm(n)}
#'  cor_mat = matrix(c(1,.49,.1, .49, 1, -.4, .1, -.4, 1), nrow = 3)
#'  Find.Order(list(f1, f1, f1), cor_mat)
#'
#'
#'
#' @export


Find.Order = function(lst, cor_mat){

  pairbounds = Compute.PairBounds(lst)
  low_bdd = pairbounds$low_bdd
  up_bdd = pairbounds$up_bdd
  prop_vec = Compute.SortProp(cor_mat, pairbounds)

  name = NULL
  upbdd_vec = NULL
  lowbdd_vec = NULL

  for (i in 1:(nrow(cor_mat) - 1)) {
    for (j in (i + 1):nrow(cor_mat)) {
      name = rbind(name, c(i, j))
      lowbdd_vec = c(lowbdd_vec, low_bdd[i,j])
      upbdd_vec = c(upbdd_vec, up_bdd[i,j])
    }
  }
  r = order(prop_vec)
  srt_pv = prop_vec[r]
  srt_name = name[r, ]
  ord = CorOrder.ToVar(srt_name)
  srt_lst = lst[ord]
  srt_low_bdd = lowbdd_vec[r]
  srt_up_bdd = upbdd_vec[r]
  srt_cor = sapply(1:length(ord), function(i){cor_mat[srt_name[i, 1], srt_name[i, 2]]})
  srt_lst = lst[ord]

  l = list(ord, srt_lst, srt_low_bdd, srt_up_bdd, srt_name, srt_pv, srt_cor)
  names(l) = c("order", "ordered list", "ordered lower bounds", "ordered upper bounds","names of ordered correlation", "ordered proportion vectors",
               "ordered correlations")
  l
}
