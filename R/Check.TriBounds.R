#' Check the bounds for the third pair's correlation in the trivariate setting
#'
#' This function checks the extra bounds in the trivariate extension for this method.
#'
#'
#' @param cor  A vector that contains the sorted correlations
#' @param pv  A vector contains the sorted proportions
#' @param low_bdd A vector contains the sorted lower bounds
#' @param up_bdd A vector contains the sorted upper bounds
#' @param ord A vector contains the order for each variable
#' @param output An option for whether gives an output
#'
#' @return An error message will be displayed when the third correlation is not within the trivariate bounds.
#'
#' @examples
#'  f1 = function(n){rnorm(n)}
#'  cor_mat = matrix(c(1,.49,.1, .49, 1, -.4, .1, -.4, 1), nrow = 3)
#'  l = Find.Order(list(f1, f1, f1), cor_mat)
#'  cor = l[[7]]
#'  pv = l[[6]]
#'  low_bdd = l[[3]]
#'  up_bdd = l[[4]]
#'  ord = l[[1]]
#'  Check.TriBounds(cor, pv, low_bdd, up_bdd, ord, TRUE)
#'
#' @export

Check.TriBounds = function(cor, pv, low_bdd, up_bdd, ord = c(1,2,3), output = T){
  a = min(pv[1], pv[2])
  b = max(pv[1], pv[2])
  c = low_bdd[3]
  d = up_bdd[3]

  if (cor[1]*cor[2] > 0) {
    low = min(c( ((1 - b)*c + a*d), (1 - a - b)*c))
    up = (1 - b)*d + a*d
  } else{
    low = (1 - b)*c + a*c
    up = max(c( ((1 - b)*d + a*c), (1 - a - b)*d))
  }
  bdd = c(low, up)

  err = paste(c("Correlation between last two variables should be between ",
                as.character(round(bdd[1], 4)), " and ",as.character(round(bdd[2], 4)),
                " in the sequence ", paste(as.character(ord), collapse = ",")), collapse = "")
  if (cor[3] > bdd[2] || cor[3] < bdd[1]) {
    stop(err)
  }
  if (output == T) {
    list(bdd,err)
  }
}


