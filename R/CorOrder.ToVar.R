#' Correlations order to the sorting order of variables
#'
#' This function returns an order of variables according the order of correlations.
#'
#' @param cor_name A list of correlation names following a certain order
#'
#' @return It will return an order for variables.
#'
#'
#'
#' @export

CorOrder.ToVar = function(cor_name) {
  indices = unique(as.vector(cor_name))
  ans = intersect(cor_name[1, ], cor_name[2, ])
  ans = c(ans, setdiff(cor_name[1, ], ans))
  ans = c(ans, setdiff(indices, ans))
  ans
}
