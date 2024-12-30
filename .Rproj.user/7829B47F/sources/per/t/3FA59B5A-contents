#' Sort the sequence of a sample according to the order of another sample
#'
#' This function sorts a sample data depending on the other reference data with the same ascending order for the positive correlation,
#'  and reversed opposite order for the negative correlation as much as targeted proportion.
#'
#' @param seq1 The first sequence which will be used as reference to sort the second sequence.
#' @param seq2 The second sequence that will be sorted following the seq1.
#' @param sort_part The proportion of data that needs to be sorted.
#' @param cor Correlation coefficient. Here it only matters whether it is positive or negative number. Any number between
#'    -1 to 1 can be entered. A positive number sort with the same order with the first sample sequence, and the negative number
#'    with reversed order.
#'
#' @return It will return a sorted sequence of the second sample, seq2.
#'
#' @examples
#'  # The first sample, seq1 has value (2,2,6,1,3,5) with three tied values,
#'  # and the second sample has (5,4,2,7,3,6).
#'  # For the positive correlation, cor=1 was entered to sort the data in an ascending order.
#'  # For the sample with six values \code{sort_part} (1:6) will sort 100 percent,
#'  # and (1:3) will sort 50 percent.
#'
#'   seq1 <- c(2,2,6,1,3,5)
#'   seq2 <- c(5,4,2,7,3,6)
#'   Rank.Sort(seq1, seq2, 1:6, 1)
#'
#'
#'
#' @export


Rank.Sort = function(seq1, seq2, sort_part, cor) {
  if (cor > 0) {
    seq2[sort_part] = sort(seq2[sort_part])[rank(seq1[sort_part], ties.method = "first")]
  }else{
    seq2[sort_part] = rev(sort(seq2[sort_part]))[rank(seq1[sort_part], ties.method = "first")]
  }
  return(seq2)
}







