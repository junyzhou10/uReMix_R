#' generate all one matrix/tensor
#' @param n,d,g tensor with dimension [n,d,g]
#' @noRd

ones <- function(n, d, g = 1) {
  if (g == 1) {
    return(matrix(1, n, d))
  }
  else{
    return(array(1, dim = c(n, d, g)))
  }
}
