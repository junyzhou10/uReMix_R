#' generate all zero matrix/tensor
#' @param n,d,g tensor with dimension [n,d,g]
#' @noRd

zeros <- function(n, d, g = 1) {
  if (g == 1) {
    return(matrix(0, n, d))
  }
  else{
    return(array(0, dim = c(n, d, g)))
  }
}
