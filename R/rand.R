#' generate all uniform random values matrix/tensor
#' @param n,d,g tensor with dimension [n,d,g]
#' @noRd

rand <- function(n, d, g = 1) {
  if (g == 1) {
    return(matrix(stats::runif(n * d), n, d))
  }
  else{
    return(array(stats::runif(n * d), dim = c(n, d, g)))
  }
}
