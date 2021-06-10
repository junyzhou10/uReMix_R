#' layout matrix repeatedly
#' @param M matrix to repeat
#' @param n,d dimension of repeat [n,d]
#' @noRd

repmat <- function(M, n, d) {
  return(kronecker(matrix(1, n, d), M))
}

