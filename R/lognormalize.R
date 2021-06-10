#' Logarithmic normalization of data
#' @param M input matrix
#' @noRd

lognormalize <- function(M) {
  if (!is.matrix(M)) {
    M <- matrix(M)
  }
  n <- nrow(M)
  d <- ncol(M)
  a <- apply(M, 1, max)
  return(M - repmat(a + log(rowSums(exp(M - repmat(a, 1, d)))), 1, d))
}
