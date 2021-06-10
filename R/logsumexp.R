#' log of sum of exponential of A by either row or column
#' @param A input matrix
#' @param margin 1 indicates rows and 2 for columns
#' @noRd

logsumexp <- function(A, margin) {
  M <- A

  if (margin == 1) {
    xstar <- apply(A, 1, max)
    M <- M - xstar %*% matrix(1, nrow = 1, ncol = ncol(M))
    A <- xstar + log(apply(exp(M), 1, sum))
  } else{
    xstar <- apply(A, 2, max)
    M <- M -  matrix(1, nrow = nrow(M), ncol = 1) %*% xstar
    A <- xstar + log(apply(exp(M), 2, sum))
  }
  return(A)
}
