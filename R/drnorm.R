#' generate all normal random values matrix
#' @param n,d matrix with dimension [n,d]
#' @param mean,sd mean and sd of normal distribution
#' @noRd

drnorm <- function(n, d, mean, sd) {
  A <- matrix(nrow = n, ncol = d)
  for (i in 1:d) {
    A[, i] <- stats::rnorm(n, mean, sd)
  }
  return(A)
}
