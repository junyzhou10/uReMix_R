#' log of determinant of matrix A by Choleski decomposition
#' @param A input matrix
#' @noRd

logdet <- function(A) {
  # log(det(A)) where A is positive-definite.
  # This is faster and more stable than using log(det(A)).

  # From Tom Minka's lightspeed toolbox

  U <- chol(A)
  y <- 2 * sum(log(diag(U)))

  return(y)
}
