#' Construct the spline regression matrix for spline regression model of order M
#' @param x predictors [nx1] column vector
#' @param M spline order (polynomial degree = M-1)
#' @param knots interior knots [Kx1] ([1xK]) (K the number of interior knots)
#' @return X: regression matrix of the model such that: Xij = [h_1(xi),..,h_j(xi),...,h_{M+K}(xi)] [nx(K+M)]
#' @noRd

splinebasis <- function(x, knots, M) {

  # X = splinebasis(x,knots,M)
  # construct the spline regression matrix  for spline regression model of
  # order M
  #
  # Inputs:
  #
  #       x: predictors [nx1] column vector
  #       M: spline order (polynomial degree = M-1)
  #       knots: interior knots [Kx1] ([1xK]) (K the number of interior knots)
  #
  # Outputs:
  #
  #       X: regression matrix of the model such that: Xij = [h_1(xi),..,h_j(xi),...,h_{M+K}(xi)] [nx(K+M)]
  #
  #######################################################################################################################

  p <- M - 1 # Polynomial degree
  K <- length(knots)

  X <- matrix(NA, length(x), K + M)
  for (ord in 0:p) {
    X[, ord + 1] <- x ^ ord # [1 t t.^2 t.^3 t.^p;......;...]
  }

  for (k in 1:K) {
    X[, k + M] <- sapply(X = x, FUN = function(x) max(x - knots[k], 0) ^ p)
  }

  return(X)
}
