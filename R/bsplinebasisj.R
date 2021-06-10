#' Compute a bspline basis function of order m.
#' @param x Vector of points at which to evaluate the b-spline basis
#' @param t Vector of knots. Normally the first and last will be repeated m times
#' @param j Scalar specifying which basis function from the basis function
#' @param M order of B-spline. Default m = 4 (cubics since degree= m-1 = 3).
#' @return Bj: Basis function evaluated at the points in x.
#' @noRd
#'

bsplinebasisj <- function(x, t, j, M = 4) {

  # bsplineBasis: compute a bspline basis function of order m.
  #
  #    Bj = bsplineBasis(x, t, j, m)
  #
  # Evaluates the j^th B-spline basis function of degree m for a given set of
  # knots t at points x using the Cox-de Boor recursion formula.
  #
  # Inputs:
  #
  #   x: Vector of points at which to evaluate the b-spline basis.
  #   t: Vector of knots.  Normally the first and last will be repeated n times.
  #   j: Scalar specifying which basis function from the basis function
  #   set.  1 <= j <= length(t) - m - 1.
  #   m: Order of basis functions.  Default m = 4 (cubics since degree= m-1 = 3).
  #
  # Outputs:
  #
  #   Bj: Basis function evaluated at the points in x.
  #
  #################################################################################

  p <- M - 1 # polynomial degree

  # There are K knots. # j=1:L = K+M
  L <- length(t) # L=K+M

  # Check validity of j.
  if ((j < 1) || (j > L - M)) {
    stop("Parameter j = ", j, " is not in range [1, ", L - M, "]!")
  }

  # Construct the b-spline recursively.
  if (p == 0) { # piecewise constant
    # Base case.
    # if(j+2 == length(t))
    #  y = ((x >= t(j+1)) & (x <= t(j+2)));
    # else
    Bj <- ((x >= t[j]) & (x < t[j + 1])) * 1

  } else {# If the two knots in the denominator are equal, we ignore the term (to not divide by zero).

    denom1 <- t[j + M - 1] - t[j]
    if (denom1 == 0) {
      Bj <- rep.int(x = 0, times = length(x))
    } else {
      Bj <- (x - t[j]) * bsplinebasisj(x, t, j, M - 1) / denom1
    }

    denom2 <- t[j + M] - t[j + 1]
    if (denom2 != 0) {
      Bj <- Bj + (t[j + M] - x) * bsplinebasisj(x, t, j + 1, M - 1) / denom2
    }
  }

  return(Bj)
}
