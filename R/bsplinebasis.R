#' Construct the B-spline regression matrix of order M
#' @param x Vector of points at which to evaluate the b-spline basis
#' @param t Vector of knots. Normally the first and last will be repeated m times
#' @param M order of B-spline (degree of poynomial pieces p = m-1)
#' @return B: the B-spline Basis (yfit = B*coeff)
#' @noRd

bsplinebasis <- function(x, t, M) {

  ############################################################################
  #
  # B = bsplinebasis(x, t, M) : construct the spline regression matrix  for Bspline regression model of
  # order M
  #
  # Inputs:
  #
  #      x: Vector of points at which to evaluate the b-spline basis
  #      t: Vector of knots.  Normally the first and last will be repeated m times
  #         here t represents the knots sequence (including the two boundary knots)
  #      M: order of B-spline (degree of poynomial pieces p = m-1)
  # Outputs:
  #
  #      B: the B-spline Basis (yfit = B*coeff)
  #
  #
  #
  #
  #
  #
  #
  # # t is tau in the document
  ##############################################################################

  m <- length(x)

  # Repeat the first and the last knots m times
  t <- c(rep(t[1], M - 1), t, rep(t[length(t)], M - 1))

  # j th basis function : 1 <= j <= length(t) - M -1 ; There are K + 2*M knots ; M = length(t) - 1;
  B <- zeros(m, length(t) - M)

  for (j in 1:(length(t) - M)) {# is the same as j=1 : K+M, K is the number of interior knots
    B[ , j] <- bsplinebasisj(x, t, j, M) # j-1 0 <= j <= length(t) - n - 2.
  }

  return(B)
}
