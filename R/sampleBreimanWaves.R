#' @export
sampleBreimanWaves <- function(n) {
  K <- 3

  Pik <- 1 / K * ones(K, 1) # mixing proportions

  dt <- 1 # periode d'echantilonnage (secondes)
  fs <- 1 / dt
  temps <- seq.int(from = 0, to = 21 - dt, by = dt) # temps;
  m <- length(temps)

  X <- zeros(n, m)
  klas <- zeros(n, 1)


  h0 <- sapply(temps, function(x) max(6 - abs(x - 11), 0))

  Dt <- 4 / dt

  h1 <- c(zeros(1, Dt), h0[1:(length(h0) - Dt)])
  h2 <- c(h0[(Dt + 1):length(h0)], zeros(1, Dt))


  for (i in 1:n) {
    mu <- rand(1, m)
    ei <- stats::rnorm(n = m)

    zik <- stats::rmultinom(n = 1, size = 1, prob = Pik)
    zi <- which(zik == 1)
    klas[i] <- zi

    if (zi == 1) {
      X[i, ] <- mu * h0 + (1 - mu) * h1 + ei
    } else if (zi == 2) {
      X[i, ] = mu * h0 + (1 - mu) * h2 + ei
    } else {
      X[i, ] <- mu * h1 + (1 - mu) * h2 + ei
    }
  }

  return(list(Y = X, Z = klas))

}
