#' Simulate an example data with 100 observations for each trajectory
#' @param n1,n2,n3 number of subjects for each group
#' @return x is the observation time points; Y is the data saved in wide form; Z stores the labels
#' @export
generateToyDataset <- function(n1 = 40, n2 = 30, n3 = 30) {

  n <- n1 + n2 + n3

  m <- 100
  x <- seq.int(from = 0, to = 1, length.out = m)

  klas <- c(ones(n1, 1), 2 * ones(n2, 1), 3 * ones(n3, 1))

  sigma <- 0.1

  Ey1 <- function(x) {
    return(0.8 + 0.5 * exp(-1.5 * x) * sin(1.3 * pi * x))
  }

  Ey2 <- function(x) {
    return(0.5 + 0.8 * exp(-0.1 * x) * sin(0.9 * pi * x))
  }

  Ey3 <- function(x) {
    return(1 + 0.5 * exp(-x) * sin(-1.2 * pi * x))
  }

  Y <- rbind(t(Ey1(x) + sapply(1:n1, function(i) sigma * stats::rnorm(n = m))),
             t(Ey2(x) + sapply(1:n2, function(i) sigma * stats::rnorm(n = m))),
             t(Ey3(x) + sapply(1:n3, function(i) sigma * stats::rnorm(n = m))))

  return(list(x = x, Y = Y, Z = klas))
}
