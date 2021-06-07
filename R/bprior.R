bprior <- function(gamma_hat) {
  m <- mean(gamma_hat)
  v <- var(gamma_hat)
  output <- (m * v + m^3) / v
  return(output)
}
