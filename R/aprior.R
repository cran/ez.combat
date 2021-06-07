aprior <- function(gamma_hat) {
  m <- mean(gamma_hat)
  v <- var(gamma_hat)
  output <- (2 * v + m^2) / v
  return(output)
}
