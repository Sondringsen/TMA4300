u <- scan(file="https://www.math.ntnu.no/emner/TMA4300/2024v/u.txt")
z <- scan(file="https://www.math.ntnu.no/emner/TMA4300/2024v/z.txt")

update_lambdas <- function(l0, l1) {
  n <- length(z)
  new_l0 <- n / sum(u * z + (1 - u) * (1 / l0 - z / (exp(l0 * z) - 1)))
  new_l1 <- n / sum((1 - u) * z + u * (1 / l1 - z / (exp(l1 * z) - 1)))
  return(c(new_l0, new_l1))
}

run_em <- function(n_iter, plot = TRUE) {
  l0 <- 1
  l1 <- 1

  l0s <- c(l0)
  l1s <- c(l1)

  for (i in 1:n_iter) {
    out <- update_lambdas(l0, l1)
    l0 <- out[1]
    l1 <- out[2]

    l0s <- append(l0s, l0)
    l1s <- append(l1s, l1)
  }

  if (plot) {
    matplot(1:(n_iter + 1), cbind(l0s, l1s), type = "l")
  }

  return(c(l0, l1))
}

out <- run_em(7)