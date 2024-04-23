u <- scan(file = "https://www.math.ntnu.no/emner/TMA4300/2024v/u.txt")
z <- scan(file = "https://www.math.ntnu.no/emner/TMA4300/2024v/z.txt")

update_lambdas <- function(l0, l1, bootstrap = FALSE) {
  # Update lambdas with the argmax from the M-Step
  n <- length(z)
  x <- u
  y <- z
  if (bootstrap) {
    x <- sample(x, replace = TRUE)
    y <- sample(y, replace = TRUE)
  }
  new_l0 <- n / sum(x * y + (1 - x) * (1 / l0 - y / (exp(l0 * y) - 1)))
  new_l1 <- n / sum((1 - x) * y + x * (1 / l1 - y / (exp(l1 * y) - 1)))

  return(c(new_l0, new_l1))
}

run_em <- function(n_iter, plot = TRUE, bootstrap = FALSE) {
  l0 <- 1
  l1 <- 1

  # Init arrays for storing lambdas per time step
  l0s <- c(l0)
  l1s <- c(l1)

  for (i in 1:n_iter) {
    # Update lambdas
    out <- update_lambdas(l0, l1, bootstrap = bootstrap)
    l0 <- out[1]
    l1 <- out[2]
    l0s <- append(l0s, l0)
    l1s <- append(l1s, l1)
  }

  # Plot lambdas per time step
  if (plot) {
    png("./plots/convergence.png")
    plot(1:(n_iter + 1), l1s, col="red", type = "b", xlab = "niter", ylab = "lambda")
    lines(1:(n_iter + 1), l0s, col="blue", type = "b")
    legend(5, 7, legend = c("Lambda 1", "Lambda 0"), 
    col = c("red", "blue"), lty = 1:2, cex = 0.8)
    dev.off()
  }

  return(c(l0, l1))
}

bootstrap <- function(n_iter) {
  samples <- c()

  # Run the EM algorithm n_iter times with bootstrap samples
  for (i in 1:n_iter){
    # Run the EM. We select 20 iters for the EM algorithm, which we
    # empirically observed to be past convergence
    samples <- append(samples, run_em(20, plot = FALSE, bootstrap = TRUE))
  }
  # Reshape samples to matrix of shape (n_iter, 2)
  samples <- t(matrix(samples, 2, n_iter))

  return(samples)
}

# Calculate bootstrap sample estimates
samples <- bootstrap(n_iter = 1000)
# Compute standard dev estimate
sds <- apply(samples, 2, sd)
# Compute correlation estimate
co <- cor(samples[, 1], samples[, 2])
# Compute mean estimate
means <- colMeans(samples)
# Compute non bootstrapped estimate
estimate <- run_em(20, plot = FALSE, bootstrap = FALSE)
# Compute bias estimate
bias <- means - estimate

out <- run_em(20, bootstrap = FALSE)
print(out)
print("***")
print("Boot standard deviations")
print(sds)
print("Boot correlation")
print(co)
print("Boot means")
print(means)
print("Estimate (non boot)")
print(estimate)
print("Boot bias")
print(bias)

calculate_mle <- function() {
  n <- length(u)
  # Calculate the lambda_0 estimate
  l0 <- (4*sum(u)-2*length(u))/(2*sum(u*z)-sum(z))
  # Calculate the lambda_1 estimate
  l1 <- (4*n*sum(u)*sum(z)-n^2*sum(z)-4*sum(u)**2*sum(z))/(n*sum(u*z)*sum(z)+2*sum(u)*sum(u*z)*sum(z)-sum(u)*sum(z)**2-2*n*sum(u*z)**2)
  return(c(l0, l1))
}

mle <- calculate_mle()
print("MLE")
print(mle)