

a_constant <- function(x) {
  n <- length(x)
  
  #MLE is again the mean
  beta_hat <- mean(x)
  
  #compute MLEs excludingn each data point
  estimates_excluding_i <- sapply(1:n, function(i) mean(x[-i]))
  
  #Computign the expression for a in 9.3.2.1 in Givens & Hoeting
  numerator <- sum((estimates_excluding_i - beta_hat)^3)
  denominator <- sum((estimates_excluding_i - beta_hat)^2)^1.5
  
  return(a <- numerator / (6 * denominator))
}


b_constant <- function(x) {
  n <- length(x)
  
  # MLE of beta is just mean
  beta_hat <- mean(x)
  
  # Evaluate bootstrap cdf at the MLE
  p <- pgamma(beta_hat, shape = n, scale = beta_hat / n)
  
  # Compute b using the inverse normal CDF
  b <- qnorm(p)
  
  return(b)
}

# Example usage:
# x is your vector of observations from an exponential distribution


calculate_interval_for_beta <- function(x, alpha) {
  #length and MLE
  n <- length(x)
  beta.hat <- mean(x)
  
  #compute constants based on functions alreay defined
  b <- b_constant(x)
  a <- a_constant(x)

  # Normal z-scores
  z_alpha <- qnorm(alpha / 2)
  z_one_minus_alpha <- qnorm(1 - alpha / 2)
  
  # Adjusted quantiles for the BCa interval
  beta_1 <- pnorm(b + (b + z_alpha) / (1 - a * (b + z_alpha)))
  beta_2 <- pnorm(b + (b + z_one_minus_alpha) / (1 - a * (b + z_one_minus_alpha)))
  
  #computing the lower and upper bound of the interval
  lower_bound = mean(x)*qchisq(beta_1, 2*n)/(2*n)
  upper_bound = mean(x)*qchisq(beta_2, 2*n)/(2*n)
  
  return(c(lower_bound, upper_bound))
}

# Function to simulate one experiment and check if the true beta is within the BCa interval
simulate_one_experiment <- function(true_beta, sample_size, alpha) {
  # Generate a sample
  x <- rexp(sample_size, 1 / true_beta)
  
  # Calculate the BCa confidence interval
  BCa_quantiles <- calculate_interval_for_beta(x, alpha)
  lower = BCa_quantiles[1]
  upper =BCa_quantiles[2]
  
  # Check if the true beta is within the BCa interval
  return(true_beta >= lower && true_beta <= upper)
}

# Set parameters
true_beta <- 1
sample_size <- 100
alpha <- 0.05
num_simulations <- 10000

# Simulate
results <- replicate(num_simulations, simulate_one_experiment(true_beta, sample_size, alpha))

# Get coverage and print
coverage <- mean(results)
cat("Estimated coverage of the BCa interval: ", coverage)

