Monte_Carlo <- function(n) {

#Sample standard normal variables using Box-Muller transform
U1 <- runif(n)
U2 <- runif(n)
Z <- sqrt(-2 * log(U1)) * cos(2 * pi * U2)

if (SHOW_HISTOGRAMS){
  # Compare with actual pdf
  hist(Z, breaks=50, probability=TRUE, col="lightblue",
       main="Histogram of Simulated Standard Normal Variables",
       xlab="Value", ylab="Density")
  curve(dnorm, col="darkred", add=TRUE, lwd=2)
  legend("topright", legend=c("Simulated Data", "Normal PDF"), 
         fill=c("lightblue", "darkred"))
}

MC_estimate <- mean(Z > 4)
print("Monte-Carlo estimate: ")
print(MC_estimate)

print("Confidence interval:")
print(t.test(Z > 4)$conf)
}

Importance_Sampling <- function(n) {
#Importance sampling estimate

pdf_proposal <- function(x) {
  ifelse(x < 4, 0, x*exp(8 - x^2 / 2))
}

# Generate samples using inversion sampling
uniform_samples <- runif(n) 
inverse_cdf <- function(x) {
  sqrt(16 - 2 * log(1 - x))
}
samples <- inverse_cdf(uniform_samples)

if (SHOW_HISTOGRAMS){
# Verify that the inversion sampling worked as expected
hist(samples, breaks=50, probability=TRUE, col="lightblue",
     main="Histogram vs. Actual PDF", xlab="Value", ylab="Density")
curve(pdf_proposal(x), from=min(samples),
      to=max(samples), add=TRUE, col="darkred", lwd=2)
legend("topright", legend=c("Generated Samples", "Actual PDF"),
       col=c("lightblue", "darkred"), lwd=2, fill=c("lightblue", "darkred"))
}

# Importance sampling estimate
weights <- dnorm(samples) / pdf_proposal(samples)  
estimate <- mean(weights * (samples > 4))

print("Importance sampling estimate:")
print(estimate)
print("Confidence interval:")
print(t.test(weights * (samples > 4))$conf)
}

antithetic_Importance_Sampling <- function(n) {

  pdf_proposal <- function(x) {
    ifelse(x < 4, 0, x*exp(8 - x^2 / 2))
  }
  
  # Generate samples using both U and 1-U as inputs to inverse cdf
  uniform_samples <- runif(n) 
  uniform_samples <- c(uniform_samples, 1-uniform_samples)
  inverse_cdf <- function(x) {
    sqrt(16 - 2 * log(1 - x))
  }
  samples <- inverse_cdf(uniform_samples)
  
  if (SHOW_HISTOGRAMS){
    # Verify that samples have the right density
    hist(samples, breaks=50, probability=TRUE, col="lightblue", 
         main="Histogram vs. Actual PDF", xlab="Value", ylab="Density")
    curve(pdf_proposal(x), from=min(samples), to=max(samples), 
          add=TRUE, col="darkred", lwd=2)
    legend("topright", legend=c("Generated Samples", "Actual PDF"), 
           col=c("lightblue", "darkred"), lwd=2, fill=c("lightblue", "darkred"))
  }
  
  # Importance sampling
  weights <- dnorm(samples) / pdf_proposal(samples)  
  estimate <- mean(weights * (samples > 4))
  
  print("Antithetic + importance sampling estimate:")
  print(estimate)
  print("Confidence interval:")
  print(t.test(weights * (samples > 4))$conf)
}

Monte_Carlo_with_chunking <- function(chunk_size, n_chunks) {
  # Initialize counters
  count_greater_than_4 <- 0
  
  # Sample standard normal variables using Box-Muller transform
  for (i in 1:n_chunks) {
    U1 <- runif(chunk_size)
    U2 <- runif(chunk_size)
    Z <- sqrt(-2 * log(U1)) * cos(2 * pi * U2)
    count_greater_than_4 <- count_greater_than_4 + sum(Z > 4)
  }
  
  n =  (chunk_size*n_chunks)
  
  # Calculate Monte-Carlo estimate
  MC_estimate <- count_greater_than_4 / n
  print("Monte-Carlo estimate: ")
  print(MC_estimate)
  
  # Using Normal Approximation for binomial proportion confidence interval
  p_hat <- MC_estimate
  z <- qnorm(0.975) 
  margin_of_error <- z * sqrt(p_hat * (1 - p_hat) / n)
  ci_lower <- p_hat - margin_of_error
  ci_upper <- p_hat + margin_of_error
  
  print("Confidence interval (95%):")
  print(c(ci_lower, ci_upper))
}

N <- 1e7 #number of samples
SHOW_HISTOGRAMS <- FALSE
antithetic_Importance_Sampling(N)


