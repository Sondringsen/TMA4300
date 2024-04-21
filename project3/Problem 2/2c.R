compute_exact_coverage <- function(n, alpha) {
  
  # Chi-square quantiles
  chi2_lower <- qchisq(alpha/2, 2*n, lower.tail = TRUE)
  chi2_upper <- qchisq(1 - alpha/2, 2*n, lower.tail = TRUE)
  
  # Get probabilities of being inside/outside the interval
  F_lower <- pchisq(4*n^2 / chi2_lower, 2*n, lower.tail = TRUE)
  F_upper <- pchisq(4*n^2 / chi2_upper, 2*n, lower.tail = TRUE)
  
  # Compute the coverage
  result <- 1 - (F_upper + (1 - F_lower))
  return(result)
}

# Test for different values of n
n_values <- c(5, 10, 20, 50, 100)  
alpha <- 0.05  

for (n in n_values) {
  result <- compute_exact_coverage(n, alpha)
  cat(sprintf("%d: %f\n", n, result))
}