0


bias <- estimate_bias(mod, sample)

bias_corrected_estimates <- mod$coef + bias

# Problem 1.d

percentile_method <- function(sample_vec, alpha) {
  # returns the quantiles P(q_1 <= theta <= q_2) = 1 - alpha
  sorted_vec = sort(sample_vec)
  q1_idx = floor(length(sample_vec)*alpha/2)
  q2_idx = ceiling(length(sample_vec)*(1-alpha)/2)
  return (c(sorted_vec[q1_idx], sorted_vec[q2_idx]))
}

alpha = 0.05
theta_0_int = percentile_method(sample[ , 1], alpha)
theta_1_int = percentile_method(sample[ , 2], alpha)

var(sample)
