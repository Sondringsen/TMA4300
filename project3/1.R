load(file=url("https://www.math.ntnu.no/emner/TMA4300/2024v/data.Rdata"))

mod <- glm(cbind(y, m - y) ~ x, family = binomial, data = data)


print(mod$coef)
print(vcov(mod))

# Problem 1.a

boot <- function(model, B = 10000) {

  out <- c()

  for (b in 1:B){
    idx <- sample.int(nrow(data), nrow(data), replace = TRUE)
    bootdata <- data[idx, ]
    model <- update(model, data = bootdata)
    out <- append(out, model$coef)
  }

  return(t(matrix(out, nrow = 2)))

}

sample <- boot(mod)

# Problem 1.b

estimate_var <- function(sample){

  means <- colMeans(sample)
  vars <- 1 / (nrow(sample) - 1) * colSums((sample - means)**2)

  return(vars)
}

var <- estimate_var(sample)

# Problem 1.c

estimate_bias <- function(model, sample) {

  return(colMeans(sample) - model$coef)
}

bias <- estimate_bias(mod, sample)

bias_corrected_estimates <- mod$coef + bias

# Problem 1.d

percentile_method <- function(sample_vec, alpha) {
  # returns the quantiles P(q_1 <= theta <= q_2) = 1 - alpha
  print(sample_vec)
  sorted_vec = sort(sample_vec)
  print(sorted_vec)
  q1_idx = floor(length(sample_vec)*alpha/2)
  q2_idx = ceiling(length(sample_vec)*(1-alpha)/2)
  return (c(sorted_vec[q1_idx], sorted_vec[q2_idx]))
}

alpha = 0.05
theta_0_int = percentile_method(sample[ , 1], alpha)
theta_1_int = percentile_method(sample[ , 2], alpha)
