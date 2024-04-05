load(file=url("https://www.math.ntnu.no/emner/TMA4300/2024v/data.Rdata"))

mod <- glm(cbind(y, m - y) ~ x, family = binomial, data = data)

# names(mod)

print(mod$coef)
print(vcov(mod))


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

estimate_var <- function(sample){

  means <- colMeans(sample)
  vars <- 1 / (nrow(sample) - 1) * colSums((sample - means)**2)

  return(vars)
}

var <- estimate_var(sample)

estimate_bias <- function(model, sample) {

  return(colMeans(sample) - model$coef)
}

bias <- estimate_bias(mod, sample)

bias_corrected_estimates <- mod$coef + bias