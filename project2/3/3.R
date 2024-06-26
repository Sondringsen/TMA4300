# Import RTMB
library("RTMB")

# Load the rain data
load(file="/Users/eiriksteen/Personal/school/compstat/projects/TMA4300/project2/data/rain.rda")

# Define the expit function
sigmoid <- function(x) {
  return(1 / (1 + exp(-x)))
}

# Define the negative log likelihood
f <- function(p) {

  y <- rain$n.rain
  n <- rain$n.years

  # Compute the binomial terms of the log likelihood
  bin_terms <- -sum(y * log(sigmoid(p$x))) -
    sum((n - y) * log(1 - sigmoid(p$x)))

  # Compute the random walk terms of the log likelihood
  rw_terms <- (length(p$x) - 1) * log(1 / (p$prec ** (1 / 2))) +
    (p$prec / 2) * sum((p$x[-1] - p$x[-length(p$x)])**2)

  return(bin_terms + rw_terms)
}

# Set the random effects to 0
x <- rain$n.rain * 0
# Set the variance
prec <- 0.05

# Invoke the autodifferentiator
obj <- MakeADFun(f, list(x = x, prec = prec), random = "x", silent = TRUE)
# Obtain the maximum likelihood estimates
opt <- nlminb(obj$par, obj$fn, obj$gr)

# sdreport(obj)