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
  T <- length(y)

  # bin_terms <- - sum(y * log(sigmoid(p$x))) - sum((rain$n.years - y) * log(1 - sigmoid(p$x)))

  # rw_terms <- (length(y) - 1)*log(p$var**(1/2)) + (1/(2*p$var)) * sum((p$x[-1] - p$x[-length(p$x)])**2)

  bin_terms <- -sum(y * log(sigmoid(p$x))) - sum((n - y) * log(1 - sigmoid(p$x)))
  rw_terms <- (T-1) * log(p$var ** (1/2)) + (1/(2*p$var)) * sum((p$x[-1] - p$x[-length(p$x)])**2)

  return(bin_terms + rw_terms)
}

# Compute the true rain fraction observed in the data
# add 1e-04 to account for the rainless days
pi <- rain$n.rain / rain$n.years + 1e-04
# Compute x with the logit function
x <- log(pi / (1 - pi))
# Set the variance
var <- 0.05

# Invoke the autodifferentiator
obj <- MakeADFun(f, list(x = x, var = var), silent = TRUE)
# Obtain the maximum likelihood estimates
opt <- nlminb(obj$par, obj$fn, obj$gr)