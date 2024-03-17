# Import RTMB
library("RTMB")

# Load the rain data
# load(file="/Users/eiriksteen/Personal/school/compstat/projects/TMA4300/project2/data/rain.rda")
load(file = "project2/data/rain.rda")

# Define the expit function
sigmoid <- function(x) {
  return(1 / (1 + exp(-x)))
}

# logit function
logit <- function(x) {1/(1+exp(-x))}

# Define the negative log likelihood
f <- function(p) {

  y <- rain$n.rain
  n <- rain$n.years

  # Compute the binomial terms of the log likelihood
  bin_terms <- -sum(y * log(sigmoid(p$x))) -
    sum((n - y) * log(1 - sigmoid(p$x)))

  # Compute the random walk terms of the log likelihood
  rw_terms <- (length(p$x) - 1) * log(p$var ** (1 / 2)) +
    (1 / (2 * p$var)) * sum((p$x[-1] - p$x[-length(p$x)])**2)

  return(bin_terms + rw_terms)
}

# Set the random effects to 0
x <- rain$n.rain
# Set the variance
var <- 0.05

# Invoke the autodifferentiator
obj <- MakeADFun(f, list(x = x, var = var), silent = TRUE)
# Obtain the maximum likelihood estimates
opt <- nlminb(obj$par, obj$fn, obj$gr)

parameters <- sdreport(obj)
t_values <- 1:366
mean_probs <- logit(parameters$par.fixed[1:366])

# creates a dataframe with time-values, mean x-values and confidence levels
df <- data.frame(t = t_values, mean=mean_probs)

# plots the mean x-values with y-values from the dataset
ggplot(df, aes(x = t)) +
  geom_line(aes(y = mean, color = "Mean probability")) +
  geom_line(aes(y = y / n, color = "Observed Ratio"), alpha=0.5) +
  labs(title = "Posterior Mean of π(x_t) with 95% Credible Intervals",
       x = "t", y = "π(xt)", color = "Legend") +
  scale_color_manual(values = c("Mean probability" = "blue", "Observed Ratio" = "red")) +
  scale_fill_manual(name = "Credible Intervals", values = "blue") +  # Added for ribbon fill color
  theme(legend.position = "bottom")
