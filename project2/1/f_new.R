#setwd("/Users/sondrerogde/Library/Mobile Documents/com~apple~CloudDocs/Documents/NTNU/8. semester/TMA4300")

# declaring parameters
alpha <- 2
beta <- 0.05
iters <- 50000

# loading data
#load(file = "project2/data/rain.rda")
y <- rain[, "n.rain"]
n <- rain[, "n.years"]

# expit function
expit <- function(x) {1/(1+exp(-x))}


mcmc_sampler <- function(alpha, beta, iters, burnin=ceiling(iters/10)){
  # declaring variables
  iters = iters + burnin
  
  # declaring initial data structures
  chain <- matrix(nrow=iters, ncol=367)
  acceptance_probs = matrix(nrow=iters-1, ncol=366)
  
  
  # calculating initial x-value
  x <- (y+1)/(n+1)
  pi <- expit(x)
  
  # sampling first sigma
  sigma_u_sq <- 1/rgamma(1, alpha, beta)
  
  # setting initial row for matrix
  chain[1,] <- c(x, sigma_u_sq)
  
  # outer loop used for sampling x-vector
  for (iter in 2:iters){
    # uniform sampling used for accepting proposals
    u <- runif(366)
    
    # vectors for new samples and previous samples respectively
    new_x = c()
    old_x = chain[iter-1, ]
    
    # Gibbs sampling for sigma
    sigma_alpha <- alpha + 365/2
    sigma_beta <- beta + 1/2*sum((old_x[2:366] - old_x[1:365])^2)
    sigma_u_sq <- 1/rgamma(1, shape=sigma_alpha, scale=sigma_beta)
    
    # inner loop used for sampling x days
    for (t in 1:366) {
      # if-statement handling the three cases where t=1, 1<t<366, and t=366
      if (t == 1){
        proposal <- rnorm(1, old_x[2], sd=sqrt(sigma_u_sq))
      } else if (t != 366) {
        proposal <- rnorm(1, (tail(new_x, n=1) + old_x[t+1])/2, sd=sqrt(sigma_u_sq/2))
      } else if (t == 366) {
        proposal <- rnorm(1, tail(new_x, n=1), sd=sqrt(sigma_u_sq))
      }
      
      # calculating new and old pi
      new_pi = expit(proposal)
      pi <- expit(old_x[t])
      
      # calculating acceptance probability (alpha)
      acceptance_prob <- min(1, dbinom(y[t], n[t], new_pi)/dbinom(y[t], n[t], pi))
      
      # accepting proposal based on acceptance probability
      new_x <- c(new_x, ifelse(u[t] < acceptance_prob, proposal, old_x[t]))
      
      # used to track observed acceptance probabilities
      acceptance_probs[iter-1, t] <- (u[t] < acceptance_prob)
    }
    
    # updates chain and time data structures
    chain[iter, ] <- c(new_x, sigma_u_sq)
    
    if (iter %% 100 == 0){print(iter/iters)}
  }
  
  # calculates mean acceptance probabilities
  acceptance_probs <- colMeans(acceptance_probs)
  
  # declares list to return
  result_list <- list(
    chain = chain,
    acceptance_probs = acceptance_probs
  )
  return (result_list)
}
  
# calls the function and measures the runtime
start_time <- proc.time()[3]
results <- mcmc_sampler(alpha, beta, iters)
end_time <- proc.time()[3]
total_time = end_time - start_time
print(total_time)

# sets variables for returned values
chain <- results$chain
acceptance_probs <- results$acceptance_probs

burnin = ceiling(iters/10)

# traceplots of x_1, x_201, x_366, and sigma
plot(chain[,1], type="l", ylab="x_1")
title("Traceplot of x_1 (including burnin)")

plot(chain[,201], type="l", ylab="x_201")
title("Traceplot of x_201 (including burnin)")

plot(chain[,366], type="l", ylab="x_366")
title("Traceplot of x_366 (including burnin)")

plot(chain[,367], type="l", ylab="σ²")
title("Traceplot of σ² (including burnin)")

# plot of acceptance probabilities for all days
plot(acceptance_probs, type="l")
abline(h = mean(acceptance_probs), col = "red", lty = 2)
legend("bottomright", 
       legend = c("acceptance probabilities", "mean"), 
       col = c("black", "red"), 
       lty = 1, 
       bty = "n")
title("Acceptance probabilities for each day")
print("Mean acceptance prob:")
print(mean(acceptance_probs))

library(ggplot2)
library(coda)

# calculates probabilities
prob_chain <- expit(chain[, 1:366])

# creates a mcmc-object of the chain (excluding sigma)
mcmc_chain <- as.mcmc(prob_chain)

# calculates confidence intervals for sigma
mcmc_sigma <- as.mcmc(chain[, 367])
conf_sigma <- summary(mcmc_sigma)$quantile
print(mean(chain[, 367]))
print(conf_sigma[c(1, 5)])


# calculates confidence intervals
conf_levels = summary(mcmc_chain)$quantile

# calculates mean values of each day
mean_pixt <- colMeans(mcmc_chain[, 1:366])
t_values <- 1:366

# creates a dataframe with time-values, mean x-values and confidence levels
df <- data.frame(t=t_values, mean=mean_pixt, lwr=conf_levels[,1], upr=conf_levels[,5])
# names(df) <- c("Day", "Mean", "Lower", "Upper")

# plots the mean x-values and confidence intervals with y-values from the dataset
ggplot(df, aes(x = t)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill="Credible Intervals"), alpha = 0.2) +
  geom_line(aes(y = mean_pixt, color = "Mean probability")) +
  geom_line(aes(y = y / n, color = "Observed Ratio"), alpha=0.5) +
  labs(title = "Posterior Mean of π(x_t) with 95% Credible Intervals",
       x = "t", y = "π(xt)", color = "Legend") +
  scale_color_manual(values = c("Mean probability" = "blue", "Observed Ratio" = "red")) +
  scale_fill_manual(name = "Credible Intervals", values = "blue") +  # Added for ribbon fill color
  theme(legend.position = "bottom")


  # plots histograms of x_1, x_201, x_366 and σ²
hist(expit(chain[burnin + 1:iters, 1]), main="Histogram of π(x_1) (excluding burnin)", xlab="π(x_1)", breaks=50)
hist(expit(chain[burnin + 1:iters, 201]), main="Histogram of π(x_201) (excluding burnin)", xlab="π(x_201)", breaks=50)
hist(expit(chain[burnin + 1:iters, 366]), main="Histogram of π(x_366) (excluding burnin)", xlab="π(x_366)", breaks=50)
hist(chain[burnin + 1:iters, 367], main="Histogram of σ² (excluding burnin)", xlab="σ²", breaks=50)


# plots autocorrelation of x_1, x_201, x_366 and σ²
acf(expit(chain[burnin + 1:iters, 1]), main="Autocorrelation of π(x_1) (excluding burnin)")
acf(expit(chain[burnin + 1:iters, 201]), main="Autocorrelation of π(x_201) (excluding burnin)")
acf(expit(chain[burnin + 1:iters, 366]), main="Autocorrelation of π(x_366) (excluding burnin)")
acf(chain[burnin + 1:iters, 367], main="Autocorrelation of σ² (excluding burnin)")

