# setwd("/Users/sondrerogde/Library/Mobile Documents/com~apple~CloudDocs/Documents/NTNU/8. semester/TMA4300")
library(coda)

alpha <- 2
beta <- 0.05
iters <- 50000
conf_level <- 0.95
load(file = "project2/data/rain.rda")
y = rain[, "n.rain"]
n = rain[, "n.years"]

logit <- function(x) {1/(1+exp(-x))}


mcmc_sampler <- function(alpha, beta, iters, conf_level, burnin=ceiling(iters/10)){
  # MH sampling
  iters = iters + burnin
  x <- (y+1)/(n+1)
  pi <- logit(x)
  n_accepted <- 0
  sigma_u_sq <- 1/rgamma(1, alpha, beta)
  chain <- matrix(nrow=iters, ncol=367)
  time <- rep(NA, iters)
  
  chain[1,] <- c(x, sigma_u_sq)
  time[1] <- proc.time()[3]
  
  
  for (iter in 2:iters){
    u <- runif(366)
    
    sigma_alpha <- alpha + 1 + 365/2
    sigma_beta <- beta + 1/2*sum((x[2:366] - x[1:365])^2)
    sigma_u_sq <- 1/rgamma(1, sigma_alpha, sigma_beta)
    
    
    first_proposal <- rnorm(1, x[2], sqrt(sigma_u_sq))
    last_proposal <- rnorm(1, x[365], sqrt(sigma_u_sq))
    middle_propsal <- rnorm(364, (x[3:366]+x[1:364])/2, sqrt(sigma_u_sq/2))
    proposal <- c(first_proposal, middle_propsal, last_proposal)
    
    new_pi = logit(proposal)
    acceptance_prob <- pmin(1, dbinom(y, n, new_pi)/dbinom(y, n, pi))
    new_x <- ifelse(u < acceptance_prob, proposal, chain[iter-1,])
    pi <- logit(new_x)
    
    n_accepted <- n_accepted + sum(u < acceptance_prob)
    
    elapsed_time <- proc.time()[3]
    chain[iter, ] <- c(new_x, sigma_u_sq)
    time[iter] <- elapsed_time
    
  }
  
  acceptance_prob <- n_accepted/(366*iters)
  
  burnin = burnin + 1  
  result_list <- list(
    chain = chain[burnin: iters,],
    time = time,
    acceptance_prob = acceptance_prob
  )
  return (result_list)
  # initiate x0
  # iterate 50,000 times
  # generate y
  # calculate alpha
  # generate unif(0,1)
  # if u < alpha: x_i = y
  # else x_i = x_{i-1}
  # return x vector (remember burnin)
  
}


results = mcmc_sampler(alpha, beta, iters, conf_level)

chain <- results$chain
time <- results$time
acceptance_prob <- results$acceptance_prob

print(acceptance_prob)

# Assuming 'chain' is your MCMC matrix with each column being a variable
# and each row being an iteration
plot(chain[,201])
plot(y, type="l")
plot(colMeans(chain), axis=0, type="l")
plot(chain[, 366])
plot(chain[, 367])


# library(ggplot2)
# 
# # Convert the matrix to an mcmc object
mcmc_chain <- as.mcmc(chain)
# print(typeof(mcmc_chain))
# # Traceplots for all variables
# traceplot(mcmc_chain, col=1)
# traceplot(mcmc_chain, col=201)
# traceplot(mcmc_chain, col=366)
# traceplot(mcmc_chain, col=367)
# 
# 
# # Histograms for a specific variable (e.g., σ^2)
# hist(mcmc_chain[, "sigma2"], main="Histogram of σ²", xlab="σ²", breaks=50)
# 
# # Estimated autocorrelation for a specific variable (e.g., σ^2)
# acf(chain[, 367], main="Autocorrelation of σ²")
# 
# # Central estimates and 95% credible intervals for a specific variable
# summary(mcmc_chain[, 1])
# 
# # Calculate and plot 95% credible intervals for π(xt)
# # Assuming that π(xt) is in a column named 'pixt'
# # Calculate the credible interval
credible_interval <- apply(mcmc_chain[, 1:366], 2, function(x) quantile(x, probs=c(0.025, 0.975)))
# 
# # Plotting posterior mean with 95% credible intervals
# # Assuming that column 'pixt' contains π(xt)
mean_pixt <- colMeans(mcmc_chain[, 1:366])
t_values <- 1:366 # Change this if your t values differ
df <- data.frame(t=t_values, mean_pixt, lwr=credible_interval[1,], upr=credible_interval[2,])

ggplot(df, aes(x=t)) +
  geom_ribbon(aes(ymin=lwr, ymax=upr), fill="blue", alpha=0.2) +
  geom_line(aes(y=mean_pixt), color="blue") +
  geom_line(aes(y=y), color="red") +
  labs(title="Posterior Mean of π(xt) with 95% Credible Intervals",
       x="t", y="π(xt)") +
  theme_minimal()
# 
# # Use proc.time to get computation time
# start_time <- proc.time()
# # ... (place your MCMC code here)
# end_time <- proc.time()
# computation_time <- end_time - start_time
# computation_time[3] # Elapsed time
