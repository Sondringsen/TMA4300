alpha <- 2
beta <- 0.05
iters <- 22000
conf_level <- 0.95
load(file = "project2/data/rain.rda")
y <- rain[, "n.rain"]
n <- rain[, "n.years"]

M <- 20

logit <- function(x) {1/(1+exp(-x))}


mcmc_sampler <- function(alpha, beta, iters, conf_level, burnin=ceiling(iters/10), M){
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
    
    new_x = c()
    old_x = chain[iter-1, ]
    
    sigma_alpha <- alpha + 1 + 365/2
    sigma_beta <- beta + 1/2*sum((old_x[2:366] - old_x[1:365])^2)
    sigma_u_sq <- 1/rgamma(1, sigma_alpha, sigma_beta)
    
    for (t in 1:ceil(366/M)) {
      if (t == 1){
        # calculate the Qs
        proposal <- rnorm(1, old_x[2], sqrt(sigma_u_sq))
      } else if (t != ceil(366/M)) {
        # calculate the Qs
        proposal <- rnorm(1, (tail(new_x, n=1) + old_x[t+1])/2, sqrt(sigma_u_sq))
      } else if (t == ceil(366/M)) {
        # calculate the Qs
        proposal <- rnorm(1, tail(new_x, n=1), sqrt(sigma_u_sq))
      }
      
      new_pi = logit(proposal)
      pi <- logit(old_x[t])
      acceptance_prob <- min(1, dbinom(y[t], n[t], new_pi)/dbinom(y[t], n[t], pi))
      new_x <- c(new_x, ifelse(u[t] < acceptance_prob, proposal, old_x[t]))
      n_accepted <- n_accepted + (u[t] < acceptance_prob)
    }
    
    
    elapsed_time <- proc.time()[3]
    chain[iter, ] <- c(new_x, sigma_u_sq)
    time[iter] <- elapsed_time
    
    if (iter %% 100 == 0){print(iter/iters)}
  }
  
  
  acceptance_prob <- n_accepted/(366*iters)
  
  burnin = burnin + 1  
  result_list <- list(
    chain = chain[burnin: iters,],
    time = time,
    acceptance_prob = acceptance_prob
  )
  return (result_list)
}


results = mcmc_sampler(alpha, beta, iters, conf_level)

chain <- results$chain
time <- results$time
acceptance_prob <- results$acceptance


# plot(chain[,1])
# plot(chain[,201])
# plot.ts(y)
# plot.ts(colMeans(chain))
# plot.ts(sliced)
# plot(chain[, 366])
# plot(chain[, 367])

library(ggplot2)

prob_chain <- logit(chain[, 1:366])
mcmc_chain <- as.mcmc(prob_chain)

# credible_interval <- apply(mcmc_chain[, 1:366], 2, function(x) quantile(x, probs=c(0.025, 0.975)))
# credible_interval <- HPDinterval(mcmc_chain[, 1:366], prob = 0.95)
conf_levels = summary(mcmc_chain)$quantile
# credible_interval = confint(lm(y~x))

mean_pixt <- colMeans(mcmc_chain[, 1:366])
t_values <- 1:366 # Change this if your t values differ
# df <- data.frame(t=t_values, mean_pixt, lwr=credible_interval[1,], upr=credible_interval[2,])
df <- data.frame(t=t_values, mean_pixt, lwr=conf_levels[,1], upr=conf_levels[,5])

ggplot(df, aes(x=t)) +
  geom_ribbon(aes(ymin=lwr, ymax=upr), fill="blue", alpha=0.2) +
  geom_line(aes(y=mean_pixt), color="blue") +
  #geom_line(aes(y=y/n))
  labs(title="Posterior Mean of π(xt) with 95% Credible Intervals",
       x="t", y="π(xt)") +
  theme_minimal()


hist(chain[, 1], main="Histogram of σ²", xlab="σ²", breaks=50)
hist(chain[, 201], main="Histogram of σ²", xlab="σ²", breaks=50)
hist(chain[, 366], main="Histogram of σ²", xlab="σ²", breaks=50)
hist(chain[, 367], main="Histogram of σ²", xlab="σ²", breaks=50)

acf(chain[, 1], main="Autocorrelation of σ²")
acf(chain[, 201], main="Autocorrelation of σ²")
acf(chain[, 366], main="Autocorrelation of σ²")
acf(chain[, 367], main="Autocorrelation of σ²")
