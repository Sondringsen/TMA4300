# setwd("/Users/sondrerogde/Library/Mobile Documents/com~apple~CloudDocs/Documents/NTNU/8. semester/TMA4300")


alpha <- 2
beta <- 0.05
iters <- 5000
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
    
    
    new_x = c()
    old_x = chain[iter-1, ]
    
    for (t in 1:366) {
      if (t == 1){
        proposal <- rnorm(1, x[2], sqrt(sigma_u_sq))
      } else if (t != 366) {
        # print(tail(new_x, n=1))
        proposal <- rnorm(1, (tail(new_x, n=1) + old_x[t+1])/2, sqrt(sigma_u_sq))
      } else if (t == 366) {
        proposal <- rnorm(1, tail(new_x, n=1), sqrt(sigma_u_sq))
      }
      
      new_pi = logit(proposal)
      pi <- logit(old_x[t])
      acceptance_prob <- min(1, dbinom(y, n, new_pi)/dbinom(y, n, pi))
      new_x <- c(new_x, ifelse(u[t] < acceptance_prob, proposal, old_x[t]))
      n_accepted <- n_accepted + (u[t] < acceptance_prob)
    }
    
    
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
}


results = mcmc_sampler(alpha, beta, iters, conf_level)

chain <- results$chain
time <- results$time
acceptance_prob <- results$acceptance

print(chain)

sliced = colMeans(chain)[10:365]

plot(chain[,201])
plot.ts(y)
plot.ts(colMeans(chain))
plot.ts(sliced)
plot(chain[, 366])
plot(chain[, 367])
