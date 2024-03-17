# setwd("/Users/sondrerogde/Library/Mobile Documents/com~apple~CloudDocs/Documents/NTNU/8. semester/TMA4300")

# declaring parameters
alpha <- 2
beta <- 0.05
iters <- 50000
M <- 30

# loading data
load(file = "project2/data/rain.rda")
y <- rain[, "n.rain"]
n <- rain[, "n.years"]

# expit function
expit <- function(x) {1/(1+exp(-x))}

# create Q
get_Q <- function() {
  Q <- matrix(0, 366, 366)
  diag(Q) <- 2
  for (i in 1:365) {
    Q[i, i+1] <- -1 
    Q[i+1, i] <- -1 
  }
  Q[1,1] <- 1
  Q[366,366] <- 1
  return(Q)
}

mcmc_sampler <- function(alpha, beta, iters, M, burnin=ceiling(iters/10)){
  # declaring variables
  iters = iters + burnin
  n_chunks <- ceiling(366/M)
  
  # declaring initial data structures
  chain <- matrix(nrow=iters, ncol=367)
  acceptance_probs = matrix(nrow=iters-1, ncol=n_chunks)
  
  
  # calculating initial x-value
  x <- (y+1)/(n+1)
  pi <- expit(x)
  
  # sampling first sigma
  sigma_u_sq <- 1/rgamma(1, alpha, beta)
  
  # setting initial row for matrix
  chain[1,] <- c(x, sigma_u_sq)
  
  # start and end indices for all chunks
  a_idx <- seq(from = 1, to = floor(366 / M) * M + 1, by = M)
  b_idx <- a_idx+M-1
  b_idx[n_chunks] <- 366
  
  # create matrix Q
  Q <- get_Q()
  
  #lists to store precomputed matrices
  Q_AA_inv_matrices <- list()
  choleskys <- list()
  Q_AA_inv_Q_AB_matrices <- list()
  
  #precompute matrices
  for (ch in 1:n_chunks) {
    A <- a_idx[ch]:b_idx[ch]
    Q_AA <- Q[A, A]
    Q_AA_inv <- solve(Q[A,A])
    Q_AA_inv_matrices[[ch]] <- Q_AA_inv
    choleskys[[ch]] <- t(chol(Q_AA_inv))
    
    Q_AB <- Q[A, setdiff(1:366, A)]
    Q_AA_inv_Q_AB_matrices[[ch]] <- Q_AA_inv%*%Q_AB
  }
  
  # outer loop used for sampling x-vector
  for (iter in 2:iters){
    # uniform sampling used for accepting proposals
    u <- runif(n_chunks)
    
    # vectors for new samples and previous samples respectively
    new_x = c()
    old_x = chain[iter-1, ]
    
    # Gibb's sampling for sigma
    sigma_alpha <- alpha + 365/2
    sigma_beta <- beta + 1/2*sum((old_x[2:366] - old_x[1:365])^2)
    sigma_u_sq <- 1/rgamma(1, shape=sigma_alpha, scale=sigma_beta)
    
    # inner loop used for sampling x days
    for (ch in 1:n_chunks) {
      Q_AA_inv_Q_AB = Q_AA_inv_Q_AB_matrices[[ch]]
      
      # indexes able to handle all three cases mentioned in the problem
      a <- a_idx[ch]
      b <- b_idx[ch]
      x_b <- old_x[setdiff(1:366, a:b)]
      
      # calculates vectors and matrices needed for drawing conidtional multivariate normal
      x_b <- matrix(x_b, nrow = length(x_b), ncol = 1)
      mu <- -1*Q_AA_inv_Q_AB %*% x_b
      cholesky <- sqrt(sigma_u_sq)*choleskys[[ch]]
      
      # draws mulitvariate normal using cholesky decomposition
      proposal <- mu + cholesky%*%matrix(rnorm(b-a+1), nrow = b-a+1, ncol = 1)
      
      # calculating new and old pi
      new_pi <- expit(proposal)
      pi <- expit(old_x[a:b])
      
      # calculating acceptance probability (alpha)
      y_ab <- y[a:b]
      n_ab <- n[a:b]
      log_likelihood_ratio <- sum(y_ab*(log(new_pi) - log(pi)) + (n_ab-y_ab)*(log(1-new_pi)-log(1-pi)))
      acceptance_prob <- min(1, exp(log_likelihood_ratio))
      
      # accepting proposal based on acceptance probability
      if (u[ch] < acceptance_prob) {
        new_x <- c(new_x, t(proposal))
      } else {
        new_x <- c(new_x, old_x[a:b])
      }
      
      # used to track observed acceptance probabilities
      acceptance_probs[iter-1, ch] <- (u[ch] < acceptance_prob)
    }
    
    
    # updates chain and time data structures
    chain[iter, ] <- c(new_x, sigma_u_sq)
    
    if (iter %% 1000 == 0){print(iter/iters)}
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
results <- mcmc_sampler(alpha, beta, iters, M)
end_time <- proc.time()[3]
total_time = end_time - start_time
print(total_time)

# sets variables for returned values
chain <- results$chain
acceptance_probs <- results$acceptance_probs

burnin = ceiling(iters/10)

# traceplots of x_1, x_201, x_366, and sigma
plot(expit(chain[,1]), type="l", ylab="π(x_1)")
title("Traceplot of π(x_1) (including burnin)")

plot(expit(chain[,201]), type="l", ylab="π(x_201)")
title("Traceplot of π(x_201) (including burnin)")

plot(expit(chain[,366]), type="l", ylab="π(x_366)")
title("Traceplot of π(x_366) (including burnin)")

plot(chain[,367], type="l", ylab="σ²")
title("Traceplot of σ² (including burnin)")

# plot of acceptance probabilities for all days
plot(acceptance_probs, type="l")
abline(h = mean(acceptance_probs), col = "red", lty = 2)
legend("topleft", 
       legend = c("acceptance probabilities", "mean"), 
       col = c("black", "red"), 
       lty = 1, 
       bty = "n")
title("Acceptance probabilities for each batch")
print(mean(acceptance_probs[1:12]))

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

