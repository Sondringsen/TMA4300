# MCMC 1g --------------------------------------------------------------------------

#setwd("/Users/sondrerogde/Library/Mobile Documents/com~apple~CloudDocs/Documents/NTNU/8. semester/TMA4300")

# declaring parameters
alpha <- 2
beta <- 0.05
iters <- 50000
M <- 30

# loading data
load(file = "project2/data/rain.rda")
y <- rain[, "n.rain"]
n <- rain[, "n.years"]

# logit function
logit <- function(x) {1/(1+exp(-x))}

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
  pi <- logit(x)
  
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
      new_pi <- logit(proposal)
      pi <- logit(old_x[a:b])
      
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
results <- mcmc_sampler(alpha, beta, iters, M)

# sets variables for returned values
chain <- results$chain

# calculates probabilities
prob_chain <- logit(chain[, 1:366])

# creates a mcmc-object of the chain (excluding sigma)
mcmc_chain <- as.mcmc(prob_chain)

# calculates mean values of each day
mean_mcmc <- colMeans(mcmc_chain[, 1:366])

# INLA --------------------------------------------------------------------------

# Import the INLA library
library("INLA")

# Set the control parameters
control_inla <- list(strategy = "simplified.laplace", int.strategy = "ccd")

# Declare the formula for part a)
formula1 = n.rain ~ -1 + f(
  day, 
  model = "rw1", 
  constr = FALSE, 
  hyper = list(prec = list(prior = "loggamma", param = c(2, 0.05))))


# Fit the inla model on formula1
mod1 <- inla(formula1,
             data=rain, 
             Ntrials=n.years, 
             control.compute=list(config = TRUE),
             family="binomial", 
             control.inla=control_inla)


mean_inla = mod1$summary.fitted.values[,1]


# RTMB --------------------------------------------------------------------------



# Import RTMB
library("RTMB")

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

mle_rtmb <- logit(parameters$par.fixed[1:366])


# Plotting --------------------------------------------------------------------------



t_values <- 1:366
df <- data.frame(t = t_values, mean_mcmc=mean_mcmc, mean_inla=mean_inla, mle_rtmb=mle_rtmb)

# compares the model from 1g and the RTMB
ggplot(df, aes(x = t)) +
  geom_line(aes(y = mean_mcmc, color = "Mean MCMC")) +
  geom_line(aes(y = mean_inla, color = "Mean INLA")) +
  geom_line(aes(y = mean_rtmb, color = "MLE RTMB")) +
  labs(title = "Posterior Mean/MLE of the three models",
       x = "t", y = "π(xt)", color = "Legend") +
  scale_color_manual(values = c("Mean MCMC" = "blue", "Mean INLA" = "green",
                                "MLE RTMB" = "orange")
                     ) +
  theme(legend.position = "bottom")
