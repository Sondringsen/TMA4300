# A function of two inputs lambda and n, 
# lambda is the rate and n is the number of draws
exponential = function(lambda, n){
  # u is a vector of n i.i.d. random uniforms
  u = runif(n)
  
  # Uses inverse transform to get an exponential distribution
  return (-log(1 - u)/lambda)
}

# Invokes the function
exp_draw = exponential(4, 1000)

# Check that the function has the correct properties.
print(mean(exp_draw))
print(var(exp_draw))

# Plots the histogram
hist(exp_draw)
