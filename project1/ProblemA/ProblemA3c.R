# A function of two inputs alpha and n, 
# alpha is the first input parameter and n is the number of draws
f = function(alpha, n){
  # u is a vector of n i.i.d. random uniforms
  u = runif(n)
  
  # Uses inverse transform to get the correct distribution
  return (1/alpha*log(u/(1-u)))
}

# Invokes the function
draw = f(4, 1000)

# Check that the function has the correct properties
print(mean(draw))
print(var(draw))

# Plots the histogram
hist(draw)
