# A function of two inputs alpha and n, 
# alpha is the first input parameter and n is the number of draws
f = function(alpha, n){
  # u is a vector of n i.i.d. random uniforms
  u = runif(n)
  
  # Defines the normalising constant
  c = alpha/(1-alpha*exp(-1))
    
  # Uses inverse transform to get the distribution
  x = ifelse (u < c/alpha, 
          (u*alpha/c)**(1/alpha),
          (-log(-u/c+1/alpha+exp(-1)))
  )
  return (x)
}

# Invokes the function
draw = f(.5, 1000)

# Check that the function has the correct properties
print(mean(draw))
print(var(draw))

# Plots the histogram
hist(draw)
