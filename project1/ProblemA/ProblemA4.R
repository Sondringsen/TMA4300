# Function using the Box Muller algorithm to draw independent
# samples from a normal distribution.
# Takes one argument n specifying the number of samples to draw
BoxMuller = function(n){
  # Draws n samples from two i.i.d. uniform distribution
  u1 = runif(n)
  u2 = runif(n)
  
  # Returns n samples from a standard normal distribution
  return ( sqrt(-2*log(u1)) * cos(2*pi*u2) )
}

# Invokes the function
Z = BoxMuller(1000)

# Check that the function has the correct properties
print(mean(Z))
print(var(Z))

# Plots the histogram
hist(Z)