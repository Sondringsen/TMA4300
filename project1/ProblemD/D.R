#Global
SHOW_HISTOGRAMS <- TRUE

# Define the unnormalized target distribution function 
f_star <- function(x) {
  return((x+2)^y1 * (1-x)^(y2+y3) * x^y4)
}

# Rejection sampling function to simulate samples from the target distribution
rejection_sampling <- function(n,C) {
  samples <- numeric(n)  
  count <- 0  
  total <- 0
  
  while(count < n) {
    total <- total + 1
    x <- runif(1)  
    u <- runif(1)  
    
    # Acceptance criterion
    if(f_star(x) / (C) > 1){print("C chosen to be too small")}
    if(u <= f_star(x) / (C)) { #proposal density is just 1
      count <- count + 1
      samples[count] <- x
    }
  }
  
  print("Random numbers generated per sample: ")
  print(total/count)
  return(samples)
}

Monte_Carlo_integration <- function(n,C){
  samples <- rejection_sampling(n,C)
  mu = mean(samples)
  
  if(SHOW_HISTOGRAMS){
    normalizing_constant = integrate(f_star, lower=0, upper=1)$value
    pdf <- function(x){return(f_star(x)/normalizing_constant)}
    hist(samples, breaks=50,  freq=FALSE, col="lightblue", main="Histogram of sampled values", xlab="Value", ylab="Density")
    curve(pdf, col="darkred", add=TRUE, lwd=2)
    abline(v=mu, col="black", lwd=2, lty=2)
    legend("topright", legend=c("Sample density"), fill=c("lightblue"))
  }
  
  print("Estimated mean:")
  print(mu)
  return(mu)
}

Numerical_integration <- function(){
  normalizing_constant = integrate(f_star, lower=0, upper=1)$value
  integrand <- function(x){return(x*f_star(x)/normalizing_constant)}
  mu = integrate(integrand, lower=0, upper=1)
  print("Estimated mean: ")
  print(mu)
  return(mu)
}

Monte_Carlo_integration_beta15prior <- function(n,C){
  samples <- rejection_sampling(n,C)
  weights <- (1-samples)^4
  
  mu = mean(samples*weights)/mean(weights)
  
  print("Estimated mean with Beta(1,5) prior:")
  print(mu)
  return(mu)
}

Numerical_integration_beta15prior <- function(){
  density <- function(x){f_star(x)*(1-x)^4}
  normalizing_constant = integrate(density, lower=0, upper=1)$value
  integrand <- function(x){return(x*density(x)/normalizing_constant)}
  mu = integrate(integrand, lower=0, upper=1)
  print("Estimated mean with Beta(1,5) prior: ")
  print(mu)
  return(mu)
}


#Estimating a usable value for C
{
y1 = 125
y2 = 18
y3 = 20
y4 = 34

a = y1 + y2 + y3 + y4
b = -y1 + 2*y2 + 2*y3 + y4
c = -2*y4
discriminant = b^2 - 4*a*c

root1 <- (-b + sqrt(discriminant)) / (2*a)
root2 <- (-b - sqrt(discriminant)) / (2*a)
roots = c(root1, root2)
validRoots <- c(roots[roots >= 0 & roots <= 1]) #Should always have at least one
}
C = 1.1 * max(sapply(validRoots, f_star)) #multiplying by 1.1 to ensure c > f_star even if rounding has occurred

Monte_Carlo_integration(10000,C)
Numerical_integration()

theoretical_alpha = integrate(f_star, lower=0, upper=1)$value/C
print("Theoretical random numbers per sample (numerically estimated):")
print(1/theoretical_alpha)

Monte_Carlo_integration_beta15prior(10000,C)
Numerical_integration_beta15prior()

