
GDensity = function(x, alpha){
    c = alpha / (1 + exp(1)/alpha)
    x = ifelse (x <= 0, 0, x)
    x = ifelse(x > 0 & x < 1, c*x**(alpha-1), x)
    x = ifelse(x >= 1, c*exp(-x), x)

    return(x)
}

# A function of two inputs alpha and n, 
# alpha is the first input parameter and n is the number of draws
GSample = function(alpha, n){
  # u is a vector of n i.i.d. random uniforms
  u = runif(n)
  
  # Defines the normalising constant
#   c = alpha/(1-alpha*exp(-1))
  c = (alpha*exp(1))/(exp(1)+alpha)

  # Uses inverse transform to get the distribution
  x = ifelse (u < c/alpha, 
          (u*alpha/c)**(1/alpha),
          (-log(-u/c+1/alpha+exp(-1)))
  )
  return (x)
}

FDensity = function(x, alpha){
    return(1/gamma(alpha))*x**(alpha-1)*exp(-x)
}

RejectionSample = function(n, alpha_g, alpha_f){

    state = integer(n)
    sample = numeric(n)

    while (sum(state==0) != 0){
        x = GSample(alpha_g, n)
        f = FDensity(x, alpha_f)
        g = GDensity(x, alpha_g)
        c = (exp(1) + alpha_g) / (alpha_g*exp(1)*gamma(alpha_f))
        # c = (1 - alpha_g*exp(-1))/(alpha_g*gamma(alpha_f))
        alpha = (1/c)*f/g
        u = runif(n)
        idx = which(alpha - u > 0)
        state[idx] = 1
        sample[idx] = x[idx]
    }

    return (sample)

}

n = 100000
alpha_f = 0.5
alpha_g = 0.5

rs = RejectionSample(n, alpha_g, alpha_f)
true_s = rgamma(n=n, alpha_f, 1.0)

print(mean(rs))
print(mean(true_s))
print("@@@")
print(sd(rs))
print(sd(true_s))