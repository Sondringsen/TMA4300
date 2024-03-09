# R code for rejection sampling from the Gamma distribution
# (with alpha between 0 and 1, and beta=1)

g = function(x, alpha){

    # Define normalizing constant
    c = 1/((1/alpha) + exp(-1))

    # Calculate indices corresponding to conditions
    idx0 = which(x<=0)
    idx1 = which(x>0&x<1)
    idx2 = which(x>=1)

    # Calculate conditional densities
    x[idx0] = 0
    x[idx1] = c*x[idx1]**(alpha-1)
    x[idx2] = c*exp(-x[idx2])

    return(x)
}

gSample = function(alpha, n){

    # Define normalizing constant
    c = 1/((1/alpha) + exp(-1))
    
    # u is a vector of n i.i.d. random uniforms
    u = runif(n)

    # Uses inverse transform to get the distribution
    x = ifelse (u < c/alpha, 
          (u*alpha/c)**(1/alpha),
          -log(-u/c+1/alpha+exp(-1))
    )
    return (x)
}

f = function(x, alpha){
    # Define Gamma density
    return(ifelse(0<x, (1/gamma(alpha))*x**(alpha-1)*exp(-x), 0))
}

RejectionSample = function(n, alpha){

    # Keep track of which samples are accepted at the current time
    state = integer(n)
    # Store actual sample
    sample = numeric(n)

    while (sum(state==0) != 0){
        
        # Sample from the proposal distribution
        x = gSample(alpha, n)
        # Define the normalizing constant
        c = (alpha**(-1) + exp(-1)) / gamma(alpha)
        # Calculate the ratio f(x)/(c*g(x)) 
        # normally referred to as alpha
        lambda = (1/c)*f(x, alpha)/g(x, alpha)
        u = runif(n)
        # Calculate indices of accepted samples
        idx = which(lambda - u > 0)
        # Update state of accepted samples
        state[idx] = 1
        # Insert accepted samples
        sample[idx] = x[idx]
    }

    return (sample)

}

n = 10000
alpha = 0.9

draw = RejectionSample(n, alpha)

png("ProblemB/Plots/hist1b.png")
hist(draw, breaks=50, probability=TRUE, col="lightblue",
     main="Histogram vs. Actual PDF", xlab="Value", ylab="Density")
curve(f(x, alpha), from=min(draw),
      to=max(draw), add=TRUE, col="darkred", lwd=2)
legend("topright", legend=c("Generated Samples", "Actual PDF"),
       col=c("lightblue", "darkred"), lwd=2, fill=c("lightblue", "darkred"))
dev.off()