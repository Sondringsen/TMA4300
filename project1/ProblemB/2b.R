# R code for sampling from a Gamma distribution 
# with alpha>1 and beta=1
# using the ratio of uniforms method.
# Implemented on log-scale

# Define log of fstar
logf = function(x, alpha){
    return(ifelse(0<x, (alpha-1)*log(x)-x, 0))
}

RatioOfUniformSample = function(n, alpha){

    # Keep track of which samples are accepted at the current time
    state = integer(n)
    # Store actual sample
    sample = numeric(n)
    # Store number of tries
    tries = 0

    # Define logarithm of a
    loga = (1/2)*(alpha-1)*(log(alpha-1)-1)
    # Define logarithm of b+
    logb = (1/2)*(alpha+1)*(log(alpha+1)-1)

    while (sum(state==0) != 0){
        # Calculate log of uniform(0, a)
        # Log(unif(0,a))=log(a*unif(0,1))=log(a)+log(unif(0,1))
        u1 = loga+log(runif(n, 0, 1))
        # Calculate log of uniform(0, b+)
        u2 = logb+log(runif(n, 0, 1))

        # Convert ratio of uniform inequality to log scale
        # Calculate the indices of the accepted samples, 
        # using the log scale inequality
        idx = which(2*u1<logf(exp(u2-u1), alpha))
        # Update state
        state[idx]=1
        # Insert samples
        # We need to apply exp to scale the sample
        # back from the log scale
        sample[idx]=(exp(u2-u1))[idx]

        tries=tries+1
    }   

    return (c(sample, tries))
}

n = 10000
alpha = 20

draw = RatioOfUniformSample(n, alpha)

png("ProblemB/Plots/hist2b.png")
hist(draw, breaks=50, probability=TRUE, col="lightblue",
     main="Histogram vs. Actual PDF", xlab="Value", ylab="Density")
curve(dgamma(x, alpha, 1), from=min(draw),
      to=max(draw), add=TRUE, col="darkred", lwd=2)
legend("topright", legend=c("Generated Samples", "Actual PDF"),
       col=c("lightblue", "darkred"), lwd=2, fill=c("lightblue", "darkred"))
dev.off()

PlotRatioOfUniforms = function(n){

    # Define alphas we want to try
    alphas = 2:2000
    tries = c()
    
    # Calculate tries for every alpha
    for (alpha in alphas){
        res = RatioOfUniformSample(n, alpha)
        tries = c(tries, res[n+1])
    }

    # Plot the result
    png("ProblemB/Plots/alphatries.png")
    plot(alphas, tries)
    dev.off()

}

# PlotRatioOfUniforms(1000)