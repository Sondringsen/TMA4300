# R code for sampling from the Gamma distribution for any alpha>0 and beta>0

# Import Gamma sampling function for 0<alpha<1
source("ProblemB/1b.R")
# Import Gamma sampling function for alpha>1
source("ProblemB/2b.R")

# Function to sample from Gamma with integer alphas
GammaSampleIntegerAlpha = function(n, alpha, beta){

    # Create matrix of uniform samples of shape (n, alpha)
    u = matrix(runif(n*alpha), nrow=n)
    # Calculate sample using the sum of
    # exponentials formula
    sample = -(1/beta)*(rowSums(log(u), dims=1))

    return(sample)
}

GammaSample = function(n, alpha, beta){

    # if 0<alpha<1, use rejection sampling from part 1
    if (alpha>0&alpha<1){
        return((1/beta)*RejectionSample(n, alpha))
    }
    # if alpha>1, use ratio of uniform from part 2
    else if (alpha>1) {
       return((1/beta)*RatioOfUniformSample(n, alpha))
    }
    # for alpha=1 use the integer alpha sampler
    else if (alpha==1) {
       return(GammaSampleIntegerAlpha(n, alpha, beta))
    }
}

fwithbeta = function(x, alpha, beta){
    # Define Gamma density
    return(ifelse(0<x, (beta**alpha/gamma(alpha))*x**(alpha-1)*exp(-beta*x), 0))
}

n = 10000
alpha = 120
beta = 230

draw = GammaSample(n, alpha, beta)

png("ProblemB/Plots/hist3.png")
hist(draw, breaks=50, probability=TRUE, col="lightblue",
     main="Histogram vs. Actual PDF", xlab="Value", ylab="Density")
curve(fwithbeta(x, alpha, beta), from=min(draw),
      to=max(draw), add=TRUE, col="darkred", lwd=2)
legend("topright", legend=c("Generated Samples", "Actual PDF"),
       col=c("lightblue", "darkred"), lwd=2, fill=c("lightblue", "darkred"))
dev.off()