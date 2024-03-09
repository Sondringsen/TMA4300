# Code for sampling from a Beta distribution, given that we can sample from a Gamma distribution
# Gamma sampler imported from Problem B2b)

source("ProblemB/3.R")

BetaSample = function(n, alpha, beta){

    # Sample x from Gamma
    x = GammaSample(n, alpha, 1)
    # Sample y from Gamma
    y = GammaSample(n, beta, 1)
    # Calculate beta sample
    z = x/(x+y)

    return(z)
}

n=10000
alpha=100
beta=14
draw = BetaSample(n,alpha,beta)

png("ProblemB/Plots/hist4b.png")
hist(draw, breaks=50, probability=TRUE, col="lightblue",
     main="Histogram vs. Actual PDF", xlab="Value", ylab="Density")
curve(dbeta(x, alpha, beta), from=min(draw),
      to=max(draw), add=TRUE, col="darkred", lwd=2)
legend("topright", legend=c("Generated Samples", "Actual PDF"),
       col=c("lightblue", "darkred"), lwd=2, fill=c("lightblue", "darkred"))
dev.off()