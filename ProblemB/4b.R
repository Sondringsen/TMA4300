source("ProblemB/2b.R")

BetaSample = function(n, alpha, beta){

    x = RatioOfUniformSample(n, alpha)
    y = RatioOfUniformSample(n, beta)
    z = x/(x+y)

    return(z)
}

sample = BetaSample(1000,2,2)

print("BETA")
print(mean(sample))
print(sd(sample))