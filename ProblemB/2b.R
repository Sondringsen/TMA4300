logf = function(x, alpha){
    return(ifelse(0<x, (alpha-1)*log(x)-x, 0))
}

RatioOfUniformSample = function(n, alpha){

    state = integer(n)
    sample = numeric(n)
    tries = 0

    loga = (1/2)*(alpha-1)*(log(alpha-1)-1)
    logb = (1/2)*(alpha+1)*(log(alpha+1)-1)

    while (sum(state==0) != 0){
        u1 = loga+log(runif(n, 0, 1))
        u2 = logb+log(runif(n, 0, 1))

        idx = which(2*u1<logf(exp(u2-u1), alpha))
        state[idx]=1
        sample[idx]=(exp(u2-u1))[idx]

        tries=tries+1
    }   

    return (c(sample, tries))
}

# n=1000
# alpha=2000
# res = RatioOfUniformSample(n, alpha)
# sample = res[1:n]
# tries = res[n+1]

# print("MEAN AND SD OF GAMMA SAMPLE")
# print(mean(sample))
# print(sd(sample))
# print("TRUE MEAN AND SD")
# true_s = rgamma(n=n, alpha, 1.0)
# print(mean(true_s))
# print(sd(true_s))
# print("NUMBER OF TRIES")
# print(tries)

PlotRatioOfUniforms = function(n){

    alphas = 2:2000
    tries = c()
    
    for (alpha in alphas){
        res = RatioOfUniformSample(n, alpha)
        tries = c(tries, res[n+1])
    }

    png("ProblemB/Plots/alphatries.png")
    plot(alphas, tries)
    dev.off()

}

PlotRatioOfUniforms(1000)