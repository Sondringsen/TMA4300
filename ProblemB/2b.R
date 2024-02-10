
f = function(x, alpha){
    return(ifelse(0<x, x**(alpha-1)*exp(-x), 0))
}

RatioOfUniformSample = function(n, alpha){

    state = integer(n)
    sample = numeric(n)
    tries = 0

    a = ((alpha-1)**(alpha-1) * exp(-(alpha-1)))**(1/2)
    b = ((alpha+1)**(alpha+1) * exp(-(alpha+1)))**(1/2)
    b_ = -b
    
    while (sum(state==0) != 0){
        u1 = runif(n, 0, a)
        u2 = runif(n, b_, b)

        idx = which(u1**2 < f(u2/u1, alpha))
        state[idx]=1
        sample[idx]=(u2/u1)[idx]

        tries=tries+1
    }   

    return (c(sample, tries))
}

res = RatioOfUniformSample(1000, 100)

sample = res[1:1000]
tries = res[1001]

print("MEAN AND SD OF GAMMA SAMPLE")
print(mean(sample))
print(sd(sample))
print("NUMBER OF TRIES")
print(tries)

PlotRatioOfUniforms = function(n){

    alphas = 1:2000

    res

    for (alpha in alphas){

    }
}