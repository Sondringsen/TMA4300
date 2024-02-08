
f = function(x, alpha){
    return(ifelse(0<x, x**(alpha-1)*exp(-x), 0))
}

RatioOfUniformSample = function(n, alpha){

    state = integer(n)
    sample = numeric(n)

    a = ((alpha-1)**(alpha-1) * exp(-(alpha-1)))**(1/2)
    b = ((alpha+1)**(alpha+1) * exp(-(alpha+1)))**(1/2)
    b_ = -b

    while (sum(state==0) != 0){
        u1 = runif(n, 0, a)
        u2 = runif(n, b_, b)

        idx = which(u1**2 < f(u2/u1, alpha))
        state[idx]=1
        sample[idx]=(u2/u1)[idx]
    }   
    return(sample)
}

sample = RatioOfUniformSample(1000, 10)

print(mean(sample))
print(sd(sample))