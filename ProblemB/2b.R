
f = function(x, alpha){
    return(ifelse(0<x, x**(alpha-1)*exp(-x), 0))
}

# RatioOfUniformSample = function(n, alpha){

#     state = integer(n)
#     sample = numeric(n)
#     tries = 0

#     a = ((alpha-1)**(alpha-1) * exp(-(alpha-1)))**(1/2)
#     b = ((alpha+1)**(alpha+1) * exp(-(alpha+1)))**(1/2)
#     b_ = -b
    
#     while (sum(state==0) != 0){
#         u1 = runif(n, 0, a)
#         u2 = runif(n, b_, b)

#         idx = which(u1**2 < f(u2/u1, alpha))
#         state[idx]=1
#         sample[idx]=(u2/u1)[idx]
#         tries=tries+1
#     }   

#     return (c(sample, tries))
# }

RatioOfUniformSample = function(n, alpha){

    state = integer(n)
    sample = numeric(n)
    tries = 0

    a = log(((alpha-1)**(alpha-1) * exp(-(alpha-1)))**(1/2))
    b = log(((alpha+1)**(alpha+1) * exp(-(alpha+1)))**(1/2))
    b_ = -b

    print("A")
    print(a)
    print("B")
    print(b)
    print("B-")
    print(b_)
    
    while (sum(state==0) != 0){
        u1 = runif(n, 0, a)
        u2 = runif(n, b_, b)

        idx = which(2*u1 <= log(f(exp(u2)/exp(u1), alpha)))
        # idx = which(2*u1 < log(f(u2/u1, alpha)))
        state[idx]=1
        sample[idx]=(exp(u2)/exp(u1))[idx]

        tries=tries+1
    }   

    return (c(sample, tries))
}

n=10000
alpha=10
res = RatioOfUniformSample(n, alpha)
sample = res[1:n]
tries = res[n+1]

print("MEAN AND SD OF GAMMA SAMPLE")
print(mean(sample))
print(sd(sample))
print("TRUE MEAN AND SD")
true_s = rgamma(n=n, alpha, 1.0)
print(mean(true_s))
print(sd(true_s))
print("NUMBER OF TRIES")
print(tries)

PlotRatioOfUniforms = function(n){

    alphas = 1:2000

    res

    for (alpha in alphas){

    }
}