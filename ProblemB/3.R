GammaSample = function(n, alpha, beta){

    u = matrix(runif(n*alpha), nrow=n)
    sum = rowSums(u, dims=1)
    sample = -(1/beta)*(rowSums(log(u), dims=1))

    return(sample)
}


# n=10000
# alpha=432
# beta=123
# s = GammaSample(n, alpha, beta)
# print(mean(s))
# print(sd(s))
# true_s = rgamma(n=n, alpha, beta)
# print(mean(true_s))
# print(sd(true_s))