exponential = function(lambda, n){
  u = runif(n)
  return (-log(1 - u, base=exp(1))/lambda)
}

exp_draw = exponential(4, 1000)
print(mean(exp_draw))
print(var(exp_draw))
hist(exp_draw)
