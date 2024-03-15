library("RTMB")

load(file="/Users/eiriksteen/Personal/school/compstat/projects/TMA4300/project2/data/rain.rda")

sigmoid <- function(x) {
  return(1 / (1 + exp(-x)))
}

f <- function(p) {

  y <- rain$n.rain
  ncr<- choose(rain$n.years, rain$n.rain)

  bin_terms <- sum(-log(ncr) - y * log(sigmoid(p$x)) - (rain$n.years - y) * log(1 - sigmoid(p$x)))
  rw_terms <- -log(1 / p$var**(1 / 2)) + sum(1 / (2 * p$var) * (p$x[-1] - p$x[-length(p$x)])**2)

  return(bin_terms + rw_terms)
}

pi <- rain$n.rain / rain$n.years + 1e-04
x <- log(pi / (1-pi))
var <- 0.05

obj <- MakeADFun(f, list(x=x, var=var), silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdreport(obj)