y <- rain[, "n.rain"]
n <- rain[, "n.years"]

plot.ts(y, main="Response as function of time", xlab="t")
plot.ts(y/n, main="y/n as function of time", xlab="t")
acf(y/n, main="ACF of y/n", xlab="lag")
