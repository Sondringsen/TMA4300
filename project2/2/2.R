# Import the INLA library
library("INLA")

# Load the rain data file
load(file="/Users/eiriksteen/Personal/school/compstat/projects/TMA4300/project2/data/rain.rda")

# Set the control parameters
control_inla <- list(strategy = "simplified.laplace", int.strategy = "ccd")

# Declare the formula for part a)
formula1 = n.rain ~ -1 + f(
    day, 
    model = "rw1", 
    constr = FALSE, 
    hyper = list(prec = list(prior = "loggamma", param = c(2, 0.05))))

# Initialize the timer for timing the computation
ptm <- proc.time()

# Fit the inla model on formula1
mod1 <- inla(formula1,
             data=rain, 
             Ntrials=n.years, 
             control.compute=list(config = TRUE),
             family="binomial", 
             control.inla=control_inla)

# Print the elapsed time
print("ELAPSED TIME SIMPLIFIED LAPLACE")
print(proc.time() - ptm)

plot(mod1, single = FALSE)
dev.print(pdf, file="plots/a.pdf")

# Declare the formula for part c)
formula2 = n.rain ~ f(
    day, 
    model="rw1", 
    constr=TRUE, 
    hyper = list(prec = list(prior = "loggamma", param = c(2, 0.05))))

# Fit the inla model on formula2
mod2 <- inla(formula2,
             data=rain, 
             Ntrials=n.years, 
             control.compute=list(config = TRUE),
             family="binomial", 
             control.inla=control_inla)

plot(mod2, single = FALSE)
dev.print(pdf, file="plots/c.pdf")

control_inla_gaussian <- list(strategy = "gaussian", int.strategy = "ccd")

ptm <- proc.time()

# Model with gaussian strategy
modg <- inla(formula1,
             data=rain,
             Ntrials=n.years,
             control.compute=list(config = TRUE),
             family="binomial",
             control.inla=control_inla_gaussian)

print("ELAPSED TIME GAUSSIAN")
print(proc.time() - ptm)

plot(mod2, single = FALSE)
dev.print(pdf, file="plots/gaussian-ccd.pdf")