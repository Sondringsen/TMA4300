# Import the INLA library
library("INLA")

# Load the rain data file
load(file="/Users/eiriksteen/Personal/school/compstat/projects/TMA4300/project2/data/rain.rda")

# Initialize the timer for timing the computation
ptm <- proc.time()

# Set the control parameters
control_inla <- list(strategy = "simplified.laplace", int.strategy = "ccd")

# Declare the formula for part a)
formula1 = n.rain ~ -1 + f(day, model = "rw1", constr = FALSE, hyper = list(prec = list(prior = "loggamma", param = c(2, 0.05))))

# Fit the inla model on formula1
mod1 <- inla(formula1,
             data=rain, 
             Ntrials=n.years, 
             control.compute=list(config = TRUE),
             family="binomial", 
             control.inla=control_inla)

# Print the elapsed time
print("ELAPSED TIME")
print(proc.time() - ptm)

formula2 = n.rain ~ f(day, model="rw1", constr=TRUE, hyper = list(prec = list(prior = "loggamma", param = c(2, 0.05))))


mod2 <- inla(formula2,
             data=rain, 
             Ntrials=n.years, 
             control.compute=list(config = TRUE),
             family="binomial", 
             control.inla=control_inla)

# # Fit the inla model on formula1 without the control parameters,
# # to assess the robustness
# mod_nocontrol <- inla(formula1,
#                       data = rain,
#                       Ntrials = n.years,
#                       control.compute = list(config = TRUE),
#                       family = "binomial")
