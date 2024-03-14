# install.packages(
#     "INLA",
#     repos=c(getOption("repos"),
#     INLA="https://inla.r-inla-download.org/R/stable"),
#     dep=TRUE)

library("INLA")

load(file="/Users/eiriksteen/Personal/school/compstat/projects/TMA4300/project2/data/rain.rda")

ptm <- proc.time()

control_inla <- list(strategy = "simplified.laplace", int.strategy = "ccd")

formula1 = n.rain ~ -1 + f(day, model = "rw1", constr = FALSE, hyper = list(prec = list(prior = "gamma", param = c(2, 0.05))))

mod1 <- inla(formula1,
             data = rain,
             Ntrials = n.years,
             control.compute = list(config = TRUE),
             family = "binomial",
             control.inla = control_inla)

print(proc.time() - ptm)

formula2 = n.rain ~ f(day, model = "rw1", constr = FALSE, hyper = list(prec = list(prior = "gamma", param = c(2, 0.05))))

mod2 <- inla(formula2,
             data = rain,
             Ntrials = n.years,
             control.compute = list(config = TRUE),
             family = "binomial", verbose = TRUE,
             control.inla = control_inla)


mod_nocontrol <- inla(formula1,
                      data = rain,
                      Ntrials = n.years,
                      control.compute = list(config = TRUE),
                      family = "binomial")
