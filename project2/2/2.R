# Import the INLA library
library("INLA")

# Load the rain data file
# load(file="/Users/eiriksteen/Personal/school/compstat/projects/TMA4300/project2/data/rain.rda")
load(file = "project2/data/rain.rda")

# Set the control parameters
control_inla <- list(strategy = "simplified.laplace", int.strategy = "ccd")

?control.inla

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


# creates a dataframe with time-values, mean x-values and confidence levels
t_values <- 1:366
df <- data.frame(t=t_values, mean=mod1$summary.fitted.values[,1], lwr=mod1$summary.fitted.values[,3], upr=mod1$summary.fitted.values[,5])

# plots the mean x-values and confidence intervals with y-values from the dataset
ggplot(df, aes(x = t)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill="Credible Intervals"), alpha = 0.2) +
  geom_line(aes(y = mean, color = "Mean probability")) +
  geom_line(aes(y = y / n, color = "Observed Ratio"), alpha=0.5) +
  labs(title = "Posterior Mean of pi(x_t) with 95% Credible Intervals",
       x = "t", y = "pi(xt)", color = "Legend") +
  scale_color_manual(values = c("Mean probability" = "blue", "Observed Ratio" = "red")) +
  scale_fill_manual(name = "Credible Intervals", values = "blue") +  # Added for ribbon fill color
  theme(legend.position = "bottom")


# b --------------------------------------------------------------------------


control_inla_gaussian <- list(strategy = "gaussian", int.strategy = "ccd")

# Model with gaussian strategy
modg <- inla(formula1,
             data=rain,
             Ntrials=n.years,
             control.compute=list(config = TRUE),
             family="binomial",
             control.inla=control_inla_gaussian)

print("ELAPSED TIME GAUSSIAN")
print(proc.time() - ptm)

# creates a dataframe with time-values, mean x-values and confidence levels
t_values <- 1:366
df <- data.frame(t=t_values, mean=modg$summary.fitted.values[,1], lwr=modg$summary.fitted.values[,3], upr=modg$summary.fitted.values[,5])

# compares the model from a and the gaussian model from b
ggplot(df, aes(x = t)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill="Credible Intervals"), alpha = 0.2) +
  geom_line(aes(y = mean, color = "Mean probability")) +
  geom_line(aes(y = y / n, color = "Observed Ratio"), alpha=0.5) +
  labs(title = "Posterior Mean of pi(x_t) with 95% Credible Intervals",
       x = "t", y = "pi(xt)", color = "Legend") +
  scale_color_manual(values = c("Mean probability" = "blue", "Observed Ratio" = "red")) +
  scale_fill_manual(name = "Credible Intervals", values = "blue") +  # Added for ribbon fill color
  theme(legend.position = "bottom")

# creates a dataframe with time-values, mean x-values and confidence levels
t_values <- 1:366
df <- data.frame(t=t_values, 
                 mean_1=mod1$summary.fitted.values[,1], lwr_1=mod1$summary.fitted.values[,3], upr_1=mod1$summary.fitted.values[,5],
                 mean_g=modg$summary.fitted.values[,1], lwr_g=modg$summary.fitted.values[,3], upr_g=modg$summary.fitted.values[,5])

# plots the mean x-values and confidence intervals with y-values from the dataset
ggplot(df, aes(x = t)) +
  geom_line(aes(y = lwr_1, color = "Lower Credible Intervals Model a"), linetype = "dashed") +
  geom_line(aes(y = upr_1, color = "Upper Credible Intervals Model a"), linetype = "dashed") +
  geom_line(aes(y = mean_1, color = "Mean probability Model a"), linetype = "dashed") +
  geom_line(aes(y = lwr_g, color = "Lower Credible Intervals Gaussian"), alpha = 0.5) +
  geom_line(aes(y = upr_g, color = "Upper Credible Intervals Gaussian"), alpha = 0.5) +
  geom_line(aes(y = mean_g, color = "Mean probability Gaussian"), alpha = 0.5) +
  #geom_line(aes(y = y / n, color = "Observed Ratio"), alpha=0.5) +
  labs(title = "Posterior Mean of pi(x_t) with 95% Credible Intervals",
       x = "t", y = "pi(xt)", color = "") +
  scale_color_manual(values = c("Mean probability Model a" = "black", "Lower Credible Intervals Model a" = "grey",
                                "Upper Credible Intervals Model a" = "grey", "Mean probability Gaussian" = "blue", 
                                "Lower Credible Intervals Gaussian" = "coral", "Upper Credible Intervals Gaussian" = "coral")) +
  theme(legend.position = "bottom")

print(mean(mod1$summary.fitted.values[,1]-modg$summary.fitted.values[,1]))








control_inla_gaussian <- list(strategy = "simplified.laplace", int.strategy = "grid")

# Model with gaussian strategy
modgrid <- inla(formula1,
             data=rain,
             Ntrials=n.years,
             control.compute=list(config = TRUE),
             family="binomial",
             control.inla=control_inla_gaussian)

print("ELAPSED TIME GAUSSIAN")
print(proc.time() - ptm)

# creates a dataframe with time-values, mean x-values and confidence levels
t_values <- 1:366
df <- data.frame(t=t_values, mean=modgrid$summary.fitted.values[,1], lwr=modgrid$summary.fitted.values[,3], upr=modgrid$summary.fitted.values[,5])

# compares the model from a and the grid strategy from b
ggplot(df, aes(x = t)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill="Credible Intervals"), alpha = 0.2) +
  geom_line(aes(y = mean, color = "Mean probability")) +
  geom_line(aes(y = y / n, color = "Observed Ratio"), alpha=0.5) +
  labs(title = "Posterior Mean of pi(x_t) with 95% Credible Intervals",
       x = "t", y = "pi(xt)", color = "Legend") +
  scale_color_manual(values = c("Mean probability" = "blue", "Observed Ratio" = "red")) +
  scale_fill_manual(name = "Credible Intervals", values = "blue") +  # Added for ribbon fill color
  theme(legend.position = "bottom")

# creates a dataframe with time-values, mean x-values and confidence levels
t_values <- 1:366
df <- data.frame(t=t_values, 
                 mean_1=mod1$summary.fitted.values[,1], lwr_1=mod1$summary.fitted.values[,3], upr_1=mod1$summary.fitted.values[,5],
                 mean_g=modgrid$summary.fitted.values[,1], lwr_g=modgrid$summary.fitted.values[,3], upr_g=modgrid$summary.fitted.values[,5])

# plots the mean x-values and confidence intervals with y-values from the dataset
ggplot(df, aes(x = t)) +
  geom_line(aes(y = lwr_1, color = "Lower Credible Intervals Model a"), linetype = "dashed") +
  geom_line(aes(y = upr_1, color = "Upper Credible Intervals Model a"), linetype = "dashed") +
  geom_line(aes(y = mean_1, color = "Mean probability Model a"), linetype = "dashed") +
  geom_line(aes(y = lwr_g, color = "Lower Credible Intervals Gaussian"), alpha = 0.5) +
  geom_line(aes(y = upr_g, color = "Upper Credible Intervals Gaussian"), alpha = 0.5) +
  geom_line(aes(y = mean_g, color = "Mean probability Gaussian"), alpha = 0.5) +
  #geom_line(aes(y = y / n, color = "Observed Ratio"), alpha=0.5) +
  labs(title = "Posterior Mean of pi(x_t) with 95% Credible Intervals",
       x = "t", y = "pi(xt)", color = "") +
  scale_color_manual(values = c("Mean probability Model a" = "black", "Lower Credible Intervals Model a" = "grey",
                                "Upper Credible Intervals Model a" = "grey", "Mean probability Gaussian" = "blue", 
                                "Lower Credible Intervals Gaussian" = "coral", "Upper Credible Intervals Gaussian" = "coral")) +
  theme(legend.position = "bottom")

print(mean(mod1$summary.fitted.values[,1]-modgrid$summary.fitted.values[,1]))




# c --------------------------------------------------------------------------




# Declare the formula for part c)
formula2 = n.rain ~ f(
    day, 
    model = "rw1", 
    constr = TRUE, 
    hyper = list(prec = list(prior = "loggamma", param = c(2, 0.05))))

# Fit the inla model on formula2
mod2 <- inla(formula2,
             data=rain, 
             Ntrials=n.years, 
             control.compute=list(config = TRUE),
             family="binomial", 
             control.inla=control_inla)

# plot(mod2, single = FALSE)
# dev.print(pdf, file="plots/c.pdf")

ptm <- proc.time()

# creates a dataframe with time-values, mean x-values and confidence levels
t_values <- 1:366
df <- data.frame(t=t_values, mean=mod2$summary.fitted.values[,1], lwr=mod2$summary.fitted.values[,3], upr=mod2$summary.fitted.values[,5])

# plots the mean x-values and confidence intervals with y-values from the dataset
ggplot(df, aes(x = t)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill="Credible Intervals"), alpha = 0.2) +
  geom_line(aes(y = mean, color = "Mean probability")) +
  geom_line(aes(y = y / n, color = "Observed Ratio"), alpha=0.5) +
  labs(title = "Posterior Mean of pi(x_t) with 95% Credible Intervals",
       x = "t", y = "pi(xt)", color = "Legend") +
  scale_color_manual(values = c("Mean probability" = "blue", "Observed Ratio" = "red")) +
  scale_fill_manual(name = "Credible Intervals", values = "blue") +  # Added for ribbon fill color
  theme(legend.position = "bottom")



