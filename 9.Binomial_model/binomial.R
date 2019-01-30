##########################################
# Code for creating your first model.    #
# Based on analysis of binomial data:    #
# a vector of 0s and 1s indicating       #
# whether or not a female gave birth     #
# within a particular 5-year period.     #
##########################################


#----------------------------#
# Load appropriate libraries #
#----------------------------#
library(rstan)
library(ggplot2)


#----------------------------#
#  Read the data into R      #
#----------------------------#
calving = read.table("data/calving.csv", header = FALSE, sep = ",")
str(calving)

#----------------------------#
#  Get a feel for the data   #
#----------------------------#
females = length(calving$V1)
females

mothers = sum(calving$V1)
mothers

mothers / females

#--- Organize and Plot the data ---#
types = c("reproductive", "not-reproductive")
counts = c(mothers, females - mothers)
data = data.frame(types, counts)
data

ggplot(data) + 
  theme_bw() +
  geom_point(aes(x = counts, y = types)) +
  labs(x = "", y = "")

# Adjust axis limits
ggplot(data) + 
  theme_bw() +
  geom_point(aes(x = counts, y = types)) + 
  labs(x = "", y = "") +
  xlim(50, 175)
  

##################################
#    Frequentist Approach        #
##################################
binom.test(x = mothers, n = females, p = 0.83)


##################################
# Bayesian Approach (Homebrew)   #
##################################

#----------------------#
#    By Hand           #
#----------------------#
#--- Propose a value for the probability of success ---#
proposal = runif(n = 1, min = 0, max = 1)
proposal

#--- Generate a random binomial data set based on this value ---#
simulated = rbinom(n = 218, size = 1, prob = proposal)

#--- See if results in same number of successes as out data ---#
sum(simulated)

# Repeat above steps many times #


#----------------------#
# Using binomialSim.R  #
#----------------------#
source("binomialSim.R")
binomialCalc(nTrials = 218, nSuccesses = 79, nSteps = 10000)  # Change nSteps to play around


##################################
#  Bayesian Approach (Stan)      #
##################################

#------------------------------#
# Save data as a list for STAN #
#------------------------------#
dataList = list (
  y = mothers,
  N = females
)

#------------------------------#
#       Define the Model       #
#------------------------------#
modelString = "
  data {
    int<lower=0> N;
    int<lower=0> y;
  }

  parameters {
    real theta;
  }

  model {
    // Likelihood
    y ~ binomial(N, theta);
    
    // Priors
    theta ~ uniform(0, 1);
  }
"
writeLines(modelString, con="model.stan")


#--------------------------#
#     run STAN             #
#--------------------------#
stanFit = stan(file = "model.stan", 
                data = dataList, 
                pars = "theta",
                warmup = 2000,
                iter = 5000, 
                chains = 3)


#-------------------------#
# Check MCMC Performance  #
#-------------------------#
print(stanFit)
print(stanFit, probs = c(0.055, 0.5, 0.945))
stan_trace(stanFit, pars = "theta", inc_warmup = TRUE)
stan_plot(stanFit, par = "theta")
stan_plot(stanFit, par = "theta", ci_level = 0.89, outer_level = 0.89)
plot(stanFit, par = "theta", show_density = TRUE, ci_level = 0.89, fill_color = "dodgerblue4")


#-----------------------#
#    Extract Data       #
#-----------------------#
posteriors = as.data.frame(stanFit)
str(posteriors)
head(posteriors)

postMean = mean(posteriors$theta)
limits = quantile(posteriors$theta, probs = c(0.055, 0.945))
limits
postHDI = as.numeric(limits)
postHDI

#-----------------------#
#    Plot the Data      #
#-----------------------#

#--- Histogram ---#
ggplot(posteriors, aes(x = theta)) +
  theme_bw() +
  geom_histogram(fill = "dodgerblue4", alpha = 0.6) +
  geom_point(aes(x = postMean, y = 0), size = 3) +
  geom_errorbarh(aes(y = 0, xmin = postHDI[1], xmax = postHDI[2]), size = 1)

#--- Density Plot #1---#
ggplot(posteriors, aes(x = theta)) +
  theme_bw() +
  geom_density(fill = "dodgerblue4", alpha = 0.6) +
  geom_point(aes(x = postMean, y = 0), size = 3) +
  geom_errorbarh(aes(y = 0, xmin = postHDI[1], xmax = postHDI[2]), size = 1)

#--- Density Plot #2---#
ggplot(posteriors, aes(x = theta)) +
  theme_bw() +
  geom_density(fill = "dodgerblue4", alpha = 0.6) +
  geom_point(aes(x = postMean, y = 0), size = 3) +
  geom_errorbarh(aes(y = 0, xmin = postHDI[1], xmax = postHDI[2]), size = 1, height = 0)

#--- Point Plot ---#
data1 = c(postMean, postHDI)
data = data.frame(c("theta", data1))
data

ggplot(data, aes(y = data[1, 1])) +
  theme_bw() +
  geom_point(aes(x = data[2, 1]), size = 4, colour = "dodgerblue4") +
  geom_errorbarh(aes(xmin = data[3, 1], xmax = data[4, 1]), size = 2, colour = "dodgerblue4", alpha = 0.5, height = 0) +
  labs(y = "", x = "")


#--------------------------#
# Check validity of priors #
#--------------------------#

#--- Generate prior values ---#
priors = runif(n = 9000, min = 0, max = 1)

#--- Combine with posteriors ---#
labels = c(rep("prior", times = 9000), rep("posterior", times = 9000))
values = c(priors, posteriors$theta)
priorCheck = data.frame(labels, values)
head(priorCheck)

#--- Plot the data ---#
ggplot(priorCheck) +
  theme_bw() +
  geom_density(aes(x = values, group = labels, fill = labels, color = labels), alpha = 0.6) +
  scale_fill_manual(values = c("dodgerblue4", "dodgerblue")) +
  scale_color_manual(values = c("dodgerblue4", "dodgerblue")) +
  xlab("theta")



#----------------------------#
# Posterior Predictive Check #
#----------------------------#

#--- Re-write model to generate random values ---#
modelString = "
  data {
    int<lower=0> N;
    int<lower=0> y;
  }

  parameters {
    real theta;
  }

  model {
    // Priors
    theta ~ uniform(0, 1);

    // Likelihood
    y ~ binomial(N, theta);
  }

  generated quantities {
    real y_pred;
    y_pred = binomial_rng(N, theta);
  }
"
writeLines(modelString, con="model.stan")


#--------------------------#
#     run STAN             #
#--------------------------#
stanFit <- stan(file = "model.stan", 
                data = dataList, 
                pars = c("theta", "y_pred"),
                warmup = 2000,
                iter = 5000, 
                chains = 3)


#-----------------------#
#    Extract Data       #
#-----------------------#
posteriors = as.data.frame(stanFit)
postMean = mean(posteriors$y_pred)
limits = quantile(posteriors$y_pred, probs = c(0.055, 0.945))
postHDI = as.numeric(limits)

#--- Point Plot ---#
data1 = data.frame(postMean, postHDI)
data2 = rep("reproductive", times = 2)
data3 = rep(mothers, times = 2)
data = data.frame(data2, data1, data3)
data

ggplot(data, aes(y = data[, 1])) +
  theme_bw() +
  geom_point(aes(x = data[, 2]), size = 5, colour = "dodgerblue4") +
  geom_errorbarh(aes(xmin = data[1, 3], xmax = data[2, 3]), size = 2, colour = "dodgerblue4", alpha = 0.5, height = 0) +
  geom_point(aes(x = data[, 4]), size = 5, colour = "brown", alpha = 0.5) +
  labs(y = "", x = "")


