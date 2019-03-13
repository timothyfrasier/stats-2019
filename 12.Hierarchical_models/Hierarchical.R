###########################################################
# Code for Bayesian analysis of data with one metric      #
# predicted variable and one nominal predictor variable.  #
# This code uses a normal distribution to evaluate the    #
# variation around the mean. It uses a hierarchical model #
# to estimate means and standard deviations for each group#
# Based on code from Kruschke (2011, 2015). Uses scripts  #
# "plotPost.R" and "HDIofMCMC.R" from Kruschke (2011).    #
#                                                         #
#                            by                           #
#                       Tim Frasier                       #
#                                                         #
#                Last Updated: 05-Mar-2019                #
###########################################################

#-----------------------#
# Read the data into R  #
#-----------------------#
library(rstan)
library(ggplot2)
source("plotPost.R")
coffee = read.table("coffee.csv", header = TRUE, sep = ",")


#------------------------------#
#         Plot the Data        #
#------------------------------#

#--- Box Plot ---#
ggplot(coffee) +
  theme_bw() +
  geom_boxplot(aes(x = shop, y = time), fill = "dodgerblue4", alpha = 0.6)

#--- Jittered Points ---#
ggplot(coffee) +
  theme_bw() +
  geom_jitter(aes(x = shop, y = time), height = 0, width = 0.3, alpha = 0.6)

#--- Points by Trip ---#
plot(coffee$trip, coffee$time, ylim = c(0, 27), xlab = "Trip", ylab = "Wait Time (min)", pch = 16)
abline(v = 30.5)
abline(v = 35.5)
segments(x0 = 0, y0 = 15, x1 = 30.5, y1 = 15, lwd = 2, lty = 2)
segments(x0 = 30.5, y0 = 5, x1 = 35.5, y1 = 5, lwd = 2, lty = 2)
segments(x0 = 35.5, y0 = 10, x1 = 52, y1 = 10, lwd = 2, lty = 2)


#----------------------------#
# Standardize the y-variable #
# (metric)                   #
#----------------------------#
y = coffee$time
yMean = mean(y)
ySD = sd(y)
zy = (y - yMean) / ySD
N = length(zy)


#---------------------------#
# Organize the categorical  #
# data                      #
#---------------------------#
x = as.numeric(coffee$shop)
xNames = levels(coffee$shop)
nxLevels = length(unique(coffee$shop))


#-----------------------------#
# Create a data list for JAGS #
#-----------------------------#
dataList = list(
  y = zy,
  x = x,
  N = N,
  nxLevels = nxLevels
)


#-----------------------------#
# Define the model and write  #
# a string for Stan           #
#-----------------------------#
modelstring = "
  data {
    int N;          // Sample size
    int nxLevels;   // Number of shop types (levels in our nominal variable)
    vector[N] y;    // Vector of wait times
    int x[N];       // Vector of indicators of which level each sample is in
  }

  parameters {
    real b1[nxLevels];              // Coefficients for effects of being in each level
    real<lower=0> sigma[nxLevels];  // Coefficients for sd of being in each level
    real shopMean;                  // Mean wait time across all shops
    real<lower=0> shopMeanSD;       // sd for mean wait time across all shops
  }

  model {
    // Definitions
    vector[N] mu;

    // Likelihood
    for (i in 1:N) {
      mu[i] = b1[x[i]];
      y[i] ~ normal(mu[i], sigma[x[i]]);
    }

    // Priors
    for (j in 1:nxLevels) {
      b1[j]    ~ normal(shopMean, shopMeanSD);
      sigma[j] ~ cauchy(1, 1);
    }

    // Hyperpriors
    shopMean   ~ normal(0, 1);
    shopMeanSD ~ normal(1, 1);
  }

  generated quantities {
    // Definitions
    vector[N] y_pred;

    for (i in 1:N) {
      y_pred[i] = normal_rng(b1[x[i]], sigma[x[i]]);
    }
  }
"
writeLines(modelstring, con = "model.stan")


#--------------------------#
#     run STAN             #
#--------------------------#
stanFit <- stan(file = "model.stan", 
                data = dataList, 
                pars = c("b1", "sigma", "shopMean", "shopMeanSD", "y_pred"),
                warmup = 2000,
                iter = 7000, 
                chains = 3)


#-------------------------#
# Check MCMC Performance  #
#-------------------------#
print(stanFit)
stan_trace(stanFit, pars = "b1", inc_warmup = TRUE)
stan_trace(stanFit, pars = "sigma", inc_warmup = TRUE)

#-----------------------------------#
#  Plot Results With Stan Functions #
#-----------------------------------#
stan_plot(stanFit, par = "b1")
stan_plot(stanFit, par = "sigma")


#--------------------#
# Extract the Data   #
#--------------------#
mcmcChains = as.data.frame(stanFit)

zshopMean = mcmcChains$shopMean
zshopMeanSD = mcmcChains$shopMeanSD

chainLength = length(mcmcChains[, 1])

zsigma = matrix(0, ncol = nxLevels, nrow = chainLength)
for (i in 1:nxLevels) {
  zsigma[, i] = mcmcChains[, paste("sigma[", i, "]", sep = "")]
}

zb1 = matrix(0, ncol = nxLevels, nrow = chainLength)
for (i in 1:nxLevels) {
  zb1[, i] = mcmcChains[, paste("b1[", i, "]", sep = "")]
}

zypred = matrix(0, ncol = N, nrow = chainLength)
for (i in 1:N) {
  zypred[, i] = mcmcChains[, paste("y_pred[", i, "]", sep = "")]
}


#--------------------------------#
# Convert back to original scale #
#--------------------------------#
sigma = zsigma * ySD
b1 = (zb1 * ySD) + yMean
ypred = (zypred * ySD) + yMean
shopMean = (zshopMean * ySD) + yMean
shopMeanSD = zshopMeanSD * ySD


#---------------------#
#   Plot Estimates    #
#---------------------#
par(mfrow = c(1, 3))

#--- Mean Wait Time ---#
histInfo = plotPost(b1[, 1], xlab = "Wait time: Second Cup")
histInfo = plotPost(b1[, 2], xlab = "Wait time: Starbucks")
histInfo = plotPost(b1[, 3], xlab = "Wait time: Tim Hortons")

#--- s.d. of Wait Time ---#
histInfo = plotPost(sigma[, 1], xlab = "s.d.: Second Cup", showMode = TRUE)
histInfo = plotPost(sigma[, 2], xlab = "s.d.: Starbucks", showMode = TRUE)
histInfo = plotPost(sigma[, 3], xlab = "s.d.: Tim Hortons", showMode = TRUE)

#--- Overall mean and s.d. ---#
par(mfrow = c(1, 2))
histInfo = plotPost(shopMean, xlab = "Mean wait time")
histInfo = plotPost(shopMeanSD, xlab = "s.d. of Mean wait time", showMode = TRUE)

#--- What does this look like? ---#
x = seq(from = 0, to = 40, length = 200)
y = dnorm(x, mean = 10.671, sd = 5.641)
plot(y ~ x, type = "l", ylim = c(0, 0.15), lwd = 3, col = rgb(0, 0, 1, 0.6))

#--- As a sampling from the posterior ---#
samps1 = seq(from = 1, to = length(shopMean), length = 20)
samps2 = round(samps1)
for (i in samps2) {
  y = dnorm(x, mean = shopMean[i], sd = shopMeanSD[i])
  lines(x, y, col = rgb(0, 0, 0, 0.3))
}


#-------------------------------#
#  Posterior predictive check   #
#-------------------------------#

#--- Mean expected value for each visit to coffee shop ---#
ypredMean = apply(ypred, 2, mean)

#--- Upper and lower expected 95% HDI for each visit ---#
ypredLow  = apply(ypred, 2, quantile, probs = 0.025)
ypredHigh = apply(ypred, 2, quantile, probs = 0.975)

#--- Plot original values ---#
par(mfrow = c(1, 1))
plot(coffee$trip, coffee$time, ylim = c(0, 27), xlab = "Trip", ylab = "Wait Time (min)", pch = 16, col = rgb(0, 0, 0, 0.6))
abline(v = 30.5)
abline(v = 35.5)
segments(x0 = 0, y0 = 15, x1 = 30.5, y1 = 15, lwd = 2, lty = 2, col = "red")
segments(x0 = 30.5, y0 = 5, x1 = 35.5, y1 = 5, lwd = 2, lty = 2, col = "red")
segments(x0 = 35.5, y0 = 10, x1 = 52, y1 = 10, lwd = 2, lty = 2, col = "red")

#--- Add mean predicted values ---#
par(new = TRUE)
plot(coffee$trip, ypredMean, ylim = c(0, 27), axes = FALSE, xlab = "", ylab = "", main = "", pch = 16, col = rgb(0, 0, 1, 0.5))

#--- Add HDI lines ---#
for (i in 1:N) {
  segments(x0 = coffee$trip[i], y0 = ypredLow[i], x1 = coffee$trip[i], y1 = ypredHigh[i], col = rgb(0, 0, 1, 0.5))
}