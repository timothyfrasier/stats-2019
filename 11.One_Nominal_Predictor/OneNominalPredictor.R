################################
# Code for conducting Bayesian #
# regressioon with a metric    #
# predicted variable and one   #
# nominal predictor variable.  #
################################

###############################
# LOAD APPROPRIATE LIBRARIES  #
###############################
library(rstan)
library(ggplot2)
library(ggridges)
source("plotPost.R")


################################
# WHEN THERE ARE TWO LEVELS IN #
# THE NOMINAL CATEGORY.        #
################################

#---------------------#
#   Load the data     #
#---------------------#
howell = read.table("Howell.csv", header = TRUE, sep = ",")
str(howell)
summary(howell)


#---------------------#
#   Plot the data     #
#---------------------#

#--- Histogram ---#
ggplot(howell) +
  theme_bw() +
  geom_histogram(aes(x = howell$height), fill = "dodgerblue4", alpha = 0.6)

#--- Box Plot ---#
ggplot(howell) +
  theme_bw() +
  geom_boxplot(aes(x = sex, y = height), fill = "dodgerblue4", alpha = 0.6)

#--- Point Plot ---#
ggplot(howell) +
  theme_bw() +
  geom_point(aes(x = sex, y = height), colour = "dodgerblue4", alpha = 0.6)

#--- Point Plot Jittered ---#
ggplot(howell) +
  theme_bw() +
  geom_jitter(aes(x = sex, y = height), height = 0, colour = "dodgerblue4", alpha = 0.6)

#--- Ridges ---#
ggplot(howell) +
  theme_bw() +
  geom_density_ridges(aes(x = height, y = sex, group = sex), fill = "dodgerblue4", alpha = 0.6)

#--- Ridges With Different Colours---#
ggplot(howell) +
  theme_bw() +
  geom_density_ridges(aes(x = height, y = sex, group = sex, fill = sex), alpha = 0.6) +
  scale_fill_manual(values = c("brown", "dodgerblue4")) +
  theme(legend.position = "none")


#-----------------------#
# Frequentist Approach  #
#-----------------------#

#--- Organize the data ---#
yFemale = howell[howell$sex == "female", 1]
yMale = howell[howell$sex == "male", 1]

#--- t-test ---#
t.test(yFemale, yMale)

#--- linear model ---#
model = lm(howell$height ~ howell$sex)
summary(model)


#-----------------------#
#   Bayesian Approach   #
#-----------------------#

#--- Prepare and Standardize the Metric Data ---#
y = howell$height
N = length(y)

yMean = mean(y)
ySD = sd(y)
zy = (y - yMean) / ySD

#--- Prepare the nominal data ---#
x = as.numeric(howell$sex)
x

xNames = levels(howell$sex)
xNames

nxLevels = length(xNames)
nxLevels


#--- Prepare the data as a list for Stan ---#
dataList = list (
  y = zy,
  N = N,
  x = x,
  nxLevels = nxLevels
)


#------------------------------#
#       Define the Model       #
#------------------------------#
modelString = "
  data {
    int N;
    int nxLevels;
    vector[N] y;
    int x[N];    // Note that indices like this can't be vectors
  }

  parameters {
    real b1[nxLevels];    // A different b1 coefficient for each level in the variable
    real<lower=0> sigma;  // A single sigma value
  }

  model {
    // Definitions
    vector[N] mu;

    // Likelihood
    for (i in 1:N) {
      mu[i] = b1[x[i]];
      y[i] ~ normal(mu[i], sigma);
    }  

    // Priors
    for (j in 1:nxLevels) {
      b1[j] ~ normal(0, 1);
    }

    sigma ~ cauchy(1, 1);
  }

  generated quantities {
    vector[N] y_pred;

    for (i in 1:N) {
      y_pred[i] = normal_rng(b1[x[i]], sigma);
    } 
  }
"
writeLines(modelString, con="model.stan")


#--------------------------#
#     run STAN             #
#--------------------------#
stanFit <- stan(file = "model.stan", 
                data = dataList, 
                pars = c("b1", "sigma", "y_pred"),
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


#-------------------------#
# Extract the data        #
#-------------------------#
mcmcChains = as.data.frame(stanFit)
zsigma = mcmcChains$sigma

chainLength = length(zsigma)
zb1 = matrix(0, ncol = nxLevels, nrow = chainLength)
for (i in 1:nxLevels) {
  zb1[, i] = mcmcChains[, paste("b1[", i, "]", sep = "")]
}

ypred = matrix(0, ncol = N, nrow = chainLength)
for (i in 1:N) {
  ypred[, i] = mcmcChains[, paste("y_pred[", i, "]", sep = "")]
}

#---------------------------------#
# Plot using Kruschke's Functions #
#---------------------------------#
histInfo = plotPost(zsigma, xlab = bquote(sigma))

par(mfrow = c(1, 2))
histInfo = plotPost(zb1[, 1], xlab = bquote(beta[females]))
histInfo = plotPost(zb1[, 2], xlab = bquote(beta[males]))

#--- Calculate and plot difference ---#
difference = zb1[, 2] - zb1[, 1]
par(mfrow = c(1, 1))
histInfo = plotPost(difference, xlab = bquote(beta[males] - beta[females]))
abline(v = 0, lwd = 2, lty = 3, col = "red")


#------------------------------#
#     Plot using ggridges      #
#------------------------------#

#--- Organize the data ---#
zbCombined = c(zb1[, 1], zb1[, 2])
sex = c(rep("female", times = length(zb1[, 1])), rep("male", times = length(zb1[, 1])))
combined = data.frame(zbCombined, sex)

#--- Plot the data ---#
ggplot(combined) +
  theme_bw() +
  geom_density_ridges(aes(x = zbCombined, y = sex, group = sex, fill = sex), alpha = 0.6) +
  scale_fill_manual(values = c("brown", "dodgerblue4")) +
  theme(legend.position = "none")


#--------------------------------#
# Convert back to original scale #
#--------------------------------#
sigma = zsigma * ySD
b1 = yMean + (zb1 * ySD)
par(mfrow = c(1, 2))
histInfo = plotPost(b1[, 1], xlab = bquote(beta[females]))
histInfo = plotPost(b1[, 2], xlab = bquote(beta[males]))

#--- Calculate and plot difference ---#
difference = b1[, 2] - b1[, 1]
par(mfrow = c(1, 1))
histInfo = plotPost(difference, xlab = bquote(beta[males] - beta[females]))
abline(v = 0, lwd = 2, lty = 3, col = "red")


#------------------------------#
#  Posterior Predictive Check  #
#------------------------------#

#--- Calculate the mean, 95% high, and 95% low expected values for each individual ---#
predMean = apply(ypred, 2, mean)
predLow = apply(ypred, 2, quantile, probs = 0.025)
predHigh = apply(ypred, 2, quantile, probs = 0.975)

howellCombined = data.frame(zy, howell$sex, predMean, predLow, predHigh)

#--- Plot the data ---#
ggplot(howellCombined) +
  theme_bw() +
  geom_jitter(aes(x = howell.sex, y = zy), height = 0, alpha = 0.6) +
  geom_pointrange(aes(x = howell.sex, y = predMean, ymin = predLow, ymax = predHigh), size = 1, colour = "brown", alpha = 0.007)


#----------------------------------------------------------------------------------#

#################################
# Analysis with Multiple Levels #
#################################

#--- Read Data into R ---#
milk = read.table("milk.csv", header = TRUE, sep = ",")
summary(milk$clade)
summary(milk$kcal.per.g)

#--- Plot the data ---#
ggplot(milk) +
  theme_bw() +
  geom_boxplot(aes(x = clade, y = kcal.per.g), fill = "dodgerblue4", alpha = 0.6)


#------------------------#
# Frequentist Approach   #
#------------------------#

#--- ANOVA ---#
anovatest = aov(milk$kcal.per.g ~ milk$clade)
print(model.tables(anovatest))
summary(anovatest)

#--- Linear Regression ---#
model = lm(milk$kcal.per.g ~ milk$clade)
summary(model)


#------------------------#
#   Bayesian Approach    #
#------------------------#

#--- Organize the data ---#
y = milk$kcal.per.g
N = length(y)

x = as.numeric(milk$clade)
xNames = levels(milk$clade)
nxLevels = length(xNames)

#--- Create as a list to send to Stan ---#
dataList = list (
  y = y,
  N = N,
  x = x,
  nxLevels = nxLevels
)

#--------------------------#
#     run STAN             #
#--------------------------#
stanFit <- stan(file = "model.stan", 
                data = dataList, 
                pars = c("b1", "sigma", "y_pred"),
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


#-------------------------#
# Extract the data        #
#-------------------------#
mcmcChains = as.data.frame(stanFit)
zsigma = mcmcChains$sigma

chainLength = length(zsigma)
zb1 = matrix(0, ncol = nxLevels, nrow = chainLength)
for (i in 1:nxLevels) {
  zb1[, i] = mcmcChains[, paste("b1[", i, "]", sep = "")]
}

ypred = matrix(0, ncol = N, nrow = chainLength)
for (i in 1:N) {
  ypred[, i] = mcmcChains[, paste("y_pred[", i, "]", sep = "")]
}

#------------------------------#
#     Plot using ggridges      #
#------------------------------#

#--- Organize the data ---#
zbCombined = c(zb1[, 1], zb1[, 2], zb1[, 3], zb1[, 4])
clade = c(rep("Ape", times = length(zb1[, 1])), rep("NWM", times = length(zb1[, 1])), rep("OWM", times = length(zb1[, 1])), rep("STP", times = length(zb1[, 1])))
combined = data.frame(zbCombined, clade)

#--- Plot the data ---#
ggplot(combined) +
  theme_bw() +
  geom_density_ridges(aes(x = zbCombined, y = clade, group = clade, fill = clade), alpha = 0.6) +
  scale_fill_manual(values = c("brown", "dodgerblue4", "chartreuse4", "coral3")) +
  theme(legend.position = "none")


#------------------------------#
#  Posterior Predictive Check  #
#------------------------------#

#--- Calculate the mean, 95% high, and 95% low expected values for each individual ---#
predMean = apply(ypred, 2, mean)
predLow = apply(ypred, 2, quantile, probs = 0.025)
predHigh = apply(ypred, 2, quantile, probs = 0.975)

primateCombined = data.frame(y, milk$clade, predMean, predLow, predHigh)

#--- Plot the data ---#
ggplot(primateCombined) +
  theme_bw() +
  geom_jitter(aes(x = milk$clade, y = y), height = 0, alpha = 0.6) +
  geom_pointrange(aes(x = milk$clade, y = predMean, ymin = predLow, ymax = predHigh), size = 1, colour = "brown", alpha = 0.2)

