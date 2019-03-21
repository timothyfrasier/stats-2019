####################################################
# Code for conducting Bayesian logistic regression #
# with two metric predictor variables.             #
####################################################


#-------------------------------#
# Load the required libraries   #
# and files.                    #
#-------------------------------#
library(rstan)
library(ggplot2)
source("plotPost.R")


#--------------------------#
# Read data into R.        #
#--------------------------#
htwt = read.table("HtWt2.csv", header = TRUE, sep = ",")


#--------------------------#
# Plot the Data            #
#--------------------------#
ggplot(htwt) +
  theme_bw() +
  geom_boxplot(aes(x = sex, y = height), fill = "dodgerblue4", alpha = 0.6)

ggplot(htwt) +
  theme_bw() +
  geom_boxplot(aes(x = sex, y = weight), fill = "dodgerblue4", alpha = 0.6)


#---------------------------#
# Organize the Data         #
#---------------------------#
y = as.numeric(htwt$sex)  # female will be 1 and male will be 2
y = y - 1                 # female will be 0 and male will be 1 (logistic regression requires...
                          # ...predicted variables to be 1 or 0)
N = length(y)

#--- height ---#
height = htwt$height
heightMean = mean(height)
heightSD = sd(height)
zheight = (height - heightMean) / heightSD

#--- Weight ---#
weight = htwt$weight
weightMean = mean(weight)
weightSD = sd(weight)
zweight = (weight - weightMean) / weightSD


#-----------------------#
# Prepare data for Stan #
#-----------------------#
dataList = list(
  y = y,
  N = N,
  height = zheight,
  weight = zweight
)

#--------------------------#
#     Define the model     #
#--------------------------#
modelstring = "
  data {
    int N;                // Sample size
    int y[N];             // Indices of 0s and 1s (predicted variable)
    vector[N] height;     // Vector of height data
    vector[N] weight;     // Vector of weight data
  }

  parameters {
    real b0;              
    real b1;    // Effect of height on prob. of sex
    real b2;    // Effect of weight on prob. of sex
  }

  model {
    // Definitions
    vector[N] mu;

    // Likelihood
    for (i in 1:N) {
      mu[i] = b0 + (b1 * height[i]) + (b2 * weight[i]);
      y[i] ~ bernoulli_logit(mu[i]);
    }

    // Priors
    b0 ~ normal(0, 1);
    b1 ~ normal(0, 1);
    b2 ~ normal(0, 1);
  }

  generated quantities {
    // Definitions
    vector[N] mu_pred;
    vector[N] y_pred;

    for (i in 1:N) {
      mu_pred[i] = b0 + (b1 * height[i]) + (b2 * weight[i]);
      y_pred[i] = bernoulli_logit_rng(mu_pred[i]);
    }
  }
"
writeLines(modelstring, con = "model.stan")


#--------------------------#
#     run STAN             #
#--------------------------#
stanFit = stan(file = "model.stan", 
                data = dataList, 
                pars = c("b0", "b1", "b2", "y_pred"),
                warmup = 2000,
                iter = 7000, 
                chains = 3)


#-------------------------#
# Check MCMC Performance  #
#-------------------------#
print(stanFit)
stan_trace(stanFit, pars = c("b0", "b1", "b2"), inc_warmup = TRUE)

#-----------------------------------#
#  Plot Results With Stan Functions #
#-----------------------------------#
stan_plot(stanFit, par = c("b0", "b1", "b2"))


#--------------------#
# Extract the Data   #
#--------------------#
mcmcChains = as.data.frame(stanFit)

zb0 = mcmcChains[, "b0"]
zb1 = mcmcChains[, "b1"]
zb2 = mcmcChains[, "b2"]


#--------------------------------#
#  Convert to Original Scale     #
#--------------------------------#
b1 = zb1 / heightSD
b2 = zb2 / weightSD
b0 = zb0 - (((zb1 * heightMean) / heightSD) + ((zb2 * weightMean) / weightSD))


#--------------------------------#
#  Plot Posterior Distributions  #
#--------------------------------#
par(mfrow = c(1, 3))
histInfo = plotPost(b0, xlab = bquote(beta[0]))
histInfo = plotPost(b1, xlab = bquote(beta[1]), main = "Height")
histInfo = plotPost(b2, xlab = bquote(beta[2]), main = "Weight")


#--------------------------------#
#   Assess model fit             #
#--------------------------------#

#--- One Predictor Variable At A Time ---#
# Height
par(mfrow = c(1, 1))
chainLength = length(b0)
plot(height, y, xlab = "Height", ylab = "Sex", pch = 16, cex = 1.5, col = rgb(0, 0, 0, 0.5))
abline(h = 0.5, lty = "dotted", lwd = 2)
xNew = floor(seq(1, chainLength, length = 30))
xRange = max(height) - min(height)
xComb = seq(min(height) - 0.1 * xRange, max(height) + 0.1 * xRange, length = 200)
for (i in xNew) {
	lines(xComb , 1 / (1 + exp(-(b0[i] + (b1[i] * xComb) + (b2[i] * mean(weight))))), lwd = 1.5, col = rgb(0, 0, 1, 0.3))
}

# Weight
plot(weight, y, xlab = "Weight", ylab = "Sex", pch = 16, cex = 1.5, col = rgb(0, 0, 0, 0.5))
abline(h = 0.5, lty = "dotted", lwd = 2)
xNew = floor(seq(1, chainLength, length = 30))
xRange = max(weight) - min(weight)
xComb = seq(min(weight) - 0.1 * xRange, max(weight) + 0.1 * xRange, length = 200)
for (i in xNew) {
	lines(xComb , 1 / (1 + exp(-(b0[i] + (b1[i] * mean(height)) + (b2[i] * xComb)))), lwd = 1.5, col = rgb(0, 0, 1, 0.3))
}


#------------------------------------#
#  Posterior Predictive Check        #
# note that no HDI bars are plotted  #
#------------------------------------#

#--- Set number of values used to test prediction ---#
nPred = 20

#--- Select observed x Values & Organize     ---#
newRows = seq(from = 1, to = NROW(htwt), length = nPred)
newRows = round(newRows)
newObs = y[newRows]

#--- Extract the posterior predictions ---#
chainLength = length(mcmcChains[, 1])
ypred = matrix(0, ncol = N, nrow = chainLength)
for (i in 1:N) {
  ypred[, i] = mcmcChains[, paste("y_pred[", i, "]", sep = "")]
}

#--- Mean expected values  ---#
ypredMean = apply(ypred, 2, mean)

#--- Upper and lower expected 95% HDI  ---#
ypredLow  = apply(ypred, 2, quantile, probs = 0.025)
ypredHigh = apply(ypred, 2, quantile, probs = 0.975)

#--- Subset these for just the 20 chosen people ---#
subypredMean = ypredMean[newRows]
subypredLow = ypredLow[newRows]
subypredHigh = ypredHigh[newRows]

#--- Plot the predicted values ---#
dotchart(subypredMean, labels = 1:nPred, xlim = c(-0.1, 1.1), xlab = "Sex")
abline(v = 0.5, lty = "dotted", lwd = 2)

#--- Plot the true values ---#
points(x = newObs, y = 1:nPred, pch = 16, col = rgb(1, 0, 0, 0.5))
