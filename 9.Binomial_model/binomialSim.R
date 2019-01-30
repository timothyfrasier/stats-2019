#######################################
#           binomialSim.R             #
#                                     #
# This function is a simple MCMC-like #
# sampler for teaching. It takes a    #
# number of trials, the number        #
# of successes, and the number of     #
# steps to take, and then estimates   #
# the probability value based on a    #
# binomial distribution and a uniform #
# Prior.                              #
#######################################


binomialCalc = function(nTrials, nSuccesses, nSteps) {

#--- Generate binary data ---#
series = c(rep(1, times = nSuccesses), rep(0, times = (nTrials - nSuccesses)))
series = sample(series, size = nTrials, replace = FALSE)

#------------------------------------#
# Create model to estimate parameter #
#------------------------------------#

#--- Create holder for results ---#
thetaEstimates = rep(NA, times = nSteps)
thetaSteps = rep(NA, times = nSteps)

for (i in 1:nSteps) {
  #--- Set priors ---#
  thetaPrior = runif(1, min = 0, max = 1)
  thetaSteps[i] = thetaPrior
  
  #--- Generate binomial distribution ---#
  simOut = rbinom(n = nTrials, size = 1, prob = thetaPrior)
  
  #--- Count Successes ---#
  simSuccess = sum(simOut > 0)
  
  #--- Store theta if it fits observed number ---#
  if (simSuccess == nSuccesses) {
    thetaEstimates[i] = thetaPrior 
  } else {
    thetaEstimates[i] = NA
  }
}

#--- Plot results ---#
par(mfrow = c(1, 2))
hist(thetaSteps, main = "", xlab = "Proposed Values", col = rgb(0, 0, 1, 0.3))
dist = hist(thetaEstimates, main = "", col = rgb(0, 0, 1, 0.3), xlab = "Accepted Values")

intervals = quantile(thetaEstimates, c(0.05, 0.95), na.rm = TRUE)
low = as.numeric(intervals[1])
high = as.numeric(intervals[2])

segments(x0 = low, y0 = (0.02 * max(dist$counts)), x1 = high, y1 = (0.02 * max(dist$counts)), lwd = 4)
meanValue = mean(thetaEstimates, na.rm = TRUE)
meanValue = round(meanValue, digits = 4)
low = round(low, digits = 3)
high = round(high, digits = 3)
text(x = meanValue, font = 2, y = max(dist$counts), labels = meanValue)
text(x = low, y = (0.1 * max(dist$counts)), labels = low)
text(x = high, y = (0.1 * max(dist$counts)), labels = high)

#--- Plot accepted steps ---#
noStep = sum(is.na(thetaEstimates))
usedSteps = nSteps - noStep

text(x = (0.8 * max(dist$breaks)), y = (0.9 * max(dist$counts)), labels = "Steps:")
text(x = (0.8 * max(dist$breaks)), y = (0.82 * max(dist$counts)), labels = bquote(.(usedSteps))) 
}
