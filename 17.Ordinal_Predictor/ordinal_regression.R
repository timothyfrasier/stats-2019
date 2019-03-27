####################################################
# Code for conducting Bayesian regression with an  #
# ordinal predicted variable.                      #
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
ord = read.table("ordinalData.csv", header = TRUE, sep = ",")


#--------------------------#
# Plot the Data            #
#--------------------------#
# Pairs
pairs(ord, pch = 16, col = rgb(0, 0, 1, 0.3))

# As a histogram of each oridinal value
par(mfrow = c(1, 1))
yTable = table(ord$Y)
yTable
yTable.df = data.frame(yTable)
yTable.df[, 1] = as.numeric(as.character(yTable.df[, 1]))
plot(yTable.df[, 1], yTable.df[, 2], type = "h", ylab = "Frequency", xlab = "response", lwd = 4, col = rgb(0, 0, 1, 0.5))


#--------------------------------------#
# Convert to the Cumulative Scale      #
#--------------------------------------#

#--- Compute Cumulative Probabilities ---#
# Get proportions
pr_y = yTable / nrow(ord)

# Get cumulative proportions
cum_pr_y = cumsum(pr_y)

# Plot
plot(yTable.df[, 1], cum_pr_y, type = "b", lwd = 2, ylab = "cumulative proportion", xlab = "response", ylim = c(0, 1), col = "blue")


#---------------------------#
# Prepare the data for Stan #
#---------------------------#
y = ord$Y
N = length(y)
nLevels = length(unique(y))

x1 = ord$X1
x2 = ord$X2

#--- Standardize x values ---#
x1Mean = mean(x1)
x1SD = sd(x1)
zx1 = (x1 - x1Mean) / x1SD

x2Mean = mean(x2)
x2SD = sd(x2)
zx2 = (x2 - x2Mean) / x2SD


#-----------------------#
# Prepare data for Stan #
#-----------------------#
dataList = list(
  N = N,
  K = nLevels,
  y = y,
  x1 = zx1,
  x2 = zx2
)

#--------------------------#
#     Define the model     #
#--------------------------#
modelstring = "
  data {
    int<lower=0> N;               // Sample size
    int<lower=2> K;               // Number of ordered outcomes
    int<lower=1, upper=K> y[N];   // Ordered predicted variable
    vector[N] x1;                 // First predictor variable
    vector[N] x2;                 // Second predictor variable
  }

  parameters {
    real b0;
    real b1;              // Effect of first predictor variable
    real b2;              // Effect of second predictor variable
    ordered[K-1] c;       // Vector of cutpoints
  }

  model {
    // Definitions
    vector[N] mu;

    // Likelihood
    for (i in 1:N) {
      mu[i] = b0 + (b1 * x1[i]) + (b2 * x2[i]);
      y[i] ~ ordered_logistic(mu[i], c);
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
      mu_pred[i] = b0 + (b1 * x1[i]) + (b2 * x2[i]);
      y_pred[i] = ordered_logistic_rng(mu_pred[i], c);
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


#--------------------------------#
#      Extract the Data          #
#--------------------------------#
mcmcChain = as.data.frame(stanFit)

zb0 = mcmcChain[, "b0"]
zb1 = mcmcChain[, "b1"]
zb2 = mcmcChain[, "b2"]


#--------------------------------#
#  Convert to Original Scale     #
#--------------------------------#
b1 = zb1 / x1SD
b2 = zb2 / x2SD
b0 = zb0 - (((zb1 * x1Mean) / x1SD) + ((zb2 * x2Mean) / x2SD))


#--------------------------------#
#  Plot Posterior Distributions  #
#--------------------------------#
par(mfrow = c(1, 3))
histInfo = plotPost(b0, xlab = bquote(beta[0]))
histInfo = plotPost(b1, xlab = bquote(beta[1]), main = "x1")
histInfo = plotPost(b2, xlab = bquote(beta[2]), main = "x2")


#--------------------------------#
#  Posterior Predictive Check    #
#--------------------------------#

#--- Get the predicted values ---#
chainLength = length(mcmcChain[, 1])

zypred = matrix(0, ncol = N, nrow = chainLength)
for (i in 1:N) {
  zypred[, i] = mcmcChain[, paste("y_pred[", i, "]", sep = "")]
}

#--- Mean expected value for each visit to coffee shop ---#
zypredMean = apply(zypred, 2, mean)

#--- Upper and lower expected 95% HDI for each visit ---#
zypredLow  = apply(zypred, 2, quantile, probs = 0.025)
zypredHigh = apply(zypred, 2, quantile, probs = 0.975)

source("HDIofMCMC.R")

# Create a matrix to hold results
yPostPred <- matrix(0, nrow = chainLength, ncol = nLevels)

# For each step in the chain...
for (i in 1:chainLength) {
  
  # Initialize holders (counters for each level)
  counter1 <- 0
  counter2 <- 0
  counter3 <- 0
  counter4 <- 0
  counter5 <- 0
  counter6 <- 0
  counter7 <- 0
  
  # For each individual...
  for (j in 1:N) {
    
    # Calculate mean
    mu = b0[i] + (b1[i] * x1[j]) + (b2[i] * x2[j])
    
    # Calculate their probability of being in each level
    levelProbs <- rep(0, times = nLevels)
    levelProbs[1] <- pnorm(alpha[i, 1], mu, sigma[i])
    levelProbs[2] <- pnorm(alpha[i, 2], mu, sigma[i]) - pnorm(alpha[i, 1], mu, sigma[i])
    levelProbs[3] <- pnorm(alpha[i, 3], mu, sigma[i]) - pnorm(alpha[i, 2], mu, sigma[i])
    levelProbs[4] <- pnorm(alpha[i, 4], mu, sigma[i]) - pnorm(alpha[i, 3], mu, sigma[i])
    levelProbs[5] <- pnorm(alpha[i, 5], mu, sigma[i]) - pnorm(alpha[i, 4], mu, sigma[i])
    levelProbs[6] <- pnorm(alpha[i, 6], mu, sigma[i]) - pnorm(alpha[i, 5], mu, sigma[i])
    levelProbs[7] <- 1 - pnorm(alpha[i, 6], mu, sigma[i])
    
    # Find item number for highest value
    levelID <- which.max(levelProbs)
    
    # Increase counter for appropriate group
    if(levelID == 1) {
      counter1 <- counter1 + 1
    } else {
      if (levelID == 2) {
        counter2 <- counter2 + 1
      } else {
        if (levelID == 3) {
          counter3 <- counter3 + 1
        } else {
          if (levelID == 4) {
            counter4 <- counter4 + 1
          } else {
            if (levelID == 5) {
              counter5 <- counter5 + 1
            } else {
              if (levelID == 6 ) {
                counter6 <- counter6 + 1
              } else {
                counter7 <- counter7 + 1
              }  
            }
          }
        }
      }
    }
    
    # Place results in results matrix
    yPostPred[i, 1] <- counter1
    yPostPred[i, 2] <- counter2
    yPostPred[i, 3] <- counter3
    yPostPred[i, 4] <- counter4
    yPostPred[i, 5] <- counter5
    yPostPred[i, 6] <- counter6
    yPostPred[i, 7] <- counter7
    
  }
}

yPredMeans = apply(yPostPred, 2, median, na.rm = TRUE)
yPredHDI = apply(yPostPred, 2, HDIofMCMC)

# Plot original data
hist(y, breaks = c(0.5, (1:nLevels + 0.5)), main = "", col = "skyblue", border = "white")

# Add predicted means
points(x = 1:nLevels, y = yPredMeans, pch = 16)

# Add HDI bars
segments(x0 = 1:nLevels, y0 = yPredHDI[1, ], x1 = 1:nLevels, y1 = yPredHDI[2, ], lwd = 2)
