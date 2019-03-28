######################################
#   Code for conducting a Bayesian   #
# analysis where the predicted       #
# variable is count-based. This      #
# specific example has two nominal   #
# predictor variables.               #
#                                    #
#                 by                 #
#             Tim Frasier            #
######################################


#-------------------------------#
# Load the required libraries   #
# and files.                    #
#-------------------------------#
library(rstan)
options(mc.cores = parallel::detectCores())
library(ggplot2)
source("plotPost.R")


#--------------------------#
# Read data into R.        #
#--------------------------#
haireye = read.table("haireye.csv", header = TRUE, sep = ",")


#---------------------------#
# Prepare the data for Stan #
#---------------------------#

# Y-Data
y = as.integer(haireye$Freq)
N = length(y)
yLogMean = log(mean(y))
yLogSD = log(sd(c(rep(0, N - 1), sum(y))))

# Eye data
eye = as.numeric(haireye$Eye)
eyeColours = levels(haireye$Eye)
nEyeColours = length(unique(eye))

# Hair data
hair = as.numeric(haireye$Hair)
hairColours = levels(haireye$Hair)
nHairColours = length(unique(hair))


#-----------------------#
# Prepare data for Stan #
#-----------------------#
dataList = list(
  y = y,
  N = N,
  yLogMean = yLogMean,
  yLogSD = yLogSD,
  eye = eye,
  hair = hair,
  nEyeColours = nEyeColours,
  nHairColours = nHairColours
)


#--------------------------#
#     Define the model     #
#--------------------------#
modelstring = "
  data {
    int N;              // Sample size
    int nEyeColours;    // Number of different eye colours in data set
    int nHairColours;   // Number of different hair colours in data set
    real yLogMean;      // The log mean of the observed y data
    real yLogSD;        // The log SD of the observed y data
    int<lower=0> y[N];  // The y data, remember they are integers, and must be defined as such
    int eye[N];         // The eye data, contains indicators of eye colour
    int hair[N];        // The hair data, contains indicators of hair colour
  }

  parameters {
    real b0;
    real b1[nEyeColours];                 // Effect of eye colour
    real b2[nHairColours];                // Effect of hair colour
    real b3[nEyeColours, nHairColours];   // Interaction effect between hair and eye colour 
  }

  model {
    // Definitions
    vector[N] lambda;

    // Likelihood
    for (i in 1:N) {
      lambda[i] = exp(b0 + b1[eye[i]] + b2[hair[i]] + b3[eye[i], hair[i]]);
      y[i] ~ poisson(lambda[i]);
    }

    // Priors
    b0 ~ normal(yLogMean, yLogSD);

    for (j in 1:nEyeColours) {
      b1[j] ~ normal(0, 1);
    }

    for (j in 1:nHairColours) {
      b2[j] ~ normal(0, 1);
    }

    for (j in 1:nEyeColours) {
      for (k in 1:nHairColours) {
        b3[j, k] ~ normal(0, 1);
      }
    }
  }

  generated quantities {
    // Definitions
    vector[N] lambda_pred;
    vector[N] y_pred;

    for (i in 1:N) {
      lambda_pred[i] = exp(b0 + b1[eye[i]] + b2[hair[i]] + b3[eye[i], hair[i]]);
      y_pred[i] = poisson_rng(lambda_pred[i]);
    }
  }
"
writeLines(modelstring, con = "model.stan")


#--------------------------#
#     run STAN             #
#--------------------------#
stanFit <- stan(file = "model.stan", 
                data = dataList, 
                pars = c("b0", "b1", "b2", "b3", "y_pred"),
                warmup = 2000,
                iter = 10000, 
                chains = 3)


#-------------------------#
# Check MCMC Performance  #
#-------------------------#
print(stanFit)
stan_trace(stanFit, pars = "b0", inc_warmup = TRUE)
stan_trace(stanFit, pars = "b1", inc_warmup = TRUE)
stan_trace(stanFit, pars = "b2", inc_warmup = TRUE)
stan_trace(stanFit, pars = "b3", inc_warmup = TRUE)


#-----------------------------------#
#  Plot Results With Stan Functions #
#-----------------------------------#
stan_plot(stanFit, par = c("b0", "b1", "b2", "b3"))


#----------------------------------#
#        Extract the data          #
#----------------------------------#
mcmcChains = as.data.frame(stanFit)

b0 = mcmcChains[, "b0"]
chainLength = length(b0)

b1 = matrix(0, nrow = chainLength, ncol = nEyeColours)
for (i in 1:nEyeColours) {
  b1[, i] = mcmcChains[, paste("b1[", i, "]", sep = "")]
}

b2 = matrix(0, nrow = chainLength, ncol = nHairColours)
for (i in 1:nHairColours) {
  b2[, i] = mcmcChains[, paste("b2[", i, "]", sep = "")]
}

b3 = matrix(0, nrow = chainLength, ncol = (nEyeColours * nHairColours))
column = 1
for (i in 1:nEyeColours) {
  for (j in 1:nHairColours) {
    b3[, column] = mcmcChains[, paste("b3[", i, ",", j, "]", sep = "")]
    column = column + 1
  }  
}


#--------------------------------#
#  Plot Posterior Distributions  #
#--------------------------------#
# b0
par(mfrow = c(1, 1))
histInfo = plotPost(b0, xlab = bquote(beta[0]))

# b1
par(mfrow = c(2, 2))
for (i in 1:nEyeColours) {
    histInfo = plotPost(b1[, i], xlab = bquote(beta * 1[.(i)]), main = paste("b1:", eyeColours[i]))
}

# b2
par(mfrow = c(2, 2))
for (i in 1:nHairColours) {
    histInfo = plotPost(b2[, i], xlab = bquote(beta * 2[.(i)]), main = paste("b2:", hairColours[i]))
}

# b3
par(mfrow = c(2, 2))
nCombos = nEyeColours * nHairColours
counter = 1
while (counter <= nCombos) {
    for (i in 1:nEyeColours) {
        for (j in 1:nHairColours) {
            histInfo = plotPost(b3[, counter], xlab = bquote(beta * 1[.(i)] * beta * 2[.(j)]), main = paste("b1:", eyeColours[i], ", b2:", hairColours[j]))
            counter = counter + 1
        }
    }
}


#------------------------------------#
#  Posterior Predictive Check        #
#------------------------------------#

#--- Extract predicted data ---#
yPred = matrix(0, nrow = chainLength, ncol = N)

for (i in 1:N) {
  yPred[, i] = mcmcChains[, paste("y_pred[", i, "]", sep = "")]
}  

#--- Calculate the mean and HDI values ---#
yPredMean = apply(yPred, 2, mean)
yPredLow = apply(yPred, 2, quantile, probs = 0.025)
yPredHigh = apply(yPred, 2, quantile, probs = 0.975)

#--- Plot predicted values ---#
par(mfrow = c(1, 1))
combination = 1:N
dotchart(x = yPredMean, labels = combination, xlim = c(min(yPredLow), max(yPredHigh)))

#--- Add HDI lines ---#
segments(x0 = yPredLow, y0 = combination, x1 = yPredHigh, y1 = combination)

#--- Add observed values ---#
points(x = y, y = combination, pch = 16, col = rgb(1, 0, 0, 0.6))
