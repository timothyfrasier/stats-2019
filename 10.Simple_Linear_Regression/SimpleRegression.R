####################################
# STAN model for simple linear     #
# regression.                      #
#                                  #
#             by                   #
#         Tim Frasier              #
#   Last updated: 06-Feb-2019      #
####################################


#-----------------------------#
# Load Appropriate Libraries  #
#-----------------------------#
library(rstan)
library(ggplot2)


#-----------------------------#
#     Read data into R        #
#-----------------------------#
htwt = read.table("HtWt.csv", header = TRUE, sep = ",")

#--- Get a feel for the data ---#
summary(htwt)
sd(htwt$height)
sd(htwt$weight)


#-----------------------------#
#    Plot the data            #
#-----------------------------#
ggplot(htwt, aes(x = height, y = weight)) +
  theme_bw() +
  geom_point(alpha = 0.5) +
  ylab("weight (in pounds)") +
  xlab("height (in inches)")


##############################
#  Frequentist approach      #
##############################

#--- Conduct the analyses ---#
model = lm(htwt$weight ~ htwt$height)
summary(model)

#--- Plot the results ---#
ggplot(htwt, aes(x = height, y = weight)) +
  theme_bw() +
  geom_point(alpha = 0.5) +
  geom_smooth(method = lm, colour = "dodgerblue4", alpha = 0.5) +
  ylab("weight (in pounds)") +
  xlab("height (in inches)")


################################
#      Bayesian Approach       #
################################

#--------------------------------------#
# Organize and Standardize the Data    #
#--------------------------------------#
y = htwt$weight
yM = mean(y)
ySD = sd(y)
zy = (y - yM) / ySD

x = htwt$height
xM = mean(x)
xSD = sd(x)
zx = (x - xM) / xSD

N = length(y)


#----------------------------#
# Plot the standardized data #
#----------------------------#
data = data.frame(zy, zx)
ggplot(data, aes(x = zx, y = zy)) +
  theme_bw() +
  geom_point(alpha = 0.5) +
  ylab("weight") +
  xlab("height")


#------------------------------#
# Save data as a list for STAN #
#------------------------------#
dataList = list (
  y = zy,
  x = zx,
  N = N
)


#------------------------------#
#       Define the Model       #
#------------------------------#
modelString = "
  data {
    int N;
    vector[N] y;
    vector[N] x;
  }

  parameters {
    real b0;
    real b1;
    real<lower=0> sigma;
  }

  model {
    // Definitions
    vector[N] mu;

    // Likelihood
    mu = b0 + (b1 * x);
    y ~ normal(mu, sigma);

    // Priors
    b0 ~ normal(0, 1);
    b1 ~ cauchy(0, 1);
    sigma ~ cauchy(1, 1);
  }

  generated quantities {
    vector[N] y_signal;
    vector[N] y_pred;
    
    for (n in 1:N) {
      y_signal[n] = b0 + (b1 * x[n]);
      y_pred[n] = normal_rng(b0 + (b1 * x[n]), sigma);
    }
  }
"
writeLines(modelString, con="model.stan")


#--------------------------#
#     run STAN             #
#--------------------------#
stanFit <- stan(file = "model.stan", 
                data = dataList, 
                pars = c("b0", "b1", "sigma", "y_signal", "y_pred"),
                warmup = 2000,
                iter = 7000, 
                chains = 3)


#-------------------------#
# Check MCMC Performance  #
#-------------------------#
print(stanFit)
stan_trace(stanFit, pars = c("b0", "b1", "sigma"), inc_warmup = TRUE)


#------------------------#
#    Plot Parameters     #
#------------------------#
stan_plot(stanFit, par = c("b0", "b1", "sigma"))
plot(stanFit, par = c("b0", "b1", "sigma"), show_density = TRUE, ci_level = 0.95, fill_color = "skyblue")


#-----------------------#
#    Extract Data       #
#-----------------------#
mcmcChains = as.data.frame(stanFit)
head(mcmcChains)

zb0 = mcmcChains$b0
zb1 = mcmcChains$b1
zsigma = mcmcChains$sigma
ysignal = mcmcChains[, 4:103]
ypred = mcmcChains[, 104:203]

#--------------------------------#
# Convert back to original scale #
#--------------------------------#
sigma = zsigma * ySD
b1 = zb1 * ySD / xSD
b0 = ((zb0 * ySD) + yM - (zb1 * ySD * xM / xSD))

#---------------------------#
#  Plot the data            #
#---------------------------#
source("plotPost.R")
par(mfrow = c(1, 3))
histInfo = plotPost(b0, xlab = bquote(beta[0]))
histInfo = plotPost(b1, xlab = bquote(beta[1]))
histInfo = plotPost(sigma, xlab = bquote(sigma))

histInfo = plotPost(b0, xlab = bquote(beta[0]), credMass = 0.89, col = "gray", showMode = TRUE)
histInfo = plotPost(b1, xlab = bquote(beta[1]), credMass = 0.89, col = "gray", showMode = TRUE)
histInfo = plotPost(sigma, xlab = bquote(sigma), credMass = 0.89, col = "gray", showMode = TRUE)

par(mfrow = c(1, 1))
histInfo = plotPost(b1, xlab = bquote(beta[1]), col = rgb(1, 0.2, 0.2, 0.5))


#-------------------------------------#
#    Assess fit of model              #
#-------------------------------------#
signalMeans = apply(ysignal, 2, mean)
signalHDIMin = apply(ysignal, 2, quantile, probs = 0.025)
signalHDIMax = apply(ysignal, 2, quantile, probs = 0.975)
signalSummary = data.frame(signalMeans, signalHDIMin, signalHDIMax, zx, zy)

ggplot(signalSummary) +
  theme_bw() +
  geom_point(aes(x = zx, y = zy), alpha = 0.5) +
  geom_point(aes(x = zx, y = signalMeans), colour = "dodgerblue4", alpha = 0.5) +
  ylab("weight") +
  xlab("height")

ggplot(signalSummary) +
  theme_bw() +
  geom_point(aes(x = zx, y = zy), alpha = 0.5) +
  geom_point(aes(x = zx, y = signalMeans), colour = "dodgerblue4", alpha = 0.5) +
  geom_ribbon(aes(x = zx, ymin = signalHDIMin, ymax = signalHDIMax), fill = "dodgerblue4", alpha = 0.3) +
  ylab("weight") +
  xlab("height")


#---------------------------------#
# Posterior Predictive Check      #
#---------------------------------#
predMeans = apply(ypred, 2, mean)
predHDIMin = apply(ypred, 2, quantile, probs = 0.025)
predHDIMax = apply(ypred, 2, quantile, probs = 0.975)
predictions = data.frame(predMeans, predHDIMin, predHDIMax, zx, zy)

ggplot(predictions) +
  theme_bw() +
  geom_point(aes(x = zx, y = zy), alpha = 0.5) +
  geom_point(aes(x = zx, y = predMeans), colour = "dodgerblue4", alpha = 0.5) +
  geom_ribbon(aes(x = zx, ymin = predHDIMin, ymax = predHDIMax), fill = "dodgerblue4", alpha = 0.3) +
  ylab("weight") +
  xlab("height")



#############################
# Predict weights for new   #
# height values.            #
#############################

#--- Generate new x (height) values ---#
xNew = runif(n = 50, min = 70, max = 90)
xNew

xNew = sort(xNew)
xNew

#--- Define a matrix to hold predicted y values ---#
postSampSize = length(b1)
yNew = matrix(0, nrow = length(xNew), ncol = postSampSize)

#--- Define a matrix to hold HDI limits of predicted values ---#
yHDIlim = matrix(0, nrow = length(xNew), ncol = 2)

#-----------------------------------#
# populate the yNew matrix by       #
# generating one predicted y value  #
# for each step in the chain        #
#-----------------------------------#
for (i in 1:length(xNew)) {
  for (j in 1:postSampSize) {
    yNew[i, j] = rnorm(n = 1, mean = b0[j] + b1[j] * xNew[i], sd = sigma[j])
  }
}

#--- Calculate Means and HDI for predicted values ---#
means = rowMeans(yNew)

for (i in 1:length(xNew)) {
  yHDIlim[i, ] = quantile(yNew[i, ], probs = c(0.025, 0.975))
}

#--- Combine the data ---$
predTable = data.frame(xNew, means, yHDIlim)

#--- Plot the data ---#
ggplot(predTable) +
  theme_bw() +
  geom_point(aes(x = xNew, y = means), alpha = 0.7) +
  geom_ribbon(aes(x = xNew, ymin = predTable[, 3], ymax = predTable[, 4]), fill = "dodgerblue4", alpha = 0.3) +
  ylab("weight") +
  xlab("height")
