lleg = rnorm(n = 100, mean = 32, sd = 6)
legnoise = rnorm(n = 100, mean = 0, sd = 0.5)
rleg = lleg + legnoise
noise = rnorm(n = 100, mean = 0, sd = 4)
height = (lleg * 2) + noise
legs = data.frame(lleg, rleg, height)
write.csv(legs, "legs.csv", row.names = FALSE, quote = FALSE)


#################################
# Code for analyzing leg length #
# data.                         #
#################################

#--------------------------#
# Load libraries and files #
#--------------------------#
library(rstan)
library(ggplot2)
source("plotPost.R")

legs = read.table("legs.csv", header = TRUE, sep = ",")


#---------------#
# Plot the data #
#---------------#
pairs(legs, pch = 16, col = rgb(0, 0, 1, 0.5))

#-----------------------#
# Check the correlation #
#-----------------------#
cor(legs)


#-------------------------------------------#
# Standardize the data and prepare for Stan #
#-------------------------------------------#

heightMean = mean(legs$height)
heightSD = sd(legs$height)
zheight = (legs$height - heightMean) / heightSD

N = length(zheight)

llegMean = mean(legs$lleg)
llegSD = sd(legs$lleg)
zlleg = (legs$lleg - llegMean) / llegSD

rlegMean = mean(legs$rleg)
rlegSD = sd(legs$rleg)
zrleg = (legs$rleg - rlegMean) / rlegSD


##############################
#   EACH LEG INDEPENDENTLY   #
##############################

#----------------------------#
#      LEFT LEG              #
#----------------------------#
dataList = list (
  height = zheight,
  N = N,
  leg = zlleg
)

#------------------------------#
#       Define the Model       #
#------------------------------#
modelString = "
  data {
    int N;              // Sample size
    vector[N] height;   // Vector of height data
    vector[N] leg;      // Vector of leg data
  }

  parameters {
    real b0;              // Intercept
    real b1;              // Coefficient for left leg effect
    real<lower=0> sigma;  // Standard deviation
  }

  model {
    // Definitions
    vector[N] mu;

    // Likelihood
    mu = b0 + (b1 * leg);
    height ~ normal(mu, sigma);

    // Priors
    b0 ~ normal(0, 1);
    b1 ~ normal(0, 1);
    sigma ~ cauchy(1, 1);
  }

  generated quantities {
    vector[N] y_pred;
    vector[N] mu_pred;

    for (i in 1:N) {
      mu_pred[i] = b0 + (b1 * leg[i]);
      y_pred[i] = normal_rng(mu_pred[i], sigma);
    }
  }
"
writeLines(modelString, con="model.stan")


#--------------------------#
#     run STAN             #
#--------------------------#
stanFit <- stan(file = "model.stan", 
                data = dataList, 
                pars = c("b0", "b1", "sigma", "y_pred"),
                warmup = 2000,
                iter = 7000, 
                chains = 3)


#-------------------------#
# Check MCMC Performance  #
#-------------------------#
print(stanFit)
stan_trace(stanFit, pars = c("b0", "b1"), inc_warmup = TRUE)
stan_trace(stanFit, pars = "sigma", inc_warmup = TRUE)


#-----------------------------------#
#  Plot Results With Stan Functions #
#-----------------------------------#
stan_plot(stanFit, par = c("b0", "b1"))
stan_plot(stanFit, par = "sigma")


#----------------------------#
#      RIGHT LEG             #
#----------------------------#
dataList = list (
  height = zheight,
  N = N,
  leg = zrleg
)


#--------------------------#
#     run STAN             #
#--------------------------#
stanFit <- stan(file = "model.stan", 
                data = dataList, 
                pars = c("b0", "b1", "sigma", "y_pred"),
                warmup = 2000,
                iter = 7000, 
                chains = 3)


#-------------------------#
# Check MCMC Performance  #
#-------------------------#
print(stanFit)
stan_trace(stanFit, pars = c("b0", "b1"), inc_warmup = TRUE)
stan_trace(stanFit, pars = "sigma", inc_warmup = TRUE)


#-----------------------------------#
#  Plot Results With Stan Functions #
#-----------------------------------#
stan_plot(stanFit, par = c("b0", "b1"))
stan_plot(stanFit, par = "sigma")


####################################
#     MODEL WITH BOTH LEGS         #
####################################
dataList = list (
  height = zheight,
  N = N,
  lleg = zlleg,
  rleg = zrleg
)


#------------------------------#
#       Define the Model       #
#------------------------------#
modelString = "
  data {
    int N;                // Sample size
    vector[N] height;     // Vector of height data
    vector[N] lleg;       // Vector of left leg data
    vector[N] rleg;       // Vector of right leg data
  }

  parameters {
    real b0;              // Intercept
    real b1;              // Coefficient for left leg effect
    real b2;              // Coefficient for right leg effect
    real<lower=0> sigma;  // Standard deviation
  }

  model {
    // Definitions
    vector[N] mu;

    // Likelihood
    mu = b0 + (b1 * lleg) + (b2 * rleg);
    height ~ normal(mu, sigma);

    // Priors
    b0 ~ normal(0, 1);
    b1 ~ normal(0, 1);
    b2 ~ normal(0, 1);
    sigma ~ cauchy(1, 1);
  }

  generated quantities {
    vector[N] y_pred;
    vector[N] mu_pred;

    for (i in 1:N) {
      mu_pred[i] = b0 + (b1 * lleg[i]) + (b2 * rleg[i]);
      y_pred[i] = normal_rng(mu_pred[i], sigma);
    }
  }
"
writeLines(modelString, con="model.stan")


#--------------------------#
#     run STAN             #
#--------------------------#
stanFit <- stan(file = "model.stan", 
                data = dataList, 
                pars = c("b0", "b1", "b2", "sigma", "y_pred"),
                warmup = 2000,
                iter = 7000, 
                chains = 3)


#-------------------------#
# Check MCMC Performance  #
#-------------------------#
print(stanFit)
stan_trace(stanFit, pars = c("b0", "b1", "b2"), inc_warmup = TRUE)
stan_trace(stanFit, pars = "sigma", inc_warmup = TRUE)


#-----------------------------------#
#  Plot Results With Stan Functions #
#-----------------------------------#
stan_plot(stanFit, par = c("b0", "b1", "b2"))
stan_plot(stanFit, par = "sigma")


#####################################
#        FREQUENTIST APPROACH       #
#####################################

model_both = lm(legs$height ~ legs$lleg + legs$rleg)
model_lleg = lm(legs$height ~ legs$lleg)
model_rleg = lm(legs$height ~ legs$rleg)
AIC(model_both, model_lleg, model_rleg)
