################################
#     mixedRegression.R        #
#------------------------------#
# Code for conducting a        #
# regression analysis with     #
# multiple predictor variables #
# of metric and nominal data   #
# types.                       #
################################


#----------------------------#
# Load appropriate libraries #
#----------------------------#
library(ggplot2)
library(rstan)
library(loo)
options(mc.cores = parallel::detectCores())
source("plotPost.R")

#----------------------#
# Read the data into R #
#----------------------#
milk = read.table("milk.csv", header = TRUE, sep = ",")

#--- Get just rows without missing data ---#
milkFull = milk[complete.cases(milk), ]


#------------------------------#
# Plot pairwise comparisons of #
# predictor variables.         #
#------------------------------#
pairs(milkFull[, c(1, 3:7)], pch = 16, col = rgb(0, 0, 1, 0.5))

#--- Correlation ---#
cor(milkFull[, 3:7])


#-----------------------------#
#     Prepare the data        #
#-----------------------------#

#--- Neocortext ---#
y = milkFull$neocortex.perc
yMean = mean(y)
ySD = sd(y)
zy = (y - yMean) / ySD
N = length(y)

#--- kcal.per.g ---#
x1 = milkFull$kcal.per.g
x1Mean = mean(x1)
x1SD = sd(x1)
zx1 = (x1 - x1Mean) / x1SD

#--- perc.fat ---#
x2 = milkFull$perc.fat
x2Mean = mean(x2)
x2SD = sd(x2)
zx2 = (x2 - x2Mean) / x2SD

#--- perc.protein ---#
x3 = milkFull$perc.protein
x3Mean = mean(x3)
x3SD = sd(x3)
zx3 = (x3 - x3Mean) / x3SD

#--- perc.lactose ---#
x4 = milkFull$perc.lactose
x4Mean = mean(x4)
x4SD = sd(x4)
zx4 = (x4 - x4Mean) / x4SD

#--- mass ---#
x5 = milkFull$mass
x5Mean = mean(x5)
x5SD = sd(x5)
zx5 = (x5 - x5Mean) / x5SD

#--- clade ---#
x6 = as.numeric(milkFull$clade)
x6Names = levels(milkFull$clade)
nx6Levels = length(unique(milkFull$clade))


#---------------------------#
# Prepare the data for Stan #
#---------------------------#
dataList = list(
  y = zy,
  x1 = zx1,
  x2 = zx2,
  x3 = zx3,
  x4 = zx4,
  x5 = zx5,
  x6 = x6,
  N = N,
  nx6Levels = nx6Levels
)

#------------------------------Model 1 = Full Model ---------------------------#

#-----------------------------#
# Define the model and write  #
# a string for Stan           #
#-----------------------------#
modelstring = "
  data {
    int N;          // Sample size
    int nx6Levels;  // Number of clade types (levels in our nominal variable)
    vector[N] y;    // Vector of neocortex percentages
    vector[N] x1;   // Vector of kcal.per.g data
    vector[N] x2;   // Vector of perc.fat data
    vector[N] x3;   // Vector of perc.protein data
    vector[N] x4;   // Vector of perc.lactose data
    vector[N] x5;   // Vector of mass data
    int x6[N];      // Vector of indicators of which clade each sample is in
  }

  parameters {
    real b0;              
    real b1;                    // Coefficient for effect of kcal.per.g
    real b2;                    // Coefficient for effect of perc.fat
    real b3;                    // Coefficient for effect of perc.protein
    real b4;                    // Coefficient for effect of perc.lactose
    real b5;                    // Coefficient for effect of mass
    real b6[nx6Levels];         // Coefficients for effects of being in each clade
    real<lower=0> sigma;        // Coefficients for overall sd
    real cladeMean;             // Mean effect across all clades
    real<lower=0> cladeMeanSD;  // sd for mean effect across all clades
  }

  model {
    // Definitions
    vector[N] mu;

    // Likelihood
    for (i in 1:N) {
      mu[i] = b0 + b1*x1[i] + b2*x2[i] + b3*x3[i] + b4*x4[i] + b5*x5[i] + b6[x6[i]];
      y[i] ~ normal(mu[i], sigma);
    }

    // Priors
    b0 ~ normal(0, 1);
    b1 ~ normal(0, 1);
    b2 ~ normal(0, 1);
    b3 ~ normal(0, 1);
    b4 ~ normal(0, 1);
    b5 ~ normal(0, 1);
    sigma ~ cauchy(0, 1);

    for (j in 1:nx6Levels) {
      b6[j] ~ normal(cladeMean, cladeMeanSD);
    }

    // Hyperpriors
    cladeMean   ~ normal(0, 1);
    cladeMeanSD ~ normal(1, 1);
  }

  generated quantities {
    // Posterior Predictive Variable Definitions
    vector[N] mu_pred;
    vector[N] y_pred;

    // WAIC Variable Definitions
    vector[N] log_lik;
    vector[N] mu_waic;

    // For Posterior Predictive Calculations
    for (i in 1:N) {
      mu_pred[i] = b0 + b1*x1[i] + b2*x2[i] + b3*x3[i] + b4*x4[i] + b5*x5[i] + b6[x6[i]];
      y_pred[i] = normal_rng(mu_pred[i], sigma);
    }

    // For WAIC Calculatons
    for (i in 1:N) {
      mu_waic[i] = b0 + b1*x1[i] + b2*x2[i] + b3*x3[i] + b4*x4[i] + b5*x5[i] + b6[x6[i]];
      log_lik[i] = normal_lpdf(y[i] | mu_waic[i], sigma);
    }
  }
"
writeLines(modelstring, con = "model_Full.stan")


#--------------------------#
#     run STAN             #
#--------------------------#
model_Full = stan(file = "model_Full.stan", 
              data = dataList, 
              pars = c("b0", "b1", "b2", "b3", "b4", "b5", "b6", "sigma", "y_pred", "log_lik"),
              warmup = 2000,
              iter = 12000, 
              chains = 3)


#-------------------------#
# Check MCMC Performance  #
#-------------------------#
print(model_Full)
stan_trace(model_Full, pars = c("b0", "b1"), inc_warmup = TRUE)
stan_trace(model_Full, pars = c("b2", "b3"), inc_warmup = TRUE)
stan_trace(model_Full, pars = c("b4", "b5"), inc_warmup = TRUE)
stan_trace(model_Full, pars = "b6", inc_warmup = TRUE)
stan_trace(model_Full, pars = "sigma", inc_warmup = TRUE)


#-----------------------------------#
#  Plot Results With Stan Functions #
#-----------------------------------#
stan_plot(model_Full, par = c("b0", "b1", "b2", "b3", "b4", "b5", "b6"))
stan_plot(model_Full, par = "sigma")


#-----------------------------------#
#  Posterior Predictive Check       #
#-----------------------------------#

#--- Extract the predictions ---#
mcmcChains = as.data.frame(model_Full)
chainLength = length(mcmcChains[, 1])

zypred = matrix(0, ncol = N, nrow = chainLength)
for (i in 1:N) {
  zypred[, i] = mcmcChains[, paste("y_pred[", i, "]", sep = "")]
}

#--- Mean expected value for record ---#
ypredMean = apply(zypred, 2, mean)

#--- Upper and lower expected 95% HDI for each visit ---#
ypredLow  = apply(zypred, 2, quantile, probs = 0.025)
ypredHigh = apply(zypred, 2, quantile, probs = 0.975)

#--- Plot mean predicted values ---#
record = 1:N
dotchart(ypredMean, xlim = c(-3, 3), xlab = "Standardized Percent Neocortex", ylab = "Record")

#--- Add HDIs ---#
segments(x0 = ypredLow, y0 = record, x1 = ypredHigh, y1 = record)

#--- Add observed values ---#
points(x = zy, y = record, pch = 16, col = rgb(1, 0, 0, 0.6))


#-----------------------------------#
#  Calculate WAIC                   #
#-----------------------------------#
loglik_Full = extract_log_lik(model_Full)
waic_Full = waic(loglik_Full)


#------------------------------Model 2 = Removing perc.lactose ---------------------------#

#-----------------------------#
# Define the model and write  #
# a string for Stan           #
#-----------------------------#
modelstring = "
  data {
    int N;          // Sample size
    int nx6Levels;  // Number of clade types (levels in our nominal variable)
    vector[N] y;    // Vector of neocortex percentages
    vector[N] x1;   // Vector of kcal.per.g data
    vector[N] x2;   // Vector of perc.fat data
    vector[N] x3;   // Vector of perc.protein data
    vector[N] x5;   // Vector of mass data
    int x6[N];      // Vector of indicators of which clade each sample is in
  }

  parameters {
    real b0;              
    real b1;                    // Coefficient for effect of kcal.per.g
    real b2;                    // Coefficient for effect of perc.fat
    real b3;                    // Coefficient for effect of perc.protein
    real b5;                    // Coefficient for effect of mass
    real b6[nx6Levels];         // Coefficients for effects of being in each clade
    real<lower=0> sigma;        // Coefficients for overall sd
    real cladeMean;             // Mean effect across all clades
    real<lower=0> cladeMeanSD;  // sd for mean effect across all clades
  }

  model {
    // Definitions
    vector[N] mu;

    // Likelihood
    for (i in 1:N) {
      mu[i] = b0 + b1*x1[i] + b2*x2[i] + b3*x3[i] + b5*x5[i] + b6[x6[i]];
      y[i] ~ normal(mu[i], sigma);
    }

    // Priors
    b0 ~ normal(0, 1);
    b1 ~ normal(0, 1);
    b2 ~ normal(0, 1);
    b3 ~ normal(0, 1);
    b5 ~ normal(0, 1);
    sigma ~ cauchy(0, 1);

    for (j in 1:nx6Levels) {
      b6[j] ~ normal(cladeMean, cladeMeanSD);
    }

    // Hyperpriors
    cladeMean   ~ normal(0, 1);
    cladeMeanSD ~ normal(1, 1);
  }

  generated quantities {
    // Posterior Predictive Variable Definitions
    vector[N] mu_pred;
    vector[N] y_pred;

    // WAIC Variable Definitions
    vector[N] log_lik;
    vector[N] mu_waic;

    // For Posterior Predictive Calculations
    for (i in 1:N) {
      mu_pred[i] = b0 + b1*x1[i] + b2*x2[i] + b3*x3[i] + b5*x5[i] + b6[x6[i]];
      y_pred[i] = normal_rng(mu_pred[i], sigma);
    }

    // For WAIC Calculatons
    for (i in 1:N) {
      mu_waic[i] = b0 + b1*x1[i] + b2*x2[i] + b3*x3[i] + b5*x5[i] + b6[x6[i]];
      log_lik[i] = normal_lpdf(y[i] | mu_waic[i], sigma);
    }
  }
"
writeLines(modelstring, con = "model_2.stan")


#--------------------------#
#     run STAN             #
#--------------------------#
model_2 = stan(file = "model_2.stan", 
                  data = dataList, 
                  pars = c("b0", "b1", "b2", "b3", "b5", "b6", "sigma", "y_pred", "log_lik"),
                  warmup = 2000,
                  iter = 12000, 
                  chains = 3)


#-------------------------#
# Check MCMC Performance  #
#-------------------------#
print(model_2)
stan_trace(model_2, pars = "b0", inc_warmup = TRUE)
stan_trace(model_2, pars = c("b1", "b2"), inc_warmup = TRUE)
stan_trace(model_2, pars = c("b3", "b5"), inc_warmup = TRUE)
stan_trace(model_2, pars = "b6", inc_warmup = TRUE)
stan_trace(model_2, pars = "sigma", inc_warmup = TRUE)


#-----------------------------------#
#  Plot Results With Stan Functions #
#-----------------------------------#
stan_plot(model_2, par = c("b0", "b1", "b2", "b3", "b5", "b6"))
stan_plot(model_2, par = "sigma")


#-----------------------------------#
#  Posterior Predictive Check       #
#-----------------------------------#

#--- Extract the predictions ---#
mcmcChains = as.data.frame(model_2)
chainLength = length(mcmcChains[, 1])

zypred = matrix(0, ncol = N, nrow = chainLength)
for (i in 1:N) {
  zypred[, i] = mcmcChains[, paste("y_pred[", i, "]", sep = "")]
}

#--- Mean expected value for record ---#
ypredMean = apply(zypred, 2, mean)

#--- Upper and lower expected 95% HDI for each visit ---#
ypredLow  = apply(zypred, 2, quantile, probs = 0.025)
ypredHigh = apply(zypred, 2, quantile, probs = 0.975)

#--- Plot mean predicted values ---#
record = 1:N
dotchart(ypredMean, xlim = c(-3, 3), xlab = "Standardized Percent Neocortex", ylab = "Record")

#--- Add HDIs ---#
segments(x0 = ypredLow, y0 = record, x1 = ypredHigh, y1 = record)

#--- Add observed values ---#
points(x = zy, y = record, pch = 16, col = rgb(1, 0, 0, 0.6))


#-----------------------------------#
#  Calculate WAIC                   #
#-----------------------------------#
loglik_2 = extract_log_lik(model_2)
waic_2 = waic(loglik_2)


#------------------------------Model 3 = Removing perc.lactose and kcal.per.g ---------------------------#

#-----------------------------#
# Define the model and write  #
# a string for Stan           #
#-----------------------------#
modelstring = "
  data {
    int N;          // Sample size
    int nx6Levels;  // Number of clade types (levels in our nominal variable)
    vector[N] y;    // Vector of neocortex percentages
    vector[N] x2;   // Vector of perc.fat data
    vector[N] x3;   // Vector of perc.protein data
    vector[N] x5;   // Vector of mass data
    int x6[N];      // Vector of indicators of which clade each sample is in
  }

  parameters {
    real b0;              
    real b2;                    // Coefficient for effect of perc.fat
    real b3;                    // Coefficient for effect of perc.protein
    real b5;                    // Coefficient for effect of mass
    real b6[nx6Levels];         // Coefficients for effects of being in each clade
    real<lower=0> sigma;        // Coefficients for overall sd
    real cladeMean;             // Mean effect across all clades
    real<lower=0> cladeMeanSD;  // sd for mean effect across all clades
  }

  model {
    // Definitions
    vector[N] mu;

    // Likelihood
    for (i in 1:N) {
      mu[i] = b0 + b2*x2[i] + b3*x3[i] + b5*x5[i] + b6[x6[i]];
      y[i] ~ normal(mu[i], sigma);
    }

    // Priors
    b0 ~ normal(0, 1);
    b2 ~ normal(0, 1);
    b3 ~ normal(0, 1);
    b5 ~ normal(0, 1);
    sigma ~ cauchy(0, 1);

    for (j in 1:nx6Levels) {
      b6[j] ~ normal(cladeMean, cladeMeanSD);
    }

    // Hyperpriors
    cladeMean   ~ normal(0, 1);
    cladeMeanSD ~ normal(1, 1);
  }

  generated quantities {
    // Posterior Predictive Variable Definitions
    vector[N] mu_pred;
    vector[N] y_pred;

    // WAIC Variable Definitions
    vector[N] log_lik;
    vector[N] mu_waic;

    // For Posterior Predictive Calculations
    for (i in 1:N) {
      mu_pred[i] = b0 + b2*x2[i] + b3*x3[i] + b5*x5[i] + b6[x6[i]];
      y_pred[i] = normal_rng(mu_pred[i], sigma);
    }

    // For WAIC Calculatons
    for (i in 1:N) {
      mu_waic[i] = b0 + b2*x2[i] + b3*x3[i] + b5*x5[i] + b6[x6[i]];
      log_lik[i] = normal_lpdf(y[i] | mu_waic[i], sigma);
    }
  }
"
writeLines(modelstring, con = "model_3.stan")


#--------------------------#
#     run STAN             #
#--------------------------#
model_3 = stan(file = "model_3.stan", 
               data = dataList, 
               pars = c("b0", "b2", "b3", "b5", "b6", "sigma", "y_pred", "log_lik"),
               warmup = 2000,
               iter = 12000, 
               chains = 3)


#-------------------------#
# Check MCMC Performance  #
#-------------------------#
print(model_3)
stan_trace(model_3, pars = c("b0", "b2"), inc_warmup = TRUE)
stan_trace(model_3, pars = c("b3", "b5"), inc_warmup = TRUE)
stan_trace(model_3, pars = "b6", inc_warmup = TRUE)
stan_trace(model_3, pars = "sigma", inc_warmup = TRUE)


#-----------------------------------#
#  Plot Results With Stan Functions #
#-----------------------------------#
stan_plot(model_3, par = c("b0", "b2", "b3", "b5", "b6"))
stan_plot(model_3, par = "sigma")


#-----------------------------------#
#  Posterior Predictive Check       #
#-----------------------------------#

#--- Extract the predictions ---#
mcmcChains = as.data.frame(model_3)
chainLength = length(mcmcChains[, 1])

zypred = matrix(0, ncol = N, nrow = chainLength)
for (i in 1:N) {
  zypred[, i] = mcmcChains[, paste("y_pred[", i, "]", sep = "")]
}

#--- Mean expected value for record ---#
ypredMean = apply(zypred, 2, mean)

#--- Upper and lower expected 95% HDI for each visit ---#
ypredLow  = apply(zypred, 2, quantile, probs = 0.025)
ypredHigh = apply(zypred, 2, quantile, probs = 0.975)

#--- Plot mean predicted values ---#
record = 1:N
dotchart(ypredMean, xlim = c(-3, 3), xlab = "Standardized Percent Neocortex", ylab = "Record")

#--- Add HDIs ---#
segments(x0 = ypredLow, y0 = record, x1 = ypredHigh, y1 = record)

#--- Add observed values ---#
points(x = zy, y = record, pch = 16, col = rgb(1, 0, 0, 0.6))


#-----------------------------------#
#  Calculate WAIC                   #
#-----------------------------------#
loglik_3 = extract_log_lik(model_3)
waic_3 = waic(loglik_3)
compare(waic_Full, waic_2, waic_3)


#------------------------------Model 4 = Removing perc.lactose and perc.fat ---------------------------#

#-----------------------------#
# Define the model and write  #
# a string for Stan           #
#-----------------------------#
modelstring = "
  data {
    int N;          // Sample size
    int nx6Levels;  // Number of clade types (levels in our nominal variable)
    vector[N] y;    // Vector of neocortex percentages
    vector[N] x1;   // Vector of kcal.per.g data
    vector[N] x3;   // Vector of perc.protein data
    vector[N] x5;   // Vector of mass data
    int x6[N];      // Vector of indicators of which clade each sample is in
  }

  parameters {
    real b0;      
    real b1;                    // Coefficient for effect of kcal.per.g
    real b3;                    // Coefficient for effect of perc.protein
    real b5;                    // Coefficient for effect of mass
    real b6[nx6Levels];         // Coefficients for effects of being in each clade
    real<lower=0> sigma;        // Coefficients for overall sd
    real cladeMean;             // Mean effect across all clades
    real<lower=0> cladeMeanSD;  // sd for mean effect across all clades
  }

  model {
    // Definitions
    vector[N] mu;

    // Likelihood
    for (i in 1:N) {
      mu[i] = b0 + b1*x1[i] + b3*x3[i] + b5*x5[i] + b6[x6[i]];
      y[i] ~ normal(mu[i], sigma);
    }

    // Priors
    b0 ~ normal(0, 1);
    b1 ~ normal(0, 1);
    b3 ~ normal(0, 1);
    b5 ~ normal(0, 1);
    sigma ~ cauchy(0, 1);

    for (j in 1:nx6Levels) {
      b6[j] ~ normal(cladeMean, cladeMeanSD);
    }

    // Hyperpriors
    cladeMean   ~ normal(0, 1);
    cladeMeanSD ~ normal(1, 1);
  }

  generated quantities {
    // Posterior Predictive Variable Definitions
    vector[N] mu_pred;
    vector[N] y_pred;

    // WAIC Variable Definitions
    vector[N] log_lik;
    vector[N] mu_waic;

    // For Posterior Predictive Calculations
    for (i in 1:N) {
      mu_pred[i] = b0 + b1*x1[i] + b3*x3[i] + b5*x5[i] + b6[x6[i]];
      y_pred[i] = normal_rng(mu_pred[i], sigma);
    }

    // For WAIC Calculatons
    for (i in 1:N) {
      mu_waic[i] = b0 + b1*x1[i] + b3*x3[i] + b5*x5[i] + b6[x6[i]];
      log_lik[i] = normal_lpdf(y[i] | mu_waic[i], sigma);
    }
  }
"
writeLines(modelstring, con = "model_4.stan")


#--------------------------#
#     run STAN             #
#--------------------------#
model_4 = stan(file = "model_4.stan", 
               data = dataList, 
               pars = c("b0", "b1", "b3", "b5", "b6", "sigma", "y_pred", "log_lik"),
               warmup = 2000,
               iter = 12000, 
               chains = 3)


#-------------------------#
# Check MCMC Performance  #
#-------------------------#
print(model_4)
stan_trace(model_4, pars = c("b0", "b1"), inc_warmup = TRUE)
stan_trace(model_4, pars = c("b3", "b5"), inc_warmup = TRUE)
stan_trace(model_4, pars = "b6", inc_warmup = TRUE)
stan_trace(model_4, pars = "sigma", inc_warmup = TRUE)


#-----------------------------------#
#  Plot Results With Stan Functions #
#-----------------------------------#
stan_plot(model_4, par = c("b0", "b1", "b3", "b5", "b6"))
stan_plot(model_4, par = "sigma")


#-----------------------------------#
#  Posterior Predictive Check       #
#-----------------------------------#

#--- Extract the predictions ---#
mcmcChains = as.data.frame(model_4)
chainLength = length(mcmcChains[, 1])

zypred = matrix(0, ncol = N, nrow = chainLength)
for (i in 1:N) {
  zypred[, i] = mcmcChains[, paste("y_pred[", i, "]", sep = "")]
}

#--- Mean expected value for record ---#
ypredMean = apply(zypred, 2, mean)

#--- Upper and lower expected 95% HDI for each visit ---#
ypredLow  = apply(zypred, 2, quantile, probs = 0.025)
ypredHigh = apply(zypred, 2, quantile, probs = 0.975)

#--- Plot mean predicted values ---#
record = 1:N
dotchart(ypredMean, xlim = c(-3, 3), xlab = "Standardized Percent Neocortex", ylab = "Record")

#--- Add HDIs ---#
segments(x0 = ypredLow, y0 = record, x1 = ypredHigh, y1 = record)

#--- Add observed values ---#
points(x = zy, y = record, pch = 16, col = rgb(1, 0, 0, 0.6))


#-----------------------------------#
#  Calculate WAIC                   #
#-----------------------------------#
loglik_4 = extract_log_lik(model_4)
waic_4 = waic(loglik_4)
compare(waic_Full, waic_2, waic_3, waic_4)


#------------------------------Model 5 = Removing kcal.per.g and perc.fat ---------------------------#
#-----------------------------#
# Define the model and write  #
# a string for Stan           #
#-----------------------------#
modelstring = "
  data {
    int N;          // Sample size
    int nx6Levels;  // Number of clade types (levels in our nominal variable)
    vector[N] y;    // Vector of neocortex percentages
    vector[N] x3;   // Vector of perc.protein data
    vector[N] x4;   // Vector of perc.lactose data
    vector[N] x5;   // Vector of mass data
    int x6[N];      // Vector of indicators of which clade each sample is in
  }

  parameters {
    real b0;              
    real b3;                    // Coefficient for effect of perc.protein
    real b4;                    // Coefficient for effect of perc.lactose
    real b5;                    // Coefficient for effect of mass
    real b6[nx6Levels];         // Coefficients for effects of being in each clade
    real<lower=0> sigma;        // Coefficients for overall sd
    real cladeMean;             // Mean effect across all clades
    real<lower=0> cladeMeanSD;  // sd for mean effect across all clades
  }

  model {
    // Definitions
    vector[N] mu;

    // Likelihood
    for (i in 1:N) {
      mu[i] = b0 + b3*x3[i] + b4*x4[i] + b5*x5[i] + b6[x6[i]];
      y[i] ~ normal(mu[i], sigma);
    }

    // Priors
    b0 ~ normal(0, 1);
    b3 ~ normal(0, 1);
    b4 ~ normal(0, 1);
    b5 ~ normal(0, 1);
    sigma ~ cauchy(0, 1);

    for (j in 1:nx6Levels) {
      b6[j] ~ normal(cladeMean, cladeMeanSD);
    }

    // Hyperpriors
    cladeMean   ~ normal(0, 1);
    cladeMeanSD ~ normal(1, 1);
  }

  generated quantities {
    // Posterior Predictive Variable Definitions
    vector[N] mu_pred;
    vector[N] y_pred;

    // WAIC Variable Definitions
    vector[N] log_lik;
    vector[N] mu_waic;

    // For Posterior Predictive Calculations
    for (i in 1:N) {
      mu_pred[i] = b0 + b3*x3[i] + b4*x4[i] + b5*x5[i] + b6[x6[i]];
      y_pred[i] = normal_rng(mu_pred[i], sigma);
    }

    // For WAIC Calculatons
    for (i in 1:N) {
      mu_waic[i] = b0 + b3*x3[i] + b4*x4[i] + b5*x5[i] + b6[x6[i]];
      log_lik[i] = normal_lpdf(y[i] | mu_waic[i], sigma);
    }
  }
"
writeLines(modelstring, con = "model_5.stan")


#--------------------------#
#     run STAN             #
#--------------------------#
model_5 = stan(file = "model_5.stan", 
               data = dataList, 
               pars = c("b0", "b3", "b4", "b5", "b6", "sigma", "y_pred", "log_lik"),
               warmup = 2000,
               iter = 12000, 
               chains = 3)


#-------------------------#
# Check MCMC Performance  #
#-------------------------#
print(model_5)
stan_trace(model_5, pars = c("b0", "b3"), inc_warmup = TRUE)
stan_trace(model_5, pars = c("b4", "b5"), inc_warmup = TRUE)
stan_trace(model_5, pars = "b6", inc_warmup = TRUE)
stan_trace(model_5, pars = "sigma", inc_warmup = TRUE)


#-----------------------------------#
#  Plot Results With Stan Functions #
#-----------------------------------#
stan_plot(model_5, par = c("b0", "b3", "b4", "b5", "b6"))
stan_plot(model_5, par = "sigma")


#-----------------------------------#
#  Posterior Predictive Check       #
#-----------------------------------#

#--- Extract the predictions ---#
mcmcChains = as.data.frame(model_5)
chainLength = length(mcmcChains[, 1])

zypred = matrix(0, ncol = N, nrow = chainLength)
for (i in 1:N) {
  zypred[, i] = mcmcChains[, paste("y_pred[", i, "]", sep = "")]
}

#--- Mean expected value for record ---#
ypredMean = apply(zypred, 2, mean)

#--- Upper and lower expected 95% HDI for each visit ---#
ypredLow  = apply(zypred, 2, quantile, probs = 0.025)
ypredHigh = apply(zypred, 2, quantile, probs = 0.975)

#--- Plot mean predicted values ---#
record = 1:N
dotchart(ypredMean, xlim = c(-3, 3), xlab = "Standardized Percent Neocortex", ylab = "Record")

#--- Add HDIs ---#
segments(x0 = ypredLow, y0 = record, x1 = ypredHigh, y1 = record)

#--- Add observed values ---#
points(x = zy, y = record, pch = 16, col = rgb(1, 0, 0, 0.6))


#-----------------------------------#
#  Calculate WAIC                   #
#-----------------------------------#
loglik_5 = extract_log_lik(model_5)
waic_5 = waic(loglik_5)
compare(waic_Full, waic_2, waic_3, waic_4, waic_5)

#------------------------------------ INTERPRETATIONS ----------------------------------------------#

####################
#    kcal.per.g    #
#     model 4      #
####################

#--- Obtain Coefficients from Model ---#
mcmcData = as.data.frame(model_4)
zb0 = mcmcData$b0
zb1 = mcmcData$b1
zb3 = mcmcData$b3
zb5 = mcmcData$b5
zb6 = matrix(0, ncol = nx6Levels, nrow = chainLength)
for (i in 1:nx6Levels) {
  zb6[, i] = mcmcChains[, paste("b6[", i, "]", sep = "")]
}

#--------------------------------------#
#  Within Minimums of Other Variables  #
#--------------------------------------#

#--- Calculate predicted y values based on minimums ---#
kcalpergMins = matrix(0, nrow = length(zb0), ncol = length(zy))

for (i in 1:length(zb0)) {
  for (j in 1:length(zy)) {
    kcalpergMins[i, j] = zb0[i] + (zb1[i] * zx1[j]) + (zb3[i] * min(zx3)) + (zb5[i] * min(zx5) + zb6[i, 4]) 
  }
}

#--- Calculate mean and HDI values ---#
kcalpergMinsMean = apply(kcalpergMins, 2, mean)
kcalpergMinsLow = apply(kcalpergMins, 2, quantile, probs = 0.025)
kcalpergMinsHigh = apply(kcalpergMins, 2, quantile, probs = 0.975)

#--- Combine data into data frame and plot ---#
kcalpergMinsSummary = data.frame(zx1, kcalpergMinsMean, kcalpergMinsLow, kcalpergMinsHigh)

ggplot(kcalpergMinsSummary) +
  theme_bw() +
  geom_point(aes(x = zx1, y = kcalpergMinsMean), alpha = 0.6) +
  geom_ribbon(aes(x = zx1, ymin = kcalpergMinsLow, ymax = kcalpergMinsHigh), fill = "dodgerblue4", alpha = 0.3)


#--------------------------------------#
#  Within Means of Other Variables     #
#--------------------------------------#

#--- Calculate predicted y values based on means ---#
kcalpergMeans = matrix(0, nrow = length(zb0), ncol = length(zy))

for (i in 1:length(zb0)) {
  for (j in 1:length(zy)) {
    kcalpergMeans[i, j] = zb0[i] + (zb1[i] * zx1[j]) + (zb3[i] * mean(zx3)) + (zb5[i] * mean(zx5) + zb6[i, 2]) 
  }
}

#--- Calculate mean and HDI values ---#
kcalpergMeansMean = apply(kcalpergMeans, 2, mean)
kcalpergMeansLow = apply(kcalpergMeans, 2, quantile, probs = 0.025)
kcalpergMeansHigh = apply(kcalpergMeans, 2, quantile, probs = 0.975)

#--- Combine data into data frame and plot ---#
kcalpergMeansSummary = data.frame(zx1, kcalpergMeansMean, kcalpergMeansLow, kcalpergMeansHigh)

ggplot(kcalpergMeansSummary) +
  theme_bw() +
  geom_point(aes(x = zx1, y = kcalpergMeansMean), alpha = 0.6) +
  geom_ribbon(aes(x = zx1, ymin = kcalpergMeansLow, ymax = kcalpergMeansHigh), fill = "dodgerblue4", alpha = 0.3)


#--------------------------------------#
#  Within Max of Other Variables       #
#--------------------------------------#

#--- Calculate predicted y values based on means ---#
kcalpergMax = matrix(0, nrow = length(zb0), ncol = length(zy))

for (i in 1:length(zb0)) {
  for (j in 1:length(zy)) {
    kcalpergMax[i, j] = zb0[i] + (zb1[i] * zx1[j]) + (zb3[i] * max(zx3)) + (zb5[i] * max(zx5) + zb6[i, 3]) 
  }
}

#--- Calculate mean and HDI values ---#
kcalpergMaxMean = apply(kcalpergMax, 2, mean)
kcalpergMaxLow = apply(kcalpergMax, 2, quantile, probs = 0.025)
kcalpergMaxHigh = apply(kcalpergMax, 2, quantile, probs = 0.975)

#--- Combine data into data frame and plot ---#
kcalpergMaxSummary = data.frame(zx1, kcalpergMaxMean, kcalpergMaxLow, kcalpergMaxHigh)

ggplot(kcalpergMaxSummary) +
  theme_bw() +
  geom_point(aes(x = zx1, y = kcalpergMaxMean), alpha = 0.6) +
  geom_ribbon(aes(x = zx1, ymin = kcalpergMaxLow, ymax = kcalpergMaxHigh), fill = "dodgerblue4", alpha = 0.3)
