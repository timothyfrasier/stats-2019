###############################
# Load Appropariate Libraries #
###############################
library(rstan)
library(ggplot2)
library(loo)


##############################
#  Load the data and examine #
##############################
sat = read.table("Guber1999data.csv", header = TRUE, sep = ",")

#--- Plot the data ---#
pairs(sat[, 2:8], pch = 16, col = rgb(0, 0, 1, 0.5))

#--- Check Correlations Among Predictors ---#
cor(sat[, 2:5])


#############################
#       FULL MODEL          #
#############################

#---------------------#
# Organize the Data   #
#---------------------#

#--- SAT Scores ---#
y     = sat$SATT
yMean = mean(y)
ySD   = sd(y)
zy    = (y - yMean) / ySD

N = length(y)

#--- Spending ---#
x1     = sat$Spend
x1Mean = mean(x1)
x1SD   = sd(x1)
zx1    = (x1 - x1Mean) / x1SD

#--- Student:Teacher Ratio ---#
x2     = sat$StuTeaRat
x2Mean = mean(x2)
x2SD   = sd(x2)
zx2    = (x2 - x2Mean) / x2SD

#--- Teacher Salary ---#
x3     = sat$Salary
x3Mean = mean(x3)
x3SD   = sd(x3)
zx3    = (x3 - x3Mean) / x3SD 

#--- Precent Students Taking SAT ---#
x4     = sat$PrcntTake
x4Mean = mean(x4)
x4SD   = sd(x4)
zx4    = (x4 - x4Mean) / x4SD


#-----------------------------#
# Prepare the data for Stan   #
#-----------------------------#
dataList = list(
  y  = zy,
  x1 = zx1,
  x2 = zx2,
  x3 = zx3,
  x4 = zx4,
  N  = N
)

#-----------------------------#
# Define the model and write  #
# a string for Stan           #
#-----------------------------#
modelstring = "
  data {
    int N;          // Sample size
    vector[N] y;    // Vector of SAT scores
    vector[N] x1;   // Vector of Spending values
    vector[N] x2;   // Vector of Student:Teacher ratios
    vector[N] x3;   // Vector of Teacher Salaries
    vector[N] x4;   // Vector of percentage of students taking exams
  }

  parameters {
    real b0;                // Coefficient for 'intercept'
    real b1;                // Coefficient for effect of spending
    real b2;                // Coefficient for effect of student:teacher ratios
    real b3;                // Coefficient for effect of teacher salaries
    real b4;                // Coefficient for effect of % students taking SAT
    real<lower=0> sigma;    // Coefficient for sd for unexplained variation
  }

  model {
    // Definitions
    vector[N] mu;

    // Likelihood
    mu = b0 + b1*x1 + b2*x2 + b3*x3 + b4*x4;
    y ~ normal(mu, sigma);

    // Priors
    b0 ~ normal(0, 1);
    b1 ~ normal(0, 1);
    b2 ~ normal(0, 1);
    b3 ~ normal(0, 1);
    b4 ~ normal(0, 1);
    sigma ~ cauchy(1, 1);
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
      mu_pred[i] = b0 + b1*x1[i] + b2*x2[i] + b3*x3[i] + b4*x4[i];
      y_pred[i] = normal_rng(mu_pred[i], sigma);
    }

    // For WAIC Calculatons
    for (i in 1:N) {
      mu_waic[i] = b0 + b1*x1[i] + b2*x2[i] + b3*x3[i] + b4*x4[i];
      log_lik[i] = normal_lpdf(y[i] | mu_waic[i], sigma);
    }
  }
"
writeLines(modelstring, con = "model1.stan")


#--------------------------#
#     run STAN             #
#--------------------------#
model1 = stan(file = "model1.stan", 
              data = dataList, 
              pars = c("b0", "b1", "b2", "b3", "b4", "sigma", "y_pred", "log_lik"),
              warmup = 2000,
              iter = 7000, 
              chains = 3)


#-------------------------#
# Check MCMC Performance  #
#-------------------------#
print(model1)
stan_trace(model1, pars = "b0", inc_warmup = TRUE)
stan_trace(model1, pars = c("b1", "b2"), inc_warmup = TRUE)
stan_trace(model1, pars = c("b3", "b4"), inc_warmup = TRUE)
stan_trace(model1, pars = "sigma", inc_warmup = TRUE)


#-----------------------------------#
#  Plot Results With Stan Functions #
#-----------------------------------#
stan_plot(model1, par = c("b0", "b1", "b2", "b3", "b4"))
stan_plot(model1, par = "sigma")


#----------------#
#    WAIC        #
#----------------#
loglik1 = extract_log_lik(model1)
waic1 = waic(loglik1)


###############################
# MODEL WITH SPENDING REMOVED #
###############################

#-----------------------------#
# Define the model and write  #
# a string for Stan           #
#-----------------------------#
modelstring = "
  data {
    int N;          // Sample size
    vector[N] y;    // Vector of SAT scores
    vector[N] x2;   // Vector of Student:Teacher ratios
    vector[N] x3;   // Vector of Teacher Salaries
    vector[N] x4;   // Vector of percentage of students taking exams
  }

  parameters {
    real b0;                // Coefficient for 'intercept'
    real b2;                // Coefficient for effect of student:teacher ratios
    real b3;                // Coefficient for effect of teacher salaries
    real b4;                // Coefficient for effect of % students taking SAT
    real<lower=0> sigma;    // Coefficient for sd for 'noise' after accounting for predictor variables
  }

  model {
    // Definitions
    vector[N] mu;

    // Likelihood
    mu = b0 + b2*x2 + b3*x3 + b4*x4;
    y ~ normal(mu, sigma);

    // Priors
    b0 ~ normal(0, 1);
    b2 ~ normal(0, 1);
    b3 ~ normal(0, 1);
    b4 ~ normal(0, 1);
    sigma ~ cauchy(1, 1);
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
      mu_pred[i] = b0 + b2*x2[i] + b3*x3[i] + b4*x4[i];
      y_pred[i] = normal_rng(mu_pred[i], sigma);
    }

    // For WAIC Calculatons
    for (i in 1:N) {
      mu_waic[i] = b0 + b2*x2[i] + b3*x3[i] + b4*x4[i];
      log_lik[i] = normal_lpdf(y[i] | mu_waic[i], sigma);
    }
  }
"
writeLines(modelstring, con = "model2.stan")


#--------------------------#
#     run STAN             #
#--------------------------#
model2 = stan(file = "model2.stan", 
               data = dataList, 
               pars = c("b0", "b2", "b3", "b4", "sigma", "y_pred", "log_lik"),
               warmup = 2000,
               iter = 7000, 
               chains = 3)


#-------------------------#
# Check MCMC Performance  #
#-------------------------#
print(model2)
stan_trace(model2, pars = c("b0", "b2"), inc_warmup = TRUE)
stan_trace(model2, pars = c("b3", "b4"), inc_warmup = TRUE)
stan_trace(model2, pars = "sigma", inc_warmup = TRUE)


#-----------------------------------#
#  Plot Results With Stan Functions #
#-----------------------------------#
stan_plot(model2, par = c("b0", "b2", "b3", "b4"))
stan_plot(model2, par = "sigma")


#----------------#
#    WAIC        #
#----------------#
loglik2 = extract_log_lik(model2)
waic2 = waic(loglik2)


###############################
# MODEL WITH SALARY REMOVED #
###############################

#-----------------------------#
# Define the model and write  #
# a string for Stan           #
#-----------------------------#
modelstring = "
  data {
    int N;          // Sample size
    vector[N] y;    // Vector of SAT scores
    vector[N] x1;   // Vector of Spending values
    vector[N] x2;   // Vector of Student:Teacher ratios
    vector[N] x4;   // Vector of percentage of students taking exams
  }

  parameters {
    real b0;                // Coefficient for 'intercept'
    real b1;                // Coefficient for effect of spending
    real b2;                // Coefficient for effect of student:teacher ratios
    real b4;                // Coefficient for effect of % students taking SAT
    real<lower=0> sigma;    // Coefficient for sd for 'noise' after accounting for predictor variables
  }

  model {
    // Definitions
    vector[N] mu;

    // Likelihood
    mu = b0 + b1*x1 + b2*x2 + b4*x4;
    y ~ normal(mu, sigma);

    // Priors
    b0 ~ normal(0, 1);
    b1 ~ normal(0, 1);
    b2 ~ normal(0, 1);
    b4 ~ normal(0, 1);
    sigma ~ cauchy(1, 1);
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
      mu_pred[i] = b0 + b1*x1[i] + b2*x2[i] + b4*x4[i];
      y_pred[i] = normal_rng(mu_pred[i], sigma);
    }

    // For WAIC Calculatons
    for (i in 1:N) {
      mu_waic[i] = b0 + b1*x1[i] + b2*x2[i] + b4*x4[i];
      log_lik[i] = normal_lpdf(y[i] | mu_waic[i], sigma);
    }
  }
"
writeLines(modelstring, con = "model3.stan")


#--------------------------#
#     run STAN             #
#--------------------------#
model3 = stan(file = "model3.stan", 
               data = dataList, 
               pars = c("b0", "b1", "b2", "b4", "sigma", "y_pred", "log_lik"),
               warmup = 2000,
               iter = 7000, 
               chains = 3)


#-------------------------#
# Check MCMC Performance  #
#-------------------------#
print(model3)
stan_trace(model3, pars = c("b0", "b1"), inc_warmup = TRUE)
stan_trace(model3, pars = c("b2", "b4"), inc_warmup = TRUE)
stan_trace(model3, pars = "sigma", inc_warmup = TRUE)


#-----------------------------------#
#  Plot Results With Stan Functions #
#-----------------------------------#
stan_plot(model3, par = c("b0", "b1", "b2", "b4"))
stan_plot(model3, par = "sigma")


#----------------#
#    WAIC        #
#----------------#
loglik3 = extract_log_lik(model3)
waic3 = waic(loglik3)
compare(waic1, waic2, waic3)


###############################
#   INTERPRETING EFFECTS      #
###############################

###############################
#         Spending            #
#   Will base on model #3     #
###############################

#--------------------------------------#
#  Within Minimums of Other Variables  #
#--------------------------------------#

#--- Obtain Coefficients from Model ---#
mcmcData = as.data.frame(model3)
zb0 = mcmcData$b0
zb1 = mcmcData$b1
zb2 = mcmcData$b2
zb4 = mcmcData$b4

#--- Calculate predicted y values based on minimums ---#
spendingMins = matrix(0, nrow = length(zb0), ncol = length(zy))

for (i in 1:length(zb0)) {
  for (j in 1:length(zy)) {
    spendingMins[i, j] = zb0[i] + (zb1[i] * zx1[j]) + (zb2[i] * min(zx2)) + (zb4[i] * min(zx4)) 
  }
}

#--- Calculate mean and HDI values ---#
spendingMinsMean = apply(spendingMins, 2, mean)
spendingMinsLow = apply(spendingMins, 2, quantile, probs = 0.025)
spendingMinsHigh = apply(spendingMins, 2, quantile, probs = 0.975)

#--- Combine data into data frame and plot ---#
spendingMinsSummary = data.frame(zx1, spendingMinsMean, spendingMinsLow, spendingMinsHigh)

ggplot(spendingMinsSummary) +
  theme_bw() +
  geom_point(aes(x = zx1, y = spendingMinsMean), alpha = 0.6) +
  geom_ribbon(aes(x = zx1, ymin = spendingMinsLow, ymax = spendingMinsHigh), fill = "dodgerblue4", alpha = 0.3)


#--------------------------------------#
#  Within Means of Other Variables     #
#--------------------------------------#

#--- Calculate predicted y values based on means ---#
spendingMeans = matrix(0, nrow = length(zb0), ncol = length(zy))

for (i in 1:length(zb0)) {
  for (j in 1:length(zy)) {
    spendingMeans[i, j] = zb0[i] + (zb1[i] * zx1[j]) + (zb2[i] * mean(zx2)) + (zb4[i] * mean(zx4)) 
  }
}

#--- Calculate mean and HDI values ---#
spendingMeansMean = apply(spendingMeans, 2, mean)
spendingMeansLow = apply(spendingMeans, 2, quantile, probs = 0.025)
spendingMeansHigh = apply(spendingMeans, 2, quantile, probs = 0.975)

#--- Combine data into data frame and plot ---#
spendingMeansSummary = data.frame(zx1, spendingMeansMean, spendingMeansLow, spendingMeansHigh)

ggplot(spendingMeansSummary) +
  theme_bw() +
  geom_point(aes(x = zx1, y = spendingMeansMean), alpha = 0.6) +
  geom_ribbon(aes(x = zx1, ymin = spendingMeansLow, ymax = spendingMeansHigh), fill = "dodgerblue4", alpha = 0.3)


#--------------------------------------#
#    Within Max of Other Variables     #
#--------------------------------------#

#--- Calculate predicted y values based on means ---#
spendingMax = matrix(0, nrow = length(zb0), ncol = length(zy))

for (i in 1:length(zb0)) {
  for (j in 1:length(zy)) {
    spendingMax[i, j] = zb0[i] + (zb1[i] * zx1[j]) + (zb2[i] * max(zx2)) + (zb4[i] * max(zx4)) 
  }
}

#--- Calculate mean and HDI values ---#
spendingMaxMean = apply(spendingMax, 2, mean)
spendingMaxLow = apply(spendingMax, 2, quantile, probs = 0.025)
spendingMaxHigh = apply(spendingMax, 2, quantile, probs = 0.975)

#--- Combine data into data frame and plot ---#
spendingMaxSummary = data.frame(zx1, spendingMaxMean, spendingMaxLow, spendingMaxHigh)

ggplot(spendingMaxSummary) +
  theme_bw() +
  geom_point(aes(x = zx1, y = spendingMaxMean), alpha = 0.6) +
  geom_ribbon(aes(x = zx1, ymin = spendingMaxLow, ymax = spendingMaxHigh), fill = "dodgerblue4", alpha = 0.3)


###############################
#         Salary              #
#   Will base on model #2     #
###############################

#--------------------------------------#
#  Within Minimums of Other Variables  #
#--------------------------------------#

#--- Obtain Coefficients from Model ---#
mcmcData = as.data.frame(model2)
zb0 = mcmcData$b0
zb1 = mcmcData$b1
zb2 = mcmcData$b2
zb3 = mcmcData$b3
zb4 = mcmcData$b4

#--- Calculate predicted y values based on minimums ---#
salaryMins = matrix(0, nrow = length(zb0), ncol = length(zy))

for (i in 1:length(zb0)) {
  for (j in 1:length(zy)) {
    salaryMins[i, j] = zb0[i] + (zb2[i] * min(zx2)) + (zb3[i] * zx3[j]) + (zb4[i] * min(zx4)) 
  }
}

#--- Calculate mean and HDI values ---#
salaryMinsMean = apply(salaryMins, 2, mean)
salaryMinsLow = apply(salaryMins, 2, quantile, probs = 0.025)
salaryMinsHigh = apply(salaryMins, 2, quantile, probs = 0.975)

#--- Combine data into data frame and plot ---#
salaryMinsSummary = data.frame(zx3, salaryMinsMean, salaryMinsLow, salaryMinsHigh)

ggplot(salaryMinsSummary) +
  theme_bw() +
  geom_point(aes(x = zx3, y = salaryMinsMean), alpha = 0.6) +
  geom_ribbon(aes(x = zx3, ymin = salaryMinsLow, ymax = salaryMinsHigh), fill = "dodgerblue4", alpha = 0.3)


#--------------------------------------#
#  Within Means of Other Variables     #
#--------------------------------------#

#--- Calculate predicted y values based on means ---#
salaryMeans = matrix(0, nrow = length(zb0), ncol = length(zy))

for (i in 1:length(zb0)) {
  for (j in 1:length(zy)) {
    salaryMeans[i, j] = zb0[i] + (zb2[i] * mean(zx2)) + (zb3[i] * zx3[j]) + (zb4[i] * mean(zx4)) 
  }
}

#--- Calculate mean and HDI values ---#
salaryMeansMean = apply(salaryMeans, 2, mean)
salaryMeansLow = apply(salaryMeans, 2, quantile, probs = 0.025)
salaryMeansHigh = apply(salaryMeans, 2, quantile, probs = 0.975)

#--- Combine data into data frame and plot ---#
salaryMeansSummary = data.frame(zx3, salaryMeansMean, salaryMeansLow, salaryMeansHigh)

ggplot(salaryMeansSummary) +
  theme_bw() +
  geom_point(aes(x = zx3, y = salaryMeansMean), alpha = 0.6) +
  geom_ribbon(aes(x = zx3, ymin = salaryMeansLow, ymax = salaryMeansHigh), fill = "dodgerblue4", alpha = 0.3)


#--------------------------------------#
#    Within Max of Other Variables     #
#--------------------------------------#

#--- Calculate predicted y values based on means ---#
salaryMax = matrix(0, nrow = length(zb0), ncol = length(zy))

for (i in 1:length(zb0)) {
  for (j in 1:length(zy)) {
    salaryMax[i, j] = zb0[i] + (zb2[i] * max(zx2)) + (zb3[i] * zx3[j]) + (zb4[i] * max(zx4)) 
  }
}

#--- Calculate mean and HDI values ---#
salaryMaxMean = apply(salaryMax, 2, mean)
salaryMaxLow = apply(salaryMax, 2, quantile, probs = 0.025)
salaryMaxHigh = apply(salaryMax, 2, quantile, probs = 0.975)

#--- Combine data into data frame and plot ---#
salaryMaxSummary = data.frame(zx3, salaryMaxMean, salaryMaxLow, salaryMaxHigh)

ggplot(salaryMaxSummary) +
  theme_bw() +
  geom_point(aes(x = zx3, y = salaryMaxMean), alpha = 0.6) +
  geom_ribbon(aes(x = zx3, ymin = salaryMaxLow, ymax = salaryMaxHigh), fill = "dodgerblue4", alpha = 0.3)


###############################
#         % Take              #
#   Will base on model #2     #
###############################

#--------------------------------------#
#  Within Minimums of Other Variables  #
#--------------------------------------#

#--- Obtain Coefficients from Model ---#
mcmcData = as.data.frame(model2)
zb0 = mcmcData$b0
zb1 = mcmcData$b1
zb2 = mcmcData$b2
zb3 = mcmcData$b3
zb4 = mcmcData$b4

#--- Calculate predicted y values based on minimums ---#
percentMins = matrix(0, nrow = length(zb0), ncol = length(zy))

for (i in 1:length(zb0)) {
  for (j in 1:length(zy)) {
    percentMins[i, j] = zb0[i] + (zb2[i] * min(zx2)) + (zb3[i] * min(zx3)) + (zb4[i] * zx4[j]) 
  }
}

#--- Calculate mean and HDI values ---#
percentMinsMean = apply(percentMins, 2, mean)
percentMinsLow = apply(percentMins, 2, quantile, probs = 0.025)
percentMinsHigh = apply(percentMins, 2, quantile, probs = 0.975)

#--- Combine data into data frame and plot ---#
percentMinsSummary = data.frame(zx4, percentMinsMean, percentMinsLow, percentMinsHigh)

ggplot(percentMinsSummary) +
  theme_bw() +
  geom_point(aes(x = zx4, y = percentMinsMean), alpha = 0.6) +
  geom_ribbon(aes(x = zx4, ymin = percentMinsLow, ymax = percentMinsHigh), fill = "dodgerblue4", alpha = 0.3)


#--------------------------------------#
#  Within Means of Other Variables     #
#--------------------------------------#

#--- Calculate predicted y values based on means ---#
percentMeans = matrix(0, nrow = length(zb0), ncol = length(zy))

for (i in 1:length(zb0)) {
  for (j in 1:length(zy)) {
    percentMeans[i, j] = zb0[i] + (zb2[i] * mean(zx2)) + (zb3[i] * mean(zx3)) + (zb4[i] * zx4[j]) 
  }
}

#--- Calculate mean and HDI values ---#
percentMeansMean = apply(percentMeans, 2, mean)
percentMeansLow = apply(percentMeans, 2, quantile, probs = 0.025)
percentMeansHigh = apply(percentMeans, 2, quantile, probs = 0.975)

#--- Combine data into data frame and plot ---#
percentMeansSummary = data.frame(zx4, percentMeansMean, percentMeansLow, percentMeansHigh)

ggplot(percentMeansSummary) +
  theme_bw() +
  geom_point(aes(x = zx4, y = percentMeansMean), alpha = 0.6) +
  geom_ribbon(aes(x = zx4, ymin = percentMeansLow, ymax = percentMeansHigh), fill = "dodgerblue4", alpha = 0.3)


#--------------------------------------#
#    Within Max of Other Variables     #
#--------------------------------------#

#--- Calculate predicted y values based on means ---#
percentMax = matrix(0, nrow = length(zb0), ncol = length(zy))

for (i in 1:length(zb0)) {
  for (j in 1:length(zy)) {
    percentMax[i, j] = zb0[i] + (zb2[i] * max(zx2)) + (zb3[i] * max(zx3)) + (zb4[i] * zx4[j]) 
  }
}

#--- Calculate mean and HDI values ---#
percentMaxMean = apply(percentMax, 2, mean)
percentMaxLow = apply(percentMax, 2, quantile, probs = 0.025)
percentMaxHigh = apply(percentMax, 2, quantile, probs = 0.975)

#--- Combine data into data frame and plot ---#
percentMaxSummary = data.frame(zx4, percentMaxMean, percentMaxLow, percentMaxHigh)

ggplot(percentMaxSummary) +
  theme_bw() +
  geom_point(aes(x = zx4, y = percentMaxMean), alpha = 0.6) +
  geom_ribbon(aes(x = zx4, ymin = percentMaxLow, ymax = percentMaxHigh), fill = "dodgerblue4", alpha = 0.3)
