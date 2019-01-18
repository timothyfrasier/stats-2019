############################################
# CODE FOR DATA TYPES AND PROBABILITY      #
# DISTRIBUTIONS.                           #
############################################

#-----------------------------#
#     Normal Distribution     #
#-----------------------------#

#--- Generate the Data ---#
x = seq(from = -4, to = 4, length = 200)
y = dnorm(x, mean = 0, sd = 1)

#--- Plot in Base R ---#
plot(x, y, type = "l", lwd = 2, xlab = "X Values", ylab = "Density")

#--- Plot using ggplot2 ---#
# Load library (only need to do once per R session)
library(ggplot2)

# Combine data into data frame
data = as.data.frame(cbind(x, y))

# plot
ggplot(data) +
  geom_line(aes(x = x, y = y))

# dark theme
ggplot(data) +
  theme_dark() +
  geom_line(aes(x = x, y = y))

# bw theme
ggplot(data) +
  theme_bw() +
  geom_line(aes(x = x, y = y))

# Make line thicker
ggplot(data) +
  theme_bw() + 
  geom_line(aes(x = x, y = y), size = 1.5)

# Make line opaque
ggplot(data) +
  theme_bw() + 
  geom_line(aes(x = x, y = y), size = 1.5, alpha = 0.6)

# Change colour of line
ggplot(data) +
  theme_bw() + 
  geom_line(aes(x = x, y = y), colour = "dodgerblue4", size = 1.5, alpha = 0.6)

# Make area plot
ggplot(data) +
  theme_bw() + 
  geom_area(aes(x = x, y = y), fill = "dodgerblue4", alpha = 0.6)


#-------------------------#
#  Uniform distribution   #
#-------------------------#

#--- Generate the data ---#
x = seq(from = -4, to = 4, length = 200)
y = dunif(x, min = -4, max = 4)

#--- Organize the data ---#
data = as.data.frame(cbind(x, y))

#--- Plot the data ---#
ggplot(data) +
  theme_bw() + 
  geom_line(aes(x = x, y = y), colour = "dodgerblue4", size = 1.5, alpha = 0.6)


#------------------------#
#  Cauchy Distribution   #
#------------------------#

#--- Generate the data ---#
x = seq(from = -4, to = 4, length = 200)
y = dcauchy(x, location = 0, scale = 1)

#--- Organize and plot the data ---#
data = as.data.frame(cbind(x, y))

ggplot(data) +
  theme_bw() + 
  geom_area(aes(x = x, y = y), fill = "dodgerblue4", alpha = 0.6)


#------------------------#
#  Beta Distribution     #
#------------------------#

#--- Generate the data ---#
x = seq(from = 0, to = 1, length = 200)
y = dbeta(x, shape1 = 2, shape2 = 3)

#--- Organize and plot the data ---#
data = as.data.frame(cbind(x, y))

ggplot(data) +
  theme_bw() + 
  geom_area(aes(x = x, y = y), fill = "dodgerblue4", alpha = 0.6)


#-----------------------------#
#  Exponential Distribution   #
#-----------------------------#

#--- Generate the data ---#
x = seq(from = 0, to = 20, length = 200)
y = dexp(x, rate = 0.5)

#--- Organize and plot the data ---#
data = as.data.frame(cbind(x, y))

ggplot(data) +
  theme_bw() + 
  geom_area(aes(x = x, y = y), fill = "dodgerblue4", alpha = 0.6)


#-----------------------------#
#  Binomial Distribution     #
#-----------------------------#

#--- Generate the data ---#
x = 1:50
y = dbinom(x, size = 50, prob = 0.5)

#--- Organize and plot the data ---#
data = as.data.frame(cbind(x, y))

ggplot(data) +
  theme_bw() + 
  geom_col(aes(x = x, y = y), fill = "dodgerblue4", alpha = 0.6)

#--- Make the width of the bars a little less ---#
ggplot(data) +
  theme_bw() + 
  geom_col(aes(x = x, y = y), fill = "dodgerblue4", alpha = 0.6, width = 0.3)


#-----------------------------#
#   Poisson Distribution      #
#-----------------------------#

#--- Generate the data ---#
x = seq(from = 0, to = 20, by = 1)
y = dpois(x, lambda = 4)

#--- Organize and plot the data ---#
data = as.data.frame(cbind(x, y))

ggplot(data) +
  theme_bw() + 
  geom_col(aes(x = x, y = y), fill = "dodgerblue4", alpha = 0.6, width = 0.3)
