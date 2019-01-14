###############################
# Code for teaching MCMC      #
# Methods                     #
#                             #
# Requires:                   #
# The number of steps to take,#
# how heated the chain should #
# be (from 0 to 1, where 1 is #
# more heated), and how long  #
# the burn-in should be.      #
#                             #
# by Tim Frasier              #
# January, 2015               #
###############################

#-----------------------------#
# Basic Chain (One Peak)      #
#-----------------------------#
onePeakMCMC <- function(nSteps, heat, burnin) {

	#--- Settings for MCMC ---#
	xmin = -10
	xmax = 10
	ymin = -10
	ymax = 10

	#--- Vectors for holding results ---#
	yOut <- rep(0, nSteps)
	xOut <- rep(0, nSteps)

	#--- Create Possible Values ---#
	yVals <- seq(from = ymin, to = ymax, by = 1)
	xVals <- seq(from = xmin, to = xmax, by = 1)
    
    #--- Create vectors of probabilities ---#
    yProb <- rep(0, nSteps)
    xProb <- rep(0, nSteps)
    
    #--- Create function for calculating normal probabilities ---#
    normalx <- function(x, mu, sig) {
        (1 / (sig * (sqrt(2 * pi)))) * exp(-(((x - mu)^2)/(2 * sig^2)))
    }
    
    #--- Set values for mean and standard deviation of peak ---#
    mu <- 3
    sig <- 2
    
    #--- Plot ellipses for peak ---#
    #--- Load needed libraries ---#
    library(MASS)
    library(cluster)
    
    #--- Create points representing density ---#
    xpeak <- rnorm(1000, mean = mu, sd = sig)
    ypeak <- rnorm(1000, mean = mu, sd = sig)
    peak <- cbind(xpeak, ypeak)
    
    #--- Find 95% highest density ---#
    fit <- cov.mve(peak, quantile.used = nrow(peak) * 0.95)
    points_in_ellipse <- peak[fit$best, ]
    ellipse_boundary <- predict(ellipsoidhull(points_in_ellipse))
    
    #--- Draw it ---#
    plot(peak, col = "white", ylim = c(min(yVals), max(yVals)), xlim = c(min(xVals), max(xVals)), ylab = "Y Values", xlab = "X Values", main = paste("N =", nSteps))
    lines(ellipse_boundary, col = "red", lwd = 2)
    
    #--- Find 75% highest density ---#
    fit <- cov.mve(peak, quantile.used = nrow(peak) * 0.75)
    points_in_ellipse <- peak[fit$best, ]
    ellipse_boundary <- predict(ellipsoidhull(points_in_ellipse))
    
    #--- Draw it ---#
    lines(ellipse_boundary, col = "red", lwd = 2)
    
    readline("Press <return> to see if chain can fine the peak")
    
    #--- Find initial values ---#
    
        #--- Pick x or y as limit ---#
        #--- x = 1, y = 2         ---#
        choice <- c(1, 2)
        maxmin <- c(1, 2)
        startaxis <- sample(choice, 1, replace = TRUE)
        startpos <- sample(maxmin, 1, replace = TRUE)
        
        #--- Choose starting position ---#
        if (startaxis == 1) {
            if (startpos == 1) {
                xOut[1] <- min(xVals)
            } else
            if (startpos == 2) {
                xOut[1] <- max(xVals)
            }
            yOut[1] <- sample(yVals, 1, replace = TRUE)
        } else
        if (startaxis == 2) {
            if (startpos == 1) {
                yOut[1] <- min(yVals)
            } else
            if (startpos == 2) {
                yOut[1] <- max(yVals)
            }
            xOut <- sample(xVals, 1, replace = TRUE)
        }

    #--- Calculate Probabilities for Initial Values ---#
    yProb[1] <- normalx(yOut[1], mu, sig)
    xProb[1] <- normalx(xOut[1], mu, sig)

    
    #--- Step through the chain ---#
    for (i in 2:nSteps) {
        
        #--- Pick next x-value ---#
        addsub <- c(1, 2)
        choice <- sample(addsub, 1, replace = TRUE)
        
        #--- Calculate potential x probabilites ---#
        xtemp1 <- normalx((xOut[i-1] + 1), mu, sig)
        xtemp2 <- normalx((xOut[i-1] - 1), mu, sig)
        
        #--- Propose new steps and decide which to take ---#
        if (choice == 1) {
            if ((xOut[i-1] + 1) <= max(xVals)) {
                if (xtemp1 >= xProb[i-1]) {
                    xOut[i] <- (xOut[i-1] + 1)
                } else
                if (runif(1, 0, 1) < (abs(xtemp1 - xProb[i-1])) + heat) {
                    xOut[i] <- (xOut[i-1] + 1)
                }
                else {
                    xOut[i] <- xOut[i-1]
                }
            }
            else {
                xOut[i] <- xOut[i-1]
            }
        } else
        
        if (choice == 2) {
            if ((xOut[i-1] - 1) >= min(xVals)) {
                if (xtemp2 >= xProb[i-1]) {
                    xOut[i] <- (xOut[i-1] - 1)
                } else
                if (runif(1, 0, 1) < (abs(xtemp2 - xProb[i-1])) + heat) {
                    xOut[i] <- (xOut[i-1] - 1)
                }
                else {
                    xOut[i] <- xOut[i-1]
                }
            }
            else {
                xOut[i] <- xOut[i-1]
            }

        }
        
        #--- Calculate new xProb value ---#
        xProb[i] <- normalx(xOut[i], mu, sig)
        
        
        #--- Pick next y-value ---#
        addsub <- c(1, 2)
        choice <- sample(addsub, 1, replace = TRUE)
        
        #--- Calculate potential y probabilites ---#
        ytemp1 <- normalx((yOut[i-1] + 1), mu, sig)
        ytemp2 <- normalx((yOut[i-1] - 1), mu, sig)
        
        #--- Propose new steps and decide which to take ---#
        if (choice == 1) {
            if ((yOut[i-1] + 1) <= max(yVals)) {
                if (ytemp1 >= yProb[i-1]) {
                    yOut[i] <- (yOut[i-1] + 1)
                } else
                if (runif(1, 0, 1) < (abs(ytemp1 - yProb[i-1])) + heat) {
                    yOut[i] <- (yOut[i-1] + 1)
                }
                else {
                    yOut[i] <- yOut[i-1]
                }
            }
            else {
                yOut[i] <- yOut[i-1]
            }
        } else
        
        if (choice == 2) {
            if ((yOut[i-1] - 1) >= min(yVals)) {
                if (ytemp2 >= yProb[i-1]) {
                    yOut[i] <- (yOut[i-1] - 1)
                } else
                if (runif(1, 0, 1) < (abs(ytemp2 - yProb[i-1])) + heat) {
                    yOut[i] <- (yOut[i-1] - 1)
                }
                else {
                    yOut[i] <- yOut[i-1]
                }
            }
            else {
                yOut[i] <- yOut[i-1]
            }
        }
        
        #--- Calculate new yProb value ---#
        yProb[i] <- normalx(yOut[i], mu, sig)
        
    }

    #--- Separate points into burn-in and recorded steps ---#
    
        #--- Create vectors to hold data ---#
        xOutBurn <- rep(0, burnin)
        yOutBurn <- rep(0, burnin)
        
        xOutRecord <- rep(0, (nSteps - burnin))
        yOutRecord <- rep(0, (nSteps - burnin))
        
        #--- Parse data ---#
        for (i in 1:nSteps) {
            
            if (burnin >= i) {
                xOutBurn[i] <- xOut[i]
                yOutBurn[i] <- yOut[i]
            } else {
                xOutRecord[i - burnin] <- xOut[i]
                yOutRecord[i - burnin] <- yOut[i]
            }
        }
        
	#--- Plot results ---#
    if (burnin == 0) {
        
        par(new = TRUE)
        plot(xOut, yOut, ylim = c(min(yVals), max(yVals)), xlim = c(min(xVals), max(xVals)), ylab = "", xlab = "", main = "", type = "o", pch = 16)
    } else {
    
        par(new = TRUE)
        plot(xOutBurn, yOutBurn, ylim = c(min(yVals), max(yVals)), xlim = c(min(xVals), max(xVals)), ylab = "Y Values", xlab = "X Values", main = paste("N =", nSteps), type = "o", pch = 16, col = "gray")
        
        par(new = TRUE)
        
        plot(xOutRecord, yOutRecord, ylim = c(min(yVals), max(yVals)), xlim = c(min(xVals), max(xVals)), axes = FALSE, ylab = "", xlab = "", main = "", type = "o", pch = 16)
    }
    
}