###############################
# Code for teaching MCMC      #
# Methods                     #
#                             #
# Requires:                   #
# The number of steps to take #
###############################

#-----------------------------#
# Basic Chain (No Peaks)      #
#-----------------------------#
basicMCMC <- function(nSteps) {

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


    #--- Step through the chain ---#
    for (i in 2:nSteps) {
        
        #--- Pick next x-value ---#
        addsub <- c(1, 2)
        choice <- sample(addsub, 1, replace = TRUE)
        
        if (choice == 1) {
            if ((xOut[i-1] + 1) <= max(xVals)) {
                xOut[i] <- (xOut[i-1] + 1)
            }
            else {
                xOut[i] <- xOut[i-1]
            }
        } else
        if (choice == 2) {
            if ((xOut[i-1] - 1) >= min(xVals)) {
                xOut[i] <- (xOut[i-1] - 1)
            }
            else {
                xOut[i] <- xOut[i-1]
            }

        }
        
        #--- Pick next y-value ---#
        addsub <- c(1, 2)
        choice <- sample(addsub, 1, replace = TRUE)
        
        if (choice == 1) {
            if ((yOut[i-1] + 0.5) <= max(yVals)) {
                yOut[i] <- (yOut[i-1] + 0.5)
            }
            else {
                yOut[i] <- yOut[i-1]
            }
        } else
        if (choice == 2) {
            if ((yOut[i-1] - 0.5) >= min(yVals)) {
                yOut[i] <- (yOut[i-1] - 0.5)
            }
            else {
                yOut[i] <- yOut[i-1]
            }
            
        }
        
    }


	#--- Plot results ---#
	plot(xOut, yOut, ylim = c(min(yVals), max(yVals)), xlim = c(min(xVals), max(xVals)), ylab = "Y Values", xlab = "X Values", main = paste("N =", nSteps), type = "o", pch = 16)
}	