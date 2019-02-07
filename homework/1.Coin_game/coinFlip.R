###########################
# Answer to homework for  #
# MCMC and probability    #
# distributions.          #
###########################

#------------------------------------------#
# Create a function that will ask the user #
# how many times they want to play the     #
# game, and output the average winnings.   #
#------------------------------------------#
fairGame = function(nTrials) {
  
  #--- Initialize Vector to Hold Earnings ---#
  earnings = rep(0, times = nTrials)
  
  for (i in 1:nTrials) {
    
    #--- Flip first coin ---#
    flip1 = rbinom(n = 1, size = 1, prob = 0.5)
  
    #--- If first flip is heads (1), flip 2 more times ---#
    if (flip1 == 1) {
      
      #--- Earnings are sum of values of next two flips ---#
      earnings[i] = sum(rbinom(n = 1, size = 2, prob = 0.5))
      
    } else {
      #-- Otherwise, earnings are zero ---#
      earnings[i] = 0  
    }
  }  
  print(sum(earnings) / nTrials)
}