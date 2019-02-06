## HOMEWORK FOR SIMPLE LINEAR REGRESSION

There are two components to the homework this week.

First, in class we did not compare the prior vs posterior distributions for our estimated parameters (b0, b1, and sigma). Therefore, the first part of your homework is to compare these by plotting the prior and posterior together for each parameter, and discuss if you think the choice of prior was appropriate based on this.

Second, we assumed that the noise around our predicted y-values was normally distributed. This is most often the case, but with weight there may be some serious outliers. Using a cauchy distribution can be better in such cases. Therefore, change the model to be "robust" (insensitive to outliers) by changing the likelihood to be based on the cauchy distribution rather than a normal distribution (note that this may have implications for the priors as well!). Discuss the characteristics of your estimates obtained from both models, and make a statement regarding if using the cauchy distribution was warranted.
