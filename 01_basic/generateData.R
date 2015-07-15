#--------------
# generateData:
#   input: n (sample size), covars.norm (normal or non-normal covariates), pscore.type (choice of treatment mechanism, 1, 2, or 3)
# 	output: full data
#------------------
generateData<- function(n, covars.norm=T, pscore.type = 1) {
  
  if (covars.norm == T) {
    W1 <- runif(n, min = 0.02, max = 0.7)  
    W2 <- rnorm(n, mean = (0.2 + 0.125*W1), sd = 1) 
  } else {
    W1 <- rgamma(n, 2, 10)
    W2 <- rbinom(n, 1, 0.55+0.05*W1)
  }
  
  if (pscore.type==1) {
    A <- rbinom(n, 1, plogis(-.5 + .5*W1 + .1*I(W1^2) + .5*I(W1*W2) - .5*W2))
  } else if (pscore.type==2) {
    A <- rbinom(n, 1, plogis(-.7 + 1.8*W1 -.1*I(W1^2) + 1.7*I(W1*W2) - 1.4*W2))
  } else if (pscore.type==3) {
    A <- rbinom(n, 1, plogis(-.3 + 2*W1 -2*I(W1^2) + 2*I(W1*W2) - 2.5*W2))
  }
  
  ## Generate outcomes
  Y.0 <- rnorm(n, mean=(-.5 + 2*poly(W1,2)  - W2), sd=1)
  Y.1 <- rnorm(n, mean=(-.5 + 2*poly(W1,2) - W2 + 1.5 + 2*poly(W1,2) - W2), sd=1)
 
  Y <-  ifelse(A==1,Y.1,Y.0)
  data.frame(W1,W2, A,Y, Y.1, Y.0)
}