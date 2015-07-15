# ------------------
# Function to calculate balance metrics
# Input: bal_met (name of desired balance metric), O (the data), varname (the name of the variable for which to calculate the balance metric; 
#  only applies to certain balance metrics), weights (the matching weights, if applicable), dw (the balance metric denominator, if applicable)
# Output: the balance metric
# ------------------

calc_bal <- function(bal_met, O, varname = NULL, weights = NULL, dw = NULL) {
  if (bal_met == "smd_cov") {
    if (is.null(weights[1])) {
      return((mean(O[O$A == 1, varname]) - mean(O[O$A == 0, varname])) / dw)
    } else {
      return((weighted.mean(O[O$A==1, varname], w = O$weights[O$A==1]) - weighted.mean(O[O$A==0, varname], w=O$weights[O$A==0])) / dw)
    }
  
  } else if (bal_met == "smd_prog") {
    if (is.null(weights[1])) {
      mod <- glm(Y ~ W1*W2 + I(W1^2) + I(W2^2), data = O[O$A==0,], weights = NULL)
      prog_score <- predict(mod, newdata = O, type = "response")
      return((mean(prog_score[O$A==1]) - mean(prog_score[O$A==0])) / dw)
    } else {
      mod <- glm(Y ~ W1*W2 + I(W1^2) + I(W2^2), data = O[O$A==0,], weights = weights)
      prog_score <- predict(mod, newdata = O, type = "response")
      return((weighted.mean(prog_score[O$A==1], w = O$weights[O$A==1]) - weighted.mean(prog_score[O$A==0], w = O$weights[O$A==0])) / dw)
    }
  
  } else if (bal_met == "smd_pscore") {
    if (is.null(weights[1])) {
      mod <- glm(A ~ W1*W2 + I(W1^2) + I(W2^2), data = O, weights = NULL, family = "binomial")
      pscore <- predict(mod, type = "response")
      return((mean(pscore[O$A==1]) - mean(pscore[O$A==0])) / dw)
    } else {
      mod <- glm(A ~ W1*W2 + I(W1^2) + I(W2^2), data = O, weights = weights, family = "binomial")
      pscore <- predict(mod, type = "response")
      return((weighted.mean(pscore[O$A==1], w = O$weights[O$A==1]) - weighted.mean(pscore[O$A==0], w = O$weights[O$A==0])) / dw)
    }
  
  } else if (bal_met == "tstat") {
    if (is.null(weights[1])) {
      return(t.test(x = O[O$A==1,varname], y = O[O$A==0,varname], alternative = "two.sided")$statistic)
    } else {
      return(wtd.t.test(x = O[O$A==1,varname], y = O[O$A==0,varname], 
                       weight = O$weights[O$A==1], 
                       weighty = O$weights[O$A==0])$coefficients[1])
    }  
  } else if (bal_met == "ks") {
    return(ks.test(x = O[O$A==1,varname], y = O[O$A==0,varname], alternative = "two.sided")$statistic) 
  }
}