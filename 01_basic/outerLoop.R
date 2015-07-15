#---------
# outer.loop: function to run each matching method, set of estimators, and balance metrics
#   input: iteration (the iteration number to set the seed), pop (the population data from which to sample), ss (the sample size to draw)
#   output: parameter estimates and balance metrics: estimands, bal
#--------

outer.loop <- function(iteration, pop, ss) {

  require("MatchIt")
  require("Hmisc")
  require("weights")
  require("SuperLearner")
  
  print(iteration)
  set.seed(iteration)

  # Make matrices to put the results in
  all<-  data.frame( matrix(NA, nrow=1, ncol = n.analysis))
  colnames(all)<- c('unadj','iptw.att.misp','iptw.att.cor','iptw.att.sl','gcomp.att.misp','gcomp.att.cor','gcomp.att.sl','tmle.att.misp.misp','tmle.att.misp.cor','tmle.att.cor.misp','tmle.att.cor.cor','tmle.att.sl')
  match.full <- match.sub <- match.gen <- match.sl <- match.nn <- match.opt<- all
  
  # Take a random (representative) sample from the population
  sample <- pop[sample(row.names(pop), size=ss, replace=F),]
  
  # Create data frames to store balance metrics
  # Calculate denominators (dw) for balance metrics as needed
  if (bal_met %in% c("smd_cov","tstat","ks")) {
    print(bal_met)
    bal <- data.frame(matrix(NA,nrow=1, ncol = 2*n.match))
    colnames(bal)<-paste0( rep(c("all.","match.nn.","match.opt.","match.sl.","match.gen.","match.sub.","match.full."),each=2), rep(c("W1","W2"), n.match) )

    if (bal_met == "smd_cov") {
      dW1 <- sqrt((var(sample$W1[sample$A==1]) + var(sample$W1[sample$A==0])) / 2)
      dW2 <- sqrt((var(sample$W2[sample$A==1]) + var(sample$W2[sample$A==0])) / 2)
    }
  
  } else if (bal_met %in% c("smd_prog","smd_pscore")) {
    print(bal_met)
    bal <- data.frame(matrix(NA,nrow=1, ncol = n.match))
    colnames(bal)<- c("all","match.nn","match.opt","match.sl","match.gen","match.sub","match.full")

    if (bal_met == "smd_prog") {
      mod <- glm(Y ~ W1*W2 + I(W1^2) + I(W2^2), data = sample[sample$A==0,])
      score <- predict(mod, newdata = sample, type = "response")
    
    } else if (bal_met == "smd_pscore") {
      mod <- glm(A ~ W1*W2 + I(W1^2) + I(W2^2), data = sample, family = "binomial")
      score <- predict(mod, type = "response")
    }
    
    dW <- sqrt(var(score))
  }
  
  ##############
  # 1. ANALYZE THE SAMPLE AS IS
  ############
  print("All data. No matching.")
  sample$weights <- 1
  all <- ESTIMATORS(sample)
  sample <- sample[,names(sample)!="weights"]
  
  if (bal_met == "smd_cov") {
    bal$all.W1 <- calc_bal(bal_met, sample, "W1", weights = NULL, dw = dW1)
    bal$all.W2 <- calc_bal(bal_met, sample, "W2", weights = NULL, dw = dW2)
  } else if (bal_met %in% c("smd_prog","smd_pscore")) {
    bal$all <- calc_bal(bal_met, sample, varname = NULL, weights = NULL, dw = dW)
  } else if (bal_met %in% c("tstat","ks")) {
    bal$all.W1 <- calc_bal(bal_met, sample, "W1", weights = NULL, dw = NULL)
    bal$all.W2 <- calc_bal(bal_met, sample, "W2", weights = NULL, dw = NULL)
  }
  
  ############# 
  # Different matching designs
  ##############
  
  #############
  # 2. Nearest neighbor matching
  #############
  print("Nearest neighbor")
  match <- Obs <- NULL
  try(match <- suppressWarnings(matchit(A ~ W1 + W2, data=sample, method="nearest", replace = T)))
  try(Obs <- match.data(match))
  
  try(match.nn<- ESTIMATORS(Obs))
  
  if (bal_met == "smd_cov") {
    bal$match.nn.W1 <- calc_bal(bal_met, Obs, "W1", weights = NULL, dw = dW1)
    bal$match.nn.W2 <- calc_bal(bal_met, Obs, "W2", weights = NULL, dw = dW2)
  } else if (bal_met %in% c("smd_prog","smd_pscore")) {
    bal$match.nn <- calc_bal(bal_met, Obs, varname = NULL, weights = NULL, dw = dW)
  } else if (bal_met %in% c("tstat","ks")) {
    bal$match.nn.W1 <- calc_bal(bal_met, Obs, "W1", weights = NULL, dw = NULL)
    bal$match.nn.W2 <- calc_bal(bal_met, Obs, "W2", weights = NULL, dw = NULL)
  }
  
  #############
  # 3. Optimal matching 
  #############
  print("Optimal")
  
  # Only works if there are more control units than treated units 
  if (sum(sample$A)< (nrow(sample)-sum(sample$A))) {
    match <- Obs <- NULL
    try(match <- suppressWarnings(matchit(A ~ W1 + W2, data=sample, method="optimal")))
    try(Obs <- match.data(match))
    
    try(match.opt<- ESTIMATORS(Obs))
    
    if (bal_met == "smd_cov") {
      bal$match.opt.W1 <- calc_bal(bal_met, Obs, "W1", weights = NULL, dw = dW1)
      bal$match.opt.W2 <- calc_bal(bal_met, Obs, "W2", weights = NULL, dw = dW2)
    } else if (bal_met %in% c("smd_prog","smd_pscore")) {
      bal$match.opt <- calc_bal(bal_met, Obs, varname = NULL, weights = NULL, dw = dW)
    } else if (bal_met %in% c("tstat","ks")) {
      bal$match.opt.W1 <- calc_bal(bal_met, Obs, "W1", weights = NULL, dw = NULL)
      bal$match.opt.W2 <- calc_bal(bal_met, Obs, "W2", weights = NULL, dw = NULL)
    }
  }
  
  #############
  # 4. Estimate propenisty score with superlearner, then nearest neighbor match
  #############
  print("Nearest neighbor with superlearner")
  SL.out <- match <- m <- Obs <- NULL
  SL.out <- SuperLearner(Y = sample$A, X = sample[,c("W1","W2")], SL.library = SL.library, family = 'binomial', cvControl = list(V=10))
  try(match <- matchit(A ~ W1 + W2, data = sample, method = "nearest", distance = SL.out$SL.predict, replace = T))
  try(Obs <- match.data(match))
  
  try(match.sl<- ESTIMATORS(Obs))
  
  if (bal_met == "smd_cov") {
    bal$match.sl.W1 <- calc_bal(bal_met, Obs, "W1", weights = NULL, dw = dW1)
    bal$match.sl.W2 <- calc_bal(bal_met, Obs, "W2", weights = NULL, dw = dW2)
  } else if (bal_met %in% c("smd_prog","smd_pscore")) {
    bal$match.sl <- calc_bal(bal_met, Obs, varname = NULL, weights = NULL, dw = dW)
  } else if (bal_met %in% c("tstat","ks")) {
    bal$match.sl.W1 <- calc_bal(bal_met, Obs, "W1", weights = NULL, dw = NULL)
    bal$match.sl.W2 <- calc_bal(bal_met, Obs, "W2", weights = NULL, dw = NULL)
  }
  
  #############
  # 5. genetic                          
  #############
  print("Genetic")
  match <- Obs <- NULL
  try(match <- suppressWarnings(matchit(A ~ W1 + W2, data = sample, method = "genetic", print.level = 4, 
                                        pop.size = 500)))
  try(Obs <- match.data(match))
  
  try(match.gen<- ESTIMATORS(Obs))
  
  # Use weights in balance metric too
  if (bal_met == "smd_cov") {
    bal$match.gen.W1 <- calc_bal(bal_met, Obs, "W1", weights = T, dw = dW1)
    bal$match.gen.W2 <- calc_bal(bal_met, Obs, "W2", weights = T, dw = dW2)
  } else if (bal_met %in% c("smd_prog","smd_pscore")) {
    bal$match.gen <- calc_bal(bal_met, Obs, varname = NULL, weights = T, dw = dW)
  } else if (bal_met %in% c("tstat","ks")) {
    bal$match.gen.W1 <- calc_bal(bal_met, Obs, "W1", weights = T, dw = NULL)
    bal$match.gen.W2 <- calc_bal(bal_met, Obs, "W2", weights = T, dw = NULL)
  }
  
  #############
  # 6. subclassification
  # Analyses within each subclass, then average
  #############
  print("Subclassification")
  match <- Obs <- NULL
  n.subclass <- 10 # Number of subclasses
  try(match <- suppressWarnings(matchit(A ~ W1 + W2, data=sample, method="subclass", subclass = n.subclass)))
  try(Obs <- match.data(match))
  
  # Estimate within each subclass and then average across the subclasses
  # Since we are estimating the ATT here, there is no need to weight this average across the subclasses, because each subclass has the same number of treated units.
  ests <- matrix(NA, nrow = n.subclass, ncol = ncol(match.sub))
  for (s in 1:n.subclass) {
    try(ests[s,]<- ESTIMATORS(Obs[Obs$subclass==s,]))  
  }

  try(match.sub<- apply(ests,2, function(x) mean(x)))

  # Use weights 
  if (bal_met == "smd_cov") {
    bal$match.sub.W1 <- calc_bal(bal_met, Obs, "W1", weights = T, dw = dW1)
    bal$match.sub.W2 <- calc_bal(bal_met, Obs, "W2", weights = T, dw = dW2)
  } else if (bal_met %in% c("smd_prog","smd_pscore")) {
    bal$match.sub <- calc_bal(bal_met, Obs, varname = NULL, weights = T, dw = dW)
  } else if (bal_met %in% c("tstat","ks")) {
    bal$match.sub.W1 <- calc_bal(bal_met, Obs, "W1", weights = T, dw = NULL)
    bal$match.sub.W2 <- calc_bal(bal_met, Obs, "W2", weights = T, dw = NULL)
  }
  
  #############
  # 7. full
  # Sample sizes within each subclass are very small
  # Estimate using weights
  #############
  print("Full matching")
  match <- Obs <- NULL
  try(match <- suppressWarnings(matchit(A ~ W1 + W2, data=sample, method="full")))
  try(Obs <- match.data(match))
  
  try(match.full<- ESTIMATORS(Obs))   
  
  # Use weights
  if (bal_met == "smd_cov") {
    bal$match.full.W1 <- calc_bal(bal_met, Obs, "W1", weights = T, dw = dW1)
    bal$match.full.W2 <- calc_bal(bal_met, Obs, "W2", weights = T, dw = dW2)
  } else if (bal_met %in% c("smd_prog","smd_pscore")) {
    bal$match.full <- calc_bal(bal_met, Obs, varname = NULL, weights = T, dw = dW)
  } else if (bal_met %in% c("tstat","ks")) {
    bal$match.full.W1 <- calc_bal(bal_met, Obs, "W1", weights = T, dw = NULL)
    bal$match.full.W2 <- calc_bal(bal_met, Obs, "W2", weights = T, dw = NULL)
  }

  estimands<- as.data.frame(t(c(all, match.nn, match.opt, match.sl, match.gen, match.sub, match.full)))
  names(estimands) <- paste0(rep(c("all.","match.nn.","match.opt.","match.sl.","match.gen.",
                                   "match.sub.","match.full."), 
                                 each = n.analysis),
                             c("unadj","iptw.att.misp","iptw.att.cor","iptw.att.sl",
                               "gcomp.att.misp","gcomp.att.cor",
                               "gcomp.att.sl","tmle.att.misp.misp","tmle.att.misp.cor",
                               "tmle.att.cor.misp","tmle.att.cor.cor","tmle.att.sl"))
  
  list(estimands = estimands, bal = bal)
}