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
  
  # Take a random (representative) sample from the population
  sample <- pop[sample(row.names(pop), size=ss, replace=F),]
  
  # Create bootstrapped samples
  boots <- list(NULL)
  for (u in 1:B) {
    boots[[u]] <- sample[sample(row.names(sample), size=ss, replace=T),]
  }
  
  # Loop over each bootstrapped sample and create estimates
  for (u in 1:B) {
    
    sample <- boots[[u]]

    # Make matrices to put the results in
    all<-  data.frame( matrix(NA, nrow=1, ncol = n.analysis))
    colnames(all)<- c('unadj','iptw.att.misp','iptw.att.cor','iptw.att.sl','gcomp.att.misp','gcomp.att.cor','gcomp.att.sl','tmle.att.misp.misp','tmle.att.misp.cor','tmle.att.cor.misp','tmle.att.cor.cor','tmle.att.sl')
    match.full <- match.sub <- match.gen <- match.sl <- match.nn <- match.opt<- all
    
    ##############
    # 1. ANALYZE THE SAMPLE AS IS
    ############
    print("All data. No matching.")
    sample$weights <- 1
    all <- ESTIMATORS(sample)
    sample <- sample[,names(sample)!="weights"]
    
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
    
    #############
    # 5. genetic                          
    #############
    print("Genetic")
    match <- Obs <- NULL
    try(match <- suppressWarnings(matchit(A ~ W1 + W2, data = sample, method = "genetic", print.level = 4, 
                                          pop.size = 500)))
    try(Obs <- match.data(match))
    
    try(match.gen<- ESTIMATORS(Obs))
    
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
  
    ####
    # Compile
    ####
    
    estimands<- as.data.frame(t(c(all, match.nn, match.opt, match.sl, match.gen, match.sub, match.full)))
    names(estimands) <- paste0(rep(c("all.","match.nn.","match.opt.","match.sl.","match.gen.",
                                     "match.sub.","match.full."), 
                                   each = n.analysis),
                               c("unadj","iptw.att.misp","iptw.att.cor","iptw.att.sl",
                                 "gcomp.att.misp","gcomp.att.cor",
                                 "gcomp.att.sl","tmle.att.misp.misp","tmle.att.misp.cor",
                                 "tmle.att.cor.misp","tmle.att.cor.cor","tmle.att.sl"))
    
    if (u==1) results <- estimands else results <- rbind(results, estimands)
  }
  
  # Collapse to calculate the variance and mean estimate for each run
  variances <- apply(results, 2, var,  na.rm=T)
  estimands <- apply(results, 2, mean, na.rm=T)
  ci.lower <- ci.upper <- NULL
  for (t in 1:ncol(results)) {
    ci.lower[t] <- estimands[t] - 1.96*variances[t]
    ci.upper[t] <- estimands[t] + 1.96*variances[t]
  }
  
  list(estimands = estimands, ci.lower = ci.lower, ci.upper = ci.upper)
}