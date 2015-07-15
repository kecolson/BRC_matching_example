#--------------
# Simulations to examine matching, balance measures, and performance of different matching and analysis options
# Ellie Colson
#----------------------
rm(list=ls())
options(scipen=5)

# Settings:
# Super Learner library: select algorithms
SL.library <- c("SL.glm", "SL.glm.interaction", "SL.mean", "SL.step", "SL.step.interaction")
# Normal or non-normal covariates?
covars.norm <- T
# Propensity score overlap (1=good, 2= medium, 3=poor)
pscore.type <- 2
# Use cluster?
cluster <- T

# ** Check n.subclass in outer.loop function

n <- 100000      # Size of population
R <- 1000        # Number of runs (simulation iterations)
index <- 1:R
ss <- 1000       # Sample size of each dataset (sampled from population)
n.analysis <- 12 # Number of analysis methods
n.match <- 7     # Number of matching methods

# Balance metric: choose from: 
# smd_cov (standardized mean difference in covariates)
# smd_prog (standardized mean difference in prognostic score)
# smd_pscore (standardized mean difference in propensity score)
# tstat (t test statistic)
# ks (KS test statistic)

bal_met <- "smd_cov"

# Cluster setup 
if (cluster==T) {
  library("Rmpi")
  
  # Spawn as many slaves as possible
  mpi.spawn.Rslaves()
  
  # In case R exits unexpectedly, have it automatically clean up
  # resources taken up by Rmpi (slaves, memory, etc...)
  .Last <- function(){
    if (is.loaded("mpi_initialize")){
      if (mpi.comm.size(1) > 0){
        print("Please use mpi.close.Rslaves() to close slaves.")
        mpi.close.Rslaves()
      }
      print("Please use mpi.quit() to quit R")
      .Call("mpi_finalize")
    }
  }
  
  # Tell all slaves to return a message identifying themselves
  mpi.remote.exec(paste("I am",mpi.comm.rank(),"of",mpi.comm.size()))
}

# Source edited version of tmle att so we can use weights properly
if (cluster == F) setwd("C:/Users/kecolson/Google Drive/computing/BRC_matching_example/01_basic")

# Source functions to generate data, apply estimators, apply matching methods, and calculate balance metrics
source("TMLE_ATT_edited.R") 
source("generateData.R")
source("estimators.R")
source("calcBal.R")
source("outerLoop.R") # Applies matching methods and calls ESTIMATORS and calc_bal for each

# Estimate and save the true value of the treatment effect (ATT)
if (covars.norm==F) { set.seed(333) } else { set.seed(123) }
O.true<- generateData(n=10000000, covars.norm = covars.norm, pscore.type = pscore.type)
ATT <- mean(O.true$Y.1[O.true$A==1] - O.true$Y.0[O.true$A==1])

write.csv(data.frame(ATT), "ATT.csv")

########################
# compare different designs and different estimators
########################

# Generate the population
if (covars.norm==F) { set.seed(333) } else { set.seed(123) }
pop <- generateData(n = n, covars.norm = covars.norm, pscore.type = pscore.type)

# Run all the analyses and matching combinations
if (cluster==T) { 
  # Send all necessary info to the slaves
  mpi.bcast.Robj2slave(ESTIMATORS)  
  mpi.bcast.Robj2slave(calc_bal)
  mpi.bcast.Robj2slave(outer.loop)
  mpi.bcast.Robj2slave(index)
  mpi.bcast.Robj2slave(pop)
  mpi.bcast.Robj2slave(ss) 
  mpi.bcast.Robj2slave(ATT)  
  mpi.bcast.Robj2slave(n.analysis)  
  mpi.bcast.Robj2slave(n.match)  
  mpi.bcast.Robj2slave(SL.library) 
  mpi.bcast.Robj2slave(bal_met)
  mpi.bcast.Robj2slave(.setColnames)
  mpi.bcast.Robj2slave(.bound)
  mpi.bcast.Robj2slave(regress)
  mpi.bcast.Robj2slave(predict.regress)
  mpi.bcast.Robj2slave(tmle.att2)
  mpi.bcast.Robj2slave(print.cte)
  mpi.bcast.Robj2slave(tmle.cte)
  mpi.bcast.Robj2slave(tmle.nde)
  
  # Run our big function in parallel
  results <- mpi.parLapply(index, outer.loop, pop=pop, ss=ss) 
  
} else {
  results <- lapply(index, outer.loop, pop=pop, ss=ss) 
}

# Collapse results
est     <- do.call(rbind, lapply(index, function(x) results[[x]]$estimands ))
balance <- do.call(rbind, lapply(index, function(x) results[[x]]$bal ))

# Save
write.csv(data.frame(lapply(est, as.character),     stringsAsFactors=FALSE), "estimands1.csv")
write.csv(data.frame(lapply(balance, as.character), stringsAsFactors=FALSE), paste0("balance1.",bal_met,".csv"))

# Turn off slaves
if (cluster==T) { 
  mpi.close.Rslaves() 
  mpi.quit() 
}




