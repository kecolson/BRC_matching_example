#---------
# ESTIMATORS: function to implement the different estimators with weights
#   input: O (the observed data)
#   output: point estimates from the different analysis methods
#--------

ESTIMATORS<- function(O) {
  
  require("SuperLearner")
  
  unadj<- iptw.att.misp<- iptw.att.cor<- iptw.att.sl<- gcomp.att.misp<- gcomp.att.cor<- gcomp.att.sl<- tmle.att.misp.misp <-
    tmle.att.misp.cor<- tmle.att.cor.misp<- tmle.att.cor.cor<- tmle.att.sl<- NA
  
  txt<- control <- O
  txt$A <-1;  control$A <- 0
  n <- nrow(O)
  tx.indices <- O$A==1
  control.indices <- O$A==0
  
  ## Unadjusted
  print("Unadjusted")
  unadj<- weighted.mean(O$Y[tx.indices], w = O$weights[tx.indices]) - weighted.mean(O$Y[control.indices], w = O$weights[control.indices])
  
  ### IPTW for ATT, misspecified
  print("IPTW parametric misspecified")
  g.model<- glm(A ~ W1 + W2, family='binomial', data=O)
  g.pred<- predict(g.model, type='resp')
  O$txweight <- ifelse(O$A==1, 1, g.pred/(1-g.pred))
  mod <- glm(Y ~ A, weights = txweight, data = O)
  iptw.att.misp <- mean(predict(mod, newdata = txt, type='response')) - 
    mean(predict(mod, newdata = control, type='response'))
  O <- O[,names(O) != "txweight"]
  
  ### IPTW for ATT, correctly specified
  print("IPTW parametric correctly specified")
  g.model<- glm(A ~ I(poly(W1,2)) + W2 + W1:W2, family='binomial', data = O)
  g.pred<- predict(g.model, type = 'resp')
  O$txweight <- ifelse(O$A==1, 1, g.pred/(1-g.pred))
  mod <- glm(Y ~ A, weights = txweight, data = O)
  iptw.att.cor <- mean(predict(mod, newdata = txt, type='response')) - 
    mean(predict(mod, newdata = control, type='response'))
  O <- O[,names(O) != "txweight"]
  
  # IPTW for ATT, with superlearner
  print("IPTW with superlearner")
  SL.out <- SuperLearner(Y = O$A, X = O[,c("W1","W2")], SL.library = SL.library, 
                         family = 'binomial', cvControl=list(V=10))
  g.pred <- SL.out$SL.predict
  O$txweight <- ifelse(O$A==1, 1, g.pred/(1-g.pred))
  mod <- glm(Y ~ A, weights = txweight, data = O)
  iptw.att.sl <- mean(predict(mod, newdata = txt, type='response')) - 
    mean(predict(mod, newdata = control, type='response'))
  O <- O[,names(O) != "txweight"]
  
  ### gcomputation with mispecified linear model
  print("Gcomp parametric misspecified")
  reg.model <- glm(Y~ A +W2+W1, data=O, weights=O$weights)
  gcomp.att.misp<- mean(predict(reg.model, newdata = txt[tx.indices,], type='response')) - mean(predict(reg.model, newdata = control[tx.indices,], type='response'))
  
  ### g-comp with correctly specified linear model
  print("Gcomp parametric correctly specified")
  reg.model <- glm(Y~ I(poly(W1,2)) + W2 + A + A:I(poly(W1,2)) + A:W2, data=O, weights=O$weights)
  gcomp.att.cor<- mean(predict(reg.model, newdata = txt[tx.indices,], type='response')) - mean(predict(reg.model, newdata = control[tx.indices,], type='response'))
 
  ### g-comp with superlearner
  print("Gcomp with SL")
  newX <- rbind(O, txt, control)[,c("A","W1","W2")]
  SL.out <- SuperLearner(Y = O$Y, X = O[,c("A","W1","W2")], SL.library = SL.library, newX = newX,
                           family = 'gaussian', cvControl=list(V=10), obsWeights = O$weights)
  QbarAW <- SL.out$SL.predict[1:n]
  Qbar1W <- SL.out$SL.predict[(n+1):(2*n)]
  Qbar0W <- SL.out$SL.predict[(2*n+1):(3*n)]
  gcomp.att.sl<- mean(Qbar1W[tx.indices]) - mean(Qbar0W[tx.indices])
  
  # ### TMLE
  O$Y.scaled <- (O$Y-min(O$Y))/ (max(O$Y)-min(O$Y)) # Scale Y to be between 0 and 1 so it will play nice with tmle.att
  
  # TMLE for ATT, mispecified treatment mechanism, mispecified outcome
  print("TMLE parametric, g misspecified, Q misspecified")
  tmle.out<- tmle.att2(Y=O$Y.scaled, A=O$A,  # Using tmle.att2 here, which is my adjusted function to make weights work properly
                     W = O[,c("W1","W2")], 
                     family = "gaussian", 
                     Delta = rep(1,nrow(O)), 
                     gDelta.method = "user",
                     gDelta.1 = 1/O$weights,
                     g.method = "glm",
                     g.formula = A ~ W1 + W2, 
                     Q.method = "glm",
                     Q.formula = Y ~ A + W1 + W2)
  tmle.att.misp.misp <- tmle.out$psi * (max(O$Y)-min(O$Y))
 
  # TMLE for ATT, mispecified treatment mechanism, correct outcome
  print("TMLE parametric, g misspecified, Q correctly specified")
  tmle.out<- tmle.att2(Y=O$Y.scaled, A=O$A, 
                       W = O[,c("W1","W2")], 
                       family = "gaussian", 
                       Delta = rep(1,nrow(O)), 
                       gDelta.method = "user",
                       gDelta.1 = 1/O$weights,
                       g.method = "glm",
                       g.formula = A ~ W1 + W2, 
                       Q.method = "glm",
                       Q.formula = Y ~ I(poly(W1,2)) + W2 + A + A:I(poly(W1,2)) + A:W2)
  tmle.att.misp.cor <- tmle.out$psi * (max(O$Y)-min(O$Y))
 
  # TMLE for ATT, correct treatment mechanism, mispecified outcome
  print("TMLE parametric, g correctly specified, Q misspecified")
  tmle.out<- tmle.att2(Y=O$Y.scaled, A=O$A,  
                       W = O[,c("W1","W2")], 
                       family = "gaussian", 
                       Delta = rep(1,nrow(O)), 
                       gDelta.method = "user",
                       gDelta.1 = 1/O$weights,
                       g.method = "glm",
                       g.formula = A ~ I(poly(W1,2)) + W2 + W1:W2, 
                       Q.method = "glm",
                       Q.formula = Y ~ A + W1 + W2)
  tmle.att.cor.misp <- tmle.out$psi * (max(O$Y)-min(O$Y))
 
  # TMLE for ATT, correct treatment mechanism, correct outcome
  print("TMLE parametric, g correctly specified, Q correctly specified")
  tmle.out<- tmle.att2(Y=O$Y.scaled, A=O$A,  
                       W = O[,c("W1","W2")], 
                       family = "gaussian", 
                       Delta = rep(1,nrow(O)), 
                       gDelta.method = "user",
                       gDelta.1 = 1/O$weights,
                       g.method = "glm",
                       g.formula = A ~ I(poly(W1,2)) + W2 + W1:W2, 
                       Q.method = "glm",
                       Q.formula = Y ~ I(poly(W1,2)) + W2 + A + A:I(poly(W1,2)) + A:W2)
  tmle.att.cor.cor <- tmle.out$psi * (max(O$Y)-min(O$Y))
 
  # TMLE for ATT, with superlearner
  print("TMLE with superlearner")
  tmle.out<- tmle.att2(Y=O$Y.scaled, A=O$A,  
                       W = O[,c("W1","W2")], 
                       family = "gaussian", 
                       Delta = rep(1,nrow(O)), 
                       gDelta.method = "user",
                       gDelta.1 = 1/O$weights,
                       g.method = "SL",
                       g.SL.library = SL.library, 
                       Q.method = "SL",
                       Q.SL.library = SL.library)
  tmle.att.sl <- tmle.out$psi * (max(O$Y)-min(O$Y))
 
  # Return the ATT estimates from all the different methods
  c(unadj, iptw.att.misp, iptw.att.cor, iptw.att.sl, gcomp.att.misp, gcomp.att.cor, gcomp.att.sl, tmle.att.misp.misp, tmle.att.misp.cor, tmle.att.cor.misp, tmle.att.cor.cor, tmle.att.sl)
  
}	
