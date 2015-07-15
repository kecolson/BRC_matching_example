# TMLE ATT functions
# Unedited from Sam Lendle's github

# Functions here are: .setColnames, .bound, regress, predict.regress, tmle.att2, print.cte, tmle.cte, tmle.nde

# Helpers ----------------------------

.setColnames <- function(x.colnames, x.ncols, firstChar){
  if(is.null(x.colnames)) {
		if(x.ncols > 1){
			x.colnames <- paste(firstChar,1:x.ncols, sep="")
		} else {
			x.colnames <- firstChar
		}
	} else {
		invalid.name <- nchar(x.colnames) == 0
		if(any(invalid.name)){
			x.colnames[invalid.name] <- paste(".internal",firstChar, which(invalid.name), sep="")
		}
	}
	return(x.colnames)
}


.bound <- function(x, bounds){
  x[x>max(bounds)] <- max(bounds)
	x[x<min(bounds)] <- min(bounds)
	return(x)
}


# Regress -------------------------------

regress <- function(resp, X, family=binomial(), method="glm", formula=resp ~ ., ...) {
  SL.installed <- "SuperLearner" %in% installed.packages()[,1]
  SL.version <- NULL
  if (method=="glm" || !SL.installed) {
    if (method=="SL") {
      warning("SuperLearner is not installed, using main terms glm", call.=FALSE)
      method <- "glm"
    }
    if (is.null(formula)) {
      formula <- resp ~ .
    } else {
      formula[[2]] <- as.name("resp") # changes the lhs of the formula to resp
    }
    fit <- glm(formula, data=data.frame(resp, X), family=family)
  } else if (method=="SL") {
    require(SuperLearner)
    SL.version <- packageVersion("SuperLearner")$major
    if (SL.version==1) {
      warning("Your version of SuperLearner is out of date. You should consider updating to version 2 if you don't have a good reason not to...")
      fit <- SuperLearner(resp, data.frame(X), family=family, ...)
    }
    else {
      fit <- SuperLearner(resp, data.frame(X), family=family, ...)
    }
  }
  res <- list(fit=fit, method=method, SL.version=SL.version)
  class(res) <- "regress"
  return(res)
}

predict.regress <- function(object, newdata, X=NULL, Y=NULL, ...) {
  if (object$method=="glm") {
    if (missing(newdata)) return(predict(object$fit, type="response"))
    return(predict(object$fit, newdata=newdata, type="response"))
  }
  else if (object$method=="SL") {
    if (any(is.null(Y), is.null(X)) & !missing(newdata)) warning("Original data needs to be passed to predict.regress when using SuperLearner and newdata.  predict may fail depending on the SL.library otherwise...")
    if (object$SL.version==1) {
      if (missing(newdata)) {
        return(predict(object$fit))
      } else {
        return(predict(object$fit, newdata=data.frame(newdata), X=data.frame(X), Y=Y)$fit)
      }
    } else {
      if(missing(newdata)) return(predict(object$fit)$pred)
      return(predict(object$fit, newdata=data.frame(newdata), X=data.frame(X), Y=Y)$pred)
    }
  }
}


# TMLE ATT -------------------------------

# CHANGED NAME OF FUNCTION TO TMLE.ATT2 SO WE CAN COMPARE TO PACKAGE FUNCTION
tmle.att2 <- function(A, W, Y, ...) {
  t <- tmle.cte(A, W, Y, a=1, ...)
  t$estimand <- "ATT"
  return(t)
}


# TMLE CTE ------------------------------

print.cte <- function(x, ...) {
  cat(x$estimand, ": ", x$psi,
      "\nEstimated Variance: ", x$var.psi,
      "\n95% Confidence interval: (", x$CI[1], ", ", x$CI[2],")",
      "\nP-value for H0: psi=0 vs. H1: psi!=0:", format.pval(x$pvalue, eps=0.001), "\n", sep="")
}


tmle.cte <- function(A, B, Y, a=0, Delta=NULL, Q.method="glm", Q.formula=NULL, Q.SL.library=NULL, Q.A1=NULL, Q.A0=NULL, g.method="glm", g.formula=NULL, g.SL.library=c("SL.glm", "SL.step", "SL.knn"), g.A1=NULL, gDelta.method="glm", gDelta.formula=NULL, gDelta.SL.library=c("SL.glm", "SL.step", "SL.knn"), gDelta.1=NULL, family=gaussian(), tol=1e-10, maxiter=100, target=TRUE, verbose=FALSE, Qbound=c(1e-10, 1-1e-10), gbound=c(1e-10, 1-1e-10), ...) {

  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  if (!(family$family %in% c("gaussian", "binomial"))) {
    stop("Currently only gaussian and binomial families are supported")
  }
  if (is.null(Q.SL.library)) {
    if (family$family=="gaussian") {
      Q.SL.library <- c("SL.glm", "SL.step")
    } else {
      Q.SL.library <- c("SL.glm", "SL.step", "SL.knn")
    }
  }

  if ((is.null(Delta) || all(Delta==1)) && all(!is.na(Y))) {
    missing.outcome <- FALSE
    Delta <- rep(1, length(Y))
  } else {
    missing.outcome <- TRUE
    if (is.null(Delta)) {
      warning("NA values found in Y, but Delta is NULL.  Setting Delta=0 if Y is NA, and 0 otherwise")
      Delta <- as.numeric(!is.na(Y))
    }
    if (length(Delta) != length(Y)) {
      warning("Delta is the wrong length. Setting Delta=0 if Y is NA, and 0 otherwise")
      Delta <- as.numeric(!is.na(Y))
    }
    if (any(is.na(Y)&(Delta==1))) {
      warning("Some values of Y are NA where Delta is 1.  Setting those values of Delta to 0")
      Delta[is.na(Y)] <- 0
    }
  }

  if (Q.method == "user" && (length(Q.A1) != length(Y) || length(Q.A0) != length(Y))) {
    stop ("The length of user specified Q.A1 and Q.A0 should be the same as the length of Y")
  }
  if (g.method == "user" && length(g.A1) != length(Y)) {
    stop ("The length of user specified g.A1 should be the same as the length of Y")
  }
  if (missing.outcome && gDelta.method == "user" && length(gDelta.1) != length(Y)) {
    stop ("The length of user specified gDelta.1 should be the same as the length of Y")
  }
  
  
  Y[Delta==0] <- NA
    
  Aa <- as.numeric(A==a)

  if (Q.method != "user") {
    Q.init.fit <- regress(Y[Delta==1],
                          data.frame(A=A, B)[Delta==1,],
                          family=family,
                          method=Q.method,
                          formula=Q.formula,
                          SL.library=Q.SL.library,
                          ...)

    Q.A1 <- predict(Q.init.fit, newdata=data.frame(A=1, B), X=data.frame(A=A, B)[Delta==1,], Y=Y[Delta==1])
    Q.A0 <- predict(Q.init.fit, newdata=data.frame(A=0, B), X=data.frame(A=A, B)[Delta==1,], Y=Y[Delta==1])
  }

  if (g.method != "user") {
    g.init.fit <- regress(A, B, family=binomial,
                          method=g.method,
                          formula=g.formula,
                          SL.library=g.SL.library,
                          ...)
  g.A1 <- predict(g.init.fit)
  }
  g.A1 <- .bound(g.A1, gbound)
  g.A0 <- 1-g.A1
  g.Aa <- a*g.A1 + (1-a)*g.A0

  if (missing.outcome) {
    if (gDelta.method != "user") {
      gDelta.fit <- regress(Delta,
                            data.frame(A=A, B),
                            method=gDelta.method,
                            formula=gDelta.formula,
                            SL.library=g.SL.library,
                            ...)
      gDelta.1 <- predict(gDelta.fit)
    }
    gDelta.1 <- .bound(gDelta.1, c(1, min(gbound)))
  } 
# THIS IS THE ONE PART OF THE CODE I EDITED  
#   else {
#    gDelta.1 <- rep(1, length(Y))
#   }

  fail=FALSE
  if (verbose) {
    cat("Variables in B:", names(B), "\n")
    if (!target) {
      cat("No targetting\n")
    } else {
      cat("Tolerance:", tol, "\nMaximum iterations:", maxiter, "\n")
    }
  }

  if (target) {
    iter=0
    done=FALSE
    prev.crit <- -1
    
    Q.A1 <- .bound(Q.A1, Qbound)
    Q.A0 <- .bound(Q.A0, Qbound)

    while(!done) {
      iter=iter+1

      if (iter > 1) {
        #only fluctuate g after Q has been fluctuated.
        
        H2 <- (Q.A1-Q.A0 - psi)
        
        g.up <- glm(Aa~-1+H2, offset=qlogis(g.Aa), family=binomial)
        g.eps <- coef(g.up)
        if (is.na(g.eps)) g.eps <- 0
        
        #Fluctuate g
        g.Aa <- plogis(qlogis(g.Aa)+g.eps*H2)
        
        #bound g away from 0 and 1 so logit(g) can be calculated, and clever covars
        g.Aa <- .bound(g.Aa, gbound)
        if (a==1) {
          g.A1 <- g.Aa
          g.A0 <- 1-g.A1
        } else {
          g.A0 <- g.Aa
          g.A1 <- 1-g.A0
        }
      } else {
        g.eps <- 99
      }

      #Calculate expected outcome under treatment recieved
      Q.AA <- A*Q.A1 + (1-A)*Q.A0
      #Calculate covariates for fluctuation
      H1.A1 <- (1/gDelta.1)*(g.Aa/g.A1)
      H1.A0 <- (1/gDelta.1)*(-g.Aa/g.A0)
      H1 <- A*H1.A1 + (1-A)*H1.A0
      
      Q.up <- suppressWarnings(glm(Y~-1+H1, offset=qlogis(Q.AA), family=binomial))
      Q.eps <- coef(Q.up)
      Q.A1 <- plogis(qlogis(Q.A1)+Q.eps[1]*H1.A1)
      Q.A0 <- plogis(qlogis(Q.A0)+Q.eps[1]*H1.A0)
      
      #bound Q away from 0 and 1 so logit(Q) can be calculated
      Q.A1 <- .bound(Q.A1, Qbound)
      Q.A0 <- .bound(Q.A0, Qbound)
      
      #Calculate target parameter at current step
      psi <- mean(Q.A1[A==a] - Q.A0[A==a])          
      
      crit <- max(abs(c(Q.eps, g.eps)))
      if (crit <= tol || iter>=maxiter) {
        done=TRUE
      }

      #If not done and critical value has not improved (probably because of bounding on Q)
      #then convergence failed
      #Note that if Q is not bounded, logit(Q)=+/-Inf so convergence fails either way...
      if(done==FALSE && crit==prev.crit) {
        done=TRUE
        fail=TRUE
        warning("Convergence criterion failed to improve between iterations")
      }
      prev.crit <- crit
      
      if (verbose) cat("iter:", iter, "eps:", Q.eps, g.eps, "crit:", crit, "psi:", psi, "\n")
    }
  }

  if (fail) {
    res <- list(psi=NA,
                var.psi=NA,
                CI=NA,
                pvalue=NA,
                Qstar=NA,
                gstar=NA
                )
  }

  if (!fail) {
    #Calculate final parameter of interest after convergence
    psi <-  mean(Q.A1[A==a] - Q.A0[A==a])
    Q.AA <- A*Q.A1 + (1-A)*Q.A0
    p.Aa <- mean(Aa)
    #To simplify calculation of IC
    Y[Delta==0] <- 0
    IC <- ((Delta/gDelta.1)*(A*g.Aa/g.A1 - (1-A)*g.Aa/g.A0)*(Y-Q.AA)+Aa*(Q.A1-Q.A0 - psi))/p.Aa
    if (verbose) {
      cat("Mean of the influence curve:", mean(IC), "\n")
    }
    var.psi <- var(IC)/length(Y)
    CI <- psi + c(-1.96, 1.96)*sqrt(var.psi)
    pvalue <- 2*(1-pnorm(abs(psi/sqrt(var.psi))))
    estimand <- paste("E(E(Y|A=1,B)-E(Y|A=0,B)|A=", a, ")", sep="")
    res <- list(psi=psi, var.psi=var.psi, CI=CI, pvalue=pvalue, IC=IC, estimand=estimand,
                Qstar=cbind("Q.A0"=Q.A0, "Q.A1"=Q.A1), gstar=g.A1)
  }
  class(res) <- c("cte")
  res
}


# TMLE NDE -----------------------------------------

tmle.nde <- function(A, WZ, Y, RCT=FALSE, ...) {
  t <- tmle.cte(A, WZ, Y, a=0, ...)
  t$estimand <- ifelse(RCT,
                       "Natural direct effect",
                       "Natural direct effect among the untreated")
  return(t)
}









