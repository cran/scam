## (c) Natalya Pya (2022-2023). Released under GPL2.
## based on (c) Simon N Wood (bfgs(mgcv))

#######################################################################
## Function to fit SCAM using quasi-Newton method, the BFGS method   ##     
#######################################################################

scam.fit1 <- function(G,sp, gamma=1, etastart=NULL, mustart=NULL, env=env, 
              null.coef=rep(0,ncol(G$X)), control=scam.control(),
             ## next are additional input arguments for BFGS:
             conv.tol=1e-6,maxNstep=3,maxSstep=2,maxHalf=30,printWarn=FALSE ) 
## Function to fit SCAM using a quasi-Newton method, the BFGS method (as an alternative to the Newton method)
## by minimizing penalized deviance wrt model coefficients

## BFGS is based on Nocedal & Wright (2006) Numerical Optimization, Springer, 
## and as applied in bfgs() of Simon Wood mgcv for smoothing parameter estimation.  
## Simon, mgcv: "In particular the step lengths are chosen to meet the Wolfe conditions
## using their algorithms 3.5 (p60) and 3.6 (p61). On p143 they recommend a post step
## adjustment to the initial Hessian. I can't understand why one would do anything
## other than adjust so that the initial Hessian would give the step taken, and
## indeed the latter adjustment seems to give faster convergence than their 
## proposal, and is therefore implemented."

## MAIN STEPS:
## 1. Initialization of parameters as in scam.fit()
## 2. shrinking towards null.coef if immediately after initialization invalid, remains as in scam.fit()
## 3. main BFGS iterations ...
## 4. define matrices at their converged values from the BFGS method..

## below is valid for scam.fit via Newton: 
##  maxit=200, devtol.fit=1e-8, steptol.fit=1e-8, trace=FALSE, print.warn=FALSE
   ## G - list of items from gam(...,fit=FALSE) needed to fit a scam
   ## sp- vector of smoothing parameters
   ## not.exp - if TRUE then notExp() function will be used in place of exp in positivity ensuring beta parameters re-parameterization
   ## null.coef - coefficients for a null model, in order to be able to check for immediate divergence. null.coef give some sort of upper bound on deviance. This allows immediate divergence problems to be controlled (from mgcv). here some of coeff need to be exponentiated. 
   ## control - control list; it includes:
      ##* maxit - a positive scalar which gives the maximum number of iterations for Newton's method
      ##* devtol.fit - a scalar giving the tolerance at which the relative penalized deviance is considered to be close enougth to 0 to terminate the algorithm
      ##* steptol.fit - a scalar giving the tolerance at which the scaled distance between two successive iterates is considered close enough to zero to terminate the algorithm 
      ##* trace turns on or off some de-bugging information.
      ##* print.warn =FALSE - when set to 'FALSE' turns off printing warning messages for step halving under non-finite exponentiated coefficients,  non-finite deviance and/or if mu or eta are out of bounds. 
##
{ score <- function(y,X,beta,offset,weights, rS,S.t, q,iv,linkinv, dev.resids,mu.eta,variance) {
  ## function to calculate penalized deviance and its gradient
    beta.t <- beta ## re-parameterized beta
    beta.t[iv] <- if (!not.exp) exp(beta[iv]) else notExp(beta[iv]) ## values of re-para beta of the shape constrained terms
    eta <- as.numeric(X%*%beta.t + as.numeric(offset))  ## linear predictor
    mu <- linkinv(eta)  ## fitted values
    dev <- sum(dev.resids(y,mu,weights)) ## deviance
    pdev <- dev + sum((rS%*%beta)^2) ## penalized deviance 
    Cdiag <- rep(1,q)
    Cdiag[iv] <- if (!not.exp) beta.t[iv] else DnotExp(beta[iv])
    tX1 <- Cdiag*t(X)
    g.deriv <- 1/mu.eta(eta)  # diag(G)
    w1 <- weights/(variance(mu)*g.deriv^2)    # diag(W1)
    grad <- -drop(tX1%*%(w1*g.deriv*(y-mu)))+S.t%*%beta # the gradient vector of the penalized deviance
    list(pdev=pdev, grad=grad)
  } ## end score()

  zoom <- function(lo,hi) { 
  ## as in bfgs() of mgcv:
  ## local function implementing Algorithm 3.6 of Nocedal & Wright
  ## (2006, p61) Numerical Optimization. Relies on R scoping rules. 
  ## alpha.lo and alpha.hi are the bracketing step lengths.
  ## This routine bisection searches for a step length that meets the
  ## Wolfe conditions. lo and hi are both objects containing fields
  ## `score', `alpha', `dscore', where `dscore' is the derivative of 
  ## the score in the current step direction, `grad' and `mustart'. 
  ## `dscore' will be NULL if the gradiant has yet to be evaluated.
    for (i in 1:40) {
      trial <- list(alpha = (lo$alpha+hi$alpha)/2)
      beta <- ibeta +step * trial$alpha
      b <- score(y=y,X=X,beta=beta,offset=offset,weights=weights,rS=rS,S.t=S.t,q=q,iv=iv,
          linkinv=linkinv,dev.resids=dev.resids,mu.eta=mu.eta,variance=variance)
      trial$pdev <- b$pdev
      
      if (trial$pdev >initial$pdev+trial$alpha*c1*initial$dscore||trial$pdev>=lo$pdev) {
        hi <- trial ## failed Wolfe 1 - insufficient decrease - step too long
      } else { ## met Wolfe 1 so check Wolve 2 - sufficiently positive second derivative?
          trial$grad <- b$grad
          trial$dscore <- sum(step*trial$grad) ## directional derivative
          if (abs(trial$dscore) <= -c2*initial$dscore) return(trial) ## met Wolfe 2

          ## failed Wolfe 2 derivative not increased enough
          if (trial$dscore*(hi$alpha-lo$alpha)>=0) {
              hi <- lo }  
          lo <- trial 
       }  
    } ## end while(TRUE)
    return(NULL) ## failed
  } ## end zoom


  y <- G$y;  X <- G$X;  S <- G$S; not.exp <- G$not.exp;
  AR1.rho <- G$AR1.rho
  attr(X,"dimnames") <- NULL
  q0 <- G$q0; q.f <- G$q.f
  p.ident <- G$p.ident; n.terms <- G$n.terms
  family <- G$family; intercept <- G$intercept; offset <- G$offset;
  weights <- G$weights;  
  n <- nobs <- NROW(y)
  q <- ncol(X)
  dg <- fix.family.link(family)
  dv <- fix.family.var(family)
  nvars <- NCOL(X)
  EMPTY <- nvars == 0
  variance <- family$variance
  linkinv <- family$linkinv
  if (!is.function(variance) || !is.function(linkinv)) 
          stop("'family' argument seems not to be a valid family object")
  dev.resids <- family$dev.resids
  aic <- family$aic
  mu.eta <- family$mu.eta
  mu.eta <- family$mu.eta
  
   if (AR1.rho!=0) {
        ld <- 1/sqrt(1-AR1.rho^2) ## leading diagonal of root inverse correlation
        sd <- -AR1.rho*ld         ## sub diagonal
        row <- c(1,rep(1:nobs,rep(2,nobs))[-c(1,2*nobs)])
        weight.r <- c(1,rep(c(sd,ld),nobs-1))
        end <- c(1,1:(nobs-1)*2+1) 
        if (!is.null(G$mf$"(AR.start)")) { ## need to correct the start of new AR sections...
                 ii <- which(G$mf$"(AR.start)"==TRUE)
                 if (length(ii)>0) {
                    if (ii[1]==1) ii <- ii[-1] ## first observation does not need any correction
                    weight.r[ii*2-2] <- 0 ## zero sub diagonal
                    weight.r[ii*2-1] <- 1 ## set leading diagonal to 1
                 }
        }
        ## apply transform...
        X <- rwMatrix(end,row,weight.r,X)
        y <- rwMatrix(end,row,weight.r,y)    
      }   

  if (!is.function(variance) || !is.function(linkinv)) 
        stop("illegal `family' argument")
  valideta <- family$valideta
  if (is.null(valideta)) 
        valideta <- function(eta) TRUE
  validmu <- family$validmu
  if (is.null(validmu)) 
        validmu <- function(mu) TRUE

  if (is.null(mustart)) {
        eval(family$initialize)
  }
  else {
        mukeep <- mustart
        eval(family$initialize)
        mustart <- mukeep
  }

   ## Added code...
   if (family$family=="gaussian"&&family$link=="identity") strictly.additive <- TRUE else
      strictly.additive <- FALSE
   ## end of added code

  if (EMPTY) {
      eta <- rep.int(0, nobs) + offset
      if (!valideta(eta)) 
            stop("Invalid linear predictor values in empty model")
      mu <- linkinv(eta)

      if (AR1.rho!=0) {
        mu <- rwMatrix(end,row,weight.r,mu)  
      }   

      if (!validmu(mu)) 
          stop("Invalid fitted means in empty model")
      dev <- sum(dev.resids(y, mu, weights))
      w <- (weights * mu.eta(eta)^2)/variance(mu) ## incorrect for Newton
      residuals <- (y - mu)/mu.eta(eta)
      good <- rep(TRUE, length(residuals))
      boundary <- conv <- TRUE
      coef <- numeric(0)
      iter <- 0
      V <- variance(mu)
      alpha <- dev
      trA <- 0
    ##  GCV <- nobs * alpha/(nobs - gamma * trA)^2
    ##  UBRE <- alpha/nobs + 2 * gamma* trA*scale/n - scale 
      scale.est <- alpha/(nobs - trA)
      aic.model <- aic(y, n, mu, weights, dev) +  2 * trA
  } ### end if (EMPTY)
  else {
      eta <- if (!is.null(etastart)) 
          etastart
            else family$linkfun(mustart)
      mu <- as.numeric(linkinv(eta))
    # if (!(validmu(mu) && valideta(eta))) 
    #        stop("Can't find valid starting values: please specify some")
      S.t <- matrix(0,q,q) # define the total sum of the penalties times sp
      n.pen <- length(S) # number of penalties 
      if (length(sp)!=n.pen) stop (paste("length of sp has to be equal to", n.pen))
      if (n.pen>0) for (j in 1:n.pen) S.t <- S.t + sp[j]*S[[j]]
      # get sqrt of total penalty matrix...
      er <- eigen(S.t,symmetric=TRUE);  er$values[er$values<0] <- 0
      rS <- crossprod(sqrt(sqrt(er$values))*t(er$vectors))
     # ii <- p.ident==1 ## set to TRUE/FALSE
     ## count <- sum(p.ident) 
     ## iv <- array(0, dim=c(count,1)) # define an index vector for the coeff-s to be exponentiated
      iv <- (1:q)[p.ident] # define an index vector for the coeff-s to be exponentiated

      ##############################################
      ## Initialization of parameters starts here... 
      beta0 <- get("start",envir=env)
      dbeta0 <- get("dbeta.start",envir=env)
      sp.old <- get("sp.last",envir=env)
      if (length(beta0)==0) {
          # list argument to pcls for initializing model coefficients
          M <- list(X=X,p=rep(0.1,q),C=matrix(0,0,0),sp=sp,y=eta-offset,w=y*0+1) 
          M$Ain <- matrix(0,q,q); diag(M$Ain) <- rep(1,q);
          M$bin <- rep(-1e+12,q); M$bin[iv] <- 1e-12
          M$off <- rep(0,n.pen); M$S <- list()
          if (n.pen>0) for (j in 1:n.pen) {M$S[[j]] <- matrix(0,q,q); M$S[[j]] <- S[[j]]}
          beta.t <- pcls(M)      # initialize model coefficients (re-parameterized beta)
          beta <- beta.t         # initialize beta
          beta[iv] <- log(beta.t[iv]) # values of beta of the constrained terms
      }
      else {
          beta <- beta0
          beta.t <- beta    ## current beta tilde
          beta.t[iv] <- if (!not.exp) exp(beta[iv]) else notExp(beta[iv]) # values of re-para beta of the constrained term
      }
      ## end of initialization of parameters  
      ##############################################

      eta <- as.numeric(X%*%beta.t + as.numeric(offset))  # define initial linear predictor
      mu <- linkinv(eta)  # define initial fitted model
    ##  dev <- sum(dev.resids(y,mu,weights)) # define initial norm/deviance
    ##  pdev <- dev + sum((rS%*%beta)^2) # define initial penalized deviance 
    ##  old.pdev <- pdev       # initialize convergence control for penalized deviance
      b <- score(y=y,X=X,beta=beta,offset=offset,weights=weights,rS=rS,S.t=S.t,q=q,iv=iv,
          linkinv=linkinv,dev.resids=dev.resids,mu.eta=mu.eta,variance=variance)
      pdev <- b$pdev ## define initial penalized deviance 
      old.pdev <- pdev  ## initialize convergence control for penalized deviance
      grad <- b$grad ## initialize gradient of the penalized deviance

      score.scale <- 1+abs(b$pdev)
      ##################################################################
      ## added code here made on May 6, 2020 for scam version 1-2-6,
      ## following Simon N Wood gam.fit3(mgcv_1.8-31))....
      ################################################################### 
        ## shrink towards null.coef if immediately invalid 
      ##  betaold <- null.coef
      ##  etaold <- null.eta <- as.numeric(X%*%null.coef + as.numeric(offset))
      ##  old.pdev <- sum(dev.resids(y, linkinv(null.eta), weights)) + sum((rS%*%null.coef)^2) 
      ##  null.eta <- as.numeric(X%*%null.coef + as.numeric(offset))    
        ii <- 0
        while (!(validmu(mu) && valideta(eta))) { ## shrink towards null.coef if immediately invalid
          ii <- ii + 1
          if (ii>20) stop("Can't find valid starting values: please specify some")
          beta <- beta * .9 + null.coef * .1
          beta.t <- beta    ## current beta tilde
          beta.t[iv] <- if (!not.exp) exp(beta[iv]) else notExp(beta[iv])
          eta <- as.numeric(X%*%beta.t + offset)   
          mu <- linkinv(eta)
        }
        betaold <- null.coef <- beta
        betaold.t <- beta   
        betaold.t[iv] <- if (!not.exp) exp(betaold[iv]) else notExp(betaold[iv])
        etaold <- as.numeric(X%*%betaold.t + offset) 
        old.pdev <- sum(dev.resids(y,linkinv(etaold),weights)) + sum((rS%*%betaold)^2)
      ##################################################################
      
      ##   E <- matrix(0,q,q)   # define diagonal matrix E- second term of the Hessian, NO NEED
   ## CHECK if var-s below are needed!!! 
      boundary <- conv <- FALSE 
      old.warn <- getOption("warn")
      if (!control$print.warn) curr.warn <- -1 else curr.warn <- old.warn

      ############################################
      ### BFGS ALGORITHM STARTS HERE...
      initial.beta <- ibeta <- beta
      initial <- list(alpha = 0,mustart=b$fitted.values,start=coef(b))
      initial$pdev <- pdev; initial$grad <- grad

      B <- diag(length(initial$grad)) ## initial Hessian

      feps <- 1e-4
      for (i in 1:length(beta)) { ## loop to FD for Hessian
         ibeta <- beta
         ibeta[i] <- ibeta[i] + feps 
         b <- score(y=y,X=X,beta=ibeta,offset=offset,weights=weights,rS=rS,S.t=S.t,q=q,iv=iv,
            linkinv=linkinv,dev.resids=dev.resids,mu.eta=mu.eta,variance=variance)
         grad1 <- b$grad ## gradient of the penalized deviance
         B[i,] <- (grad1-grad)/feps 
         rm(b)
       } ## end of FD Hessian loop
       ## force initial Hessian to +ve def and invert... 
       B <- (B+t(B))/2
       eb <- eigen(B,symmetric=TRUE)
       eb$values <- abs(eb$values)
       thresh <- max(eb$values) * 1e-4
       eb$values[eb$values<thresh] <- thresh
       B <- eb$vectors%*%(t(eb$vectors)/eb$values) ## inverse of B
       ibeta <- as.numeric(beta)
       max.step <- 200

       c1 <- 1e-4;c2 <- .9 ## Wolfe condition constants

       pdev.hist <- rep(NA,max.step+1)
       pdev.hist[1] <- pdev

       check.derivs <- FALSE; ## actually I have not implemented checking derivatives here
       eps <- 1e-5

       uconv.ind <- rep(TRUE,ncol(B))
       rolled.back <- FALSE 

       ## MAIN BFGS LOOP...
      for (i in 1:max.step) { ## the main BFGS loop
   
        ## get the trial step ...
        step <- initial$grad*0
        step[uconv.ind] <- -B[uconv.ind,uconv.ind]%*%initial$grad[uconv.ind]

        ## following tends to have lower directional grad than above (or full version commented out below)
        #step <- -drop(B%*%initial$grad)
        ## following line would mess up conditions under which Wolfe guarantees update,
        ## *if* based only on grad and not grad and hess...  
        #step[!uconv.ind] <- 0 ## don't move if apparently converged 
    
        if (sum(step*initial$grad)>=0) { ## step not descending!
          ## Following would really be in the positive definite space... 
          ##step[uconv.ind] <- -solve(chol2inv(chol(B))[uconv.ind,uconv.ind],initial$grad[uconv.ind])
          step <- -diag(B)*initial$grad ## simple scaled steepest descent 
          step[!uconv.ind] <- 0 ## don't move if apparently converged 
        }

        ms <- max(abs(step))
        trial <- list()
        if (ms>maxNstep) { 
          trial$alpha <- maxNstep/ms
          alpha.max <- trial$alpha*1.05
          ## step <- maxNstep * step/ms
          #alpha.max <- 1 ## was 50 in place of 1 here and below
        } else {
             trial$alpha <- 1 
             alpha.max <- min(2,maxNstep/ms) ## 1*maxNstep/ms
        }
        initial$dscore <- sum(step*initial$grad)
        prev <- initial

        deriv <- 1 ## only get derivatives immediately for initial step length   
        while(TRUE) { ## step length control Alg 3.5 of N&W (2006, p60)
          beta <- as.numeric(ibeta) + trial$alpha*step
          b <- score(y=y,X=X,beta=beta,offset=offset,weights=weights,rS=rS,S.t=S.t,q=q,iv=iv,
            linkinv=linkinv,dev.resids=dev.resids,mu.eta=mu.eta,variance=variance)
          trial$pdev <- b$pdev
     
          if (deriv>0) {
              trial$grad <- b$grad
              trial$dscore <- sum(trial$grad*step)
              deriv <- 0 
           } else trial$grad <- trial$dscore <- NULL
      
          # rm(b)
           Wolfe2 <- TRUE
           ## check the first Wolfe condition (sufficient decrease)...
           if (trial$pdev > initial$pdev+c1*trial$alpha*initial$dscore||(deriv==0&&trial$pdev>=prev$pdev)) {
             trial <- zoom(prev,trial) ## Wolfe 1 not met so backtracking
             break
           } 

           if (is.null(trial$dscore)) { ## getting gradients
                 b <- score(y=y,X=X,beta=beta,offset=offset,weights=weights,rS=rS,S.t=S.t,q=q,iv=iv,
                       linkinv=linkinv,dev.resids=dev.resids,mu.eta=mu.eta,variance=variance)
                 trial$grad <- b$grad
                 trial$dscore <- sum(trial$grad*step)
            }
      
            ## Note that written this way so that we can pass on to next test when appropriate...
     
            if (abs(trial$dscore) <= -c2*initial$dscore) break; ## `trial' is ok. (2nd Wolfe condition met).
            Wolfe2 <- FALSE

            if (trial$dscore>=0) { ## increase at end of trial step
              trial <- zoom(trial,prev)
              Wolfe2 <- if (is.null(trial)) FALSE else TRUE
              break
            }
      
            prev <- trial
            if (trial$alpha == alpha.max) break ## { trial <- NULL;break;} ## step failed
            trial <- list(alpha = min(prev$alpha*1.3, alpha.max)) ## increase trial step to try to meet Wolfe 2
        } ## end of while(TRUE)

        ## Now `trial' contains a suitable step, or is NULL on complete failure to meet Wolfe,
        ## or contains a step that fails to meet Wolfe2, so that B can not be updated  
        if (is.null(trial)) { ## step failed
          beta <- ibeta
          if (rolled.back) break ## failed to move, so nothing more can be done.
          ## check for working infinite beta params... 
          uconv.ind <- abs(initial$grad) > score.scale*conv.tol*.1 
          if (sum(!uconv.ind)==0) break ## nothing to move back so nothing more can be done.
          trial <- initial ## reset to allow roll back
          converged <- TRUE ## only to signal that roll back should be tried
        } else { ## update the Hessian etc...
     
          yg <- trial$grad-initial$grad
          step <- step*trial$alpha
          rho <- sum(yg*step)
          if (rho>0) { #Wolfe2) { ## only update if Wolfe2 is met, otherwise B can fail to be +ve def.
            if (i==1) { ## initial step --- adjust Hessian as p143 of N&W
              B <- B * trial$alpha ## this is Simon's version 
              ## B <- B * sum(yg*step)/sum(yg*yg) ## this is N&W
            }
            rho <- 1/rho # sum(yg*step)
            B <- B - rho*step%*%(t(yg)%*%B)

            ## Note that Wolfe 2 guarantees that rho>0 and updated B is 
            ## +ve definite ...
            B <- B - rho*(B%*%yg)%*%t(step) + rho*step%*%t(step)
          }

          pdev.hist[i+1] <- trial$pdev

          beta <- ibeta <- ibeta + step 

          ## test for convergence
          converged <- TRUE
          score.scale <- .1+abs(trial$pdev)  ##abs(trial$scale.est) + abs(trial$score)   
          uconv.ind <- abs(trial$grad) > score.scale*conv.tol 
          if (sum(uconv.ind)) converged <- FALSE
          ## following must be tighter than convergence...
          uconv.ind <- abs(trial$grad) > score.scale*conv.tol*.1 
          if (abs(initial$pdev-trial$pdev) > score.scale*conv.tol) { 
            if (!sum(uconv.ind)) uconv.ind <- uconv.ind | TRUE ## otherwise can't progress
            converged <- FALSE      
          }
        } ## end of else { update the Hessian etc...

        ## roll back any `infinite' beta parameters to the point at
        ## which pen.deviance carries some information about them and continue 
        ## optimization. Guards against early long steps missing shallow minimum. 
        if (converged) { ## try roll back for `working inf' betas...
            if (sum(!uconv.ind)==0||rolled.back) break
            rolled.back <- TRUE
            counter <- 0
            uconv.ind0 <- uconv.ind 
            while (sum(!uconv.ind0)>0&&counter<5) {
              ## shrink towards initial values...
              beta[!uconv.ind0] <- beta[!uconv.ind0]*.8 + initial.beta[!uconv.ind0]*.2
              b <- score(y=y,X=X,beta=beta,offset=offset,weights=weights,rS=rS,S.t=S.t,q=q,iv=iv,
                      linkinv=linkinv,dev.resids=dev.resids,mu.eta=mu.eta,variance=variance)
              trial$pdev <- b$pdev
              trial$grad <- b$grad
              trial$dscore <- sum(trial$grad*step)
              rm(b);counter <- counter + 1
              ## note that following rolls back until there is clear signal in derivs...
              uconv.ind0 <- abs(trial$grad) > score.scale*conv.tol*20        
              uconv.ind0 <- uconv.ind0 | uconv.ind ## make sure we don't start rolling back unproblematic betas 
            }
            uconv.ind <- uconv.ind | TRUE
            ibeta <- as.numeric(beta)
        }
    
        initial <- trial
        initial$alpha <- 0
      
      } ## end of iteration loop


       if (is.null(trial)) { 
         ct <- "step failed"
         beta <- ibeta
         trial <- initial
       }
       else if (i==max.step) ct <- "iteration limit reached" 
       else ct <- "full convergence"
       ## final pen. deviance calculation
       b <- score(y=y,X=X,beta=beta,offset=offset,weights=weights,rS=rS,S.t=S.t,q=q,iv=iv,
                      linkinv=linkinv,dev.resids=dev.resids,mu.eta=mu.eta,variance=variance)
       pdev <- b$pdev;
       grad <- b$grad
       ## get approximate Hessian and its inverse...
       ev <- eigen(B,symmetric=TRUE)
       ind <- ev$values>max(ev$values)*.Machine$double.eps^.9
       ev$values[ind] <- 1/ev$values[ind]
       ev$values[!ind] <- 0
       inv.B <- ev$vectors %*% (ev$values*t(ev$vectors)) ## inverse of B
 
       ## now define things at their converged values from the quasi-Newton method...
       beta.t <- beta
       beta.t[iv] <- if (!not.exp) exp(beta[iv]) else notExp(beta[iv]) ## values of re-para beta of the shape constrained terms     
       eta <- as.numeric(X%*%beta.t + offset)    
       mu <- linkinv(eta)   
       dev <- sum(dev.resids(y,mu,weights))
   
       Cdiag <- rep(1,q)
       Cdiag[iv] <- if (!not.exp) beta.t[iv] else DnotExp(beta[iv])
       X1 <- t(Cdiag*t(X)) 
       g.deriv <- 1/mu.eta(eta)   ## diag(G)
       w1 <- weights/(variance(mu)*g.deriv^2)   ## diag(W1)
       alpha <- 1+(y-mu)*(dv$dvar(mu)/variance(mu)+dg$d2link(mu)/g.deriv) ## alpha elements of W
       w <- w1*alpha    ## diag(W)
       abs.w <- abs(w)  ## absolute values of the diag(W) (maybe no need to output?)
  
       normgrad.Dp <- max(abs(grad)) 

       ## calculating edf and tr(A)...
       ## ISSUE: how to get approximate expected Hessian of unpenalized deviance ??????
       ## calculation of the edf below is not correct
       ## B.unpen <- B-S.t ## approximate Hessian of unpenalized deviance 
       ##  ## force approx Hessian to be +ve def ... 
       ## B.unpen <- (B.unpen+t(B.unpen))/2
       ## ebu <- eigen(B.unpen,symmetric=TRUE)
       ## ebu2 <- ifelse(ebu$values < 0, 0, ebu$values)
       ## B.unpen <- ebu$vectors%*%(t(ebu$vectors)*ebu2)
       ## B.unpen <- B.unpen/sqrt(tcrossprod(diag(B.unpen)))
       ## edf <- rowSums(inv.B*t(B.unpen)) ## diagonal elements of inv.B%*%B.unpen
       ## trA <- sum(edf)

       ## calculating edf and tr(A) using analytical Hessian (taking 2nd order deriv of the log likelihood)
       ## as it is in scam.fit()...
       ## tr(A), dev, sig2 are needed for sp selection for efsudr.scam2()...
       I.plus <- rep(1,nobs)   # define diagonal elements of the matrix I^{+}
       I.plus[w<0] <- -1
    ##   L <- c(1/alpha)    # define diagonal elements of L=diag(1/alpha)
       ## question: should set alpha=1 so use expected Hessian?
       L <- 1
       wX1<-sqrt(w1)[1:nobs]*X1
       wX11<-rbind(wX1,rS)
       Q <- qr(wX11,LAPACK=TRUE) 
       R <- qr.R(Q)
       rp <- 1:ncol(R)
       rp[Q$pivot] <- rp ## reverse pivot X=Q%*%R[,rp]
       R.out <- R[,rp]  ## unpivoted R, needed for summary function
       rank <- Rrank(R)
       if (rank==ncol(R)) { ## no need to truncate, can just use QR
             P <- backsolve(R,diag(ncol(R)))[rp,]
             K <- qr.Q(Q)[1:nobs,]
       } else { ## need SVD step
             R <- R[,rp] ## unpivoted R
             s1 <- svd(R)
             d.inv1 <- rep(0,q)
             good1 <- s1$d >= max(s1$d)*.Machine$double.eps^.5
             d.inv1[good1] <- 1/s1$d[good1]
             P <- t(d.inv1*t(s1$v))
             K <- qr.qy(Q,rbind(s1$u,matrix(0,nobs,q)))[1:nobs,]
       }
       ## matrices wX1, P are needed for efsudr.scam2() for calculating derivatives 
       ## of deviance and tr(A) w.r.t. log(sp) ...
       KtILQ1R <- crossprod(L*I.plus*K,wX1) ## t(L*I.plus*K)%*%wX1 
       edf <- rowSums(P*t(KtILQ1R))
       trA <- sum(edf)

       scale.est <- dev/(nobs-trA)  #  scale estimate
       ## re-transforming 'mu' back on original scale and hence residuals, in case of correlated errors...
       residuals <- rep.int(NA, nobs)
       residuals <- (y-mu)*g.deriv

       ## some more of the return values....
       dlink.mu  <- 1/mu.eta(eta); Var<- variance(mu)
    #   link <- family$linkfun(mu); d2link.mu <- dg$d2link(mu)
    #   dvar.mu <- dv$dvar(mu); d2var.mu <- dv$d2var(mu)
    #   d3link.mu <- dg$d3link(mu)
    #   z <- g.deriv*(y-mu)+X1%*%beta
       #############################
   
       ## derivatives of beta wrt rho[j] by the Implicit Function Theorem...
       dbeta.rho <- matrix(0,q,n.pen) # define matrix of the parameters derivatives
       if (n.pen>0) for (j in 1:n.pen) {
                     dbeta.rho[,j] <- -sp[j]*inv.B%*%(S[[j]]%*%beta) 
                    }
     
       aic.model <- aic(y, n, mu, weights, dev) +  2*trA
       if (AR1.rho!=0) { ## correct aic for AR1 transform
          df <- 1 ## if (getARs) sum(b$model$"(AR.start)") else 1
          aic.model <- aic.model - 2*(n-df)*log(ld) 
       }
       assign("start",beta,envir=env)
       assign("dbeta.start",dbeta.rho,envir=env)
       assign("sp.last",sp,envir=env)
  } ### end if (!EMPTY) 

  #list(L=L,C1diag=C1diag,E=E,iter=iter, 
  #    P=P,K=K, C2diag=C2diag, KtILQ1R= KtILQ1R, KtIQ1R=KtIQ1R, 
  #    dlink.mu=dlink.mu,Var=Var, abs.w=drop(abs.w),
  #    link=link,w=as.numeric(w),w1=drop(w1),d2link.mu=d2link.mu,I.plus=I.plus,
  #    dvar.mu=dvar.mu,d2var.mu=d2var.mu,
  #    ok1=ok1,alpha=as.numeric(alpha),d3link.mu=d3link.mu,eta=eta,iter=iter,
  #    Dp.gnorm=Dp.gnorm, Dp.g=Dp.g,d=d, conv=conv, illcond=illcond,R=R.out, edf=edf,trA=trA,
  #    residuals=residuals,z=z,dbeta.rho=dbeta.rho, aic=aic.model,rank=rank) 

 list(pdev=pdev,grad.Dp=grad,iter=i,conv =ct,old.beta=ibeta, ## B=B,inv.B=inv.B, 
      gcv=dev*nobs/(nobs-gamma*trA)^2, sp=sp, Var=Var,dlink.mu=dlink.mu, 
      mu=mu,X=G$X,y=drop(G$y), X1=X1,beta=beta,beta.t=beta.t,iv=iv,S=S,S.t=S.t,rS=rS, 
      P=P,K=K,L=L,wX1=wX1, rank=rank, F=P%*%(KtILQ1R), w=as.numeric(w),w1=drop(w1),
      deviance=dev,scale.est=scale.est, normgrad.Dp=normgrad.Dp,edf=edf,trA=trA,
      residuals=residuals,dbeta.rho=dbeta.rho, aic=aic.model,pdev.hist=pdev.hist[!is.na(pdev.hist)]
 )   
} ## end of scam.fit1


####################################################################################
## function to get null deviance and covariance matrices after quasi-Newton fit   ##
####################################################################################

scam.fit.post1<- function(G, object) ##,sig2,offset,intercept, weights,scale.known, not.exp)
{  ## Function to compute null deviance and covariance matrices after a scam fit by scam.fit1().
   ## object - object from estimate.scam()
   y <- G$y; X <- G$X;
   sig2 <- G$sig2; offset <- G$offset; intercept <- G$intercept; 
   weights <- G$weights; scale.known <- G$scale.known
   nobs <- NROW(y) # number of observations
   q <- ncol(X)
   if (G$AR1.rho!=0) {
        ld <- 1/sqrt(1-G$AR1.rho^2) ## leading diagonal of root inverse correlation
        sd <- -G$AR1.rho*ld         ## sub diagonal
        row <- c(1,rep(1:nobs,rep(2,nobs))[-c(1,2*nobs)])
        weight.r <- c(1,rep(c(sd,ld),nobs-1))
        end <- c(1,1:(nobs-1)*2+1) 
        if (!is.null(G$AR.start)) { ## need to correct the start of new AR sections...
                 ii <- which(G$AR.start==TRUE)
                 if (length(ii)>0) {
                    if (ii[1]==1) ii <- ii[-1] ## first observation does not need any correction
                    weight.r[ii*2-2] <- 0 ## zero sub diagonal
                    weight.r[ii*2-1] <- 1 ## set leading diagonal to 1
                 }
        }
        ## apply transform...
        X <- rwMatrix(end,row,weight.r,X)
        y <- rwMatrix(end,row,weight.r,y)    
      }   
  
   linkinv <- object$family$linkinv
   dev.resids <- object$family$dev.resids
   dg <- fix.family.link(object$family)
   dv <- fix.family.var(object$family)

   eta <- as.numeric(X%*%object$beta.t + offset) # linear predictor
   mu <- linkinv(eta)          # fitted values
   dev <- sum(dev.resids(y,mu,weights)) # deviance of the final model
   
   wtdmu <- if (intercept) sum(weights * G$y)/sum(weights) 
              else linkinv(offset)

   null.dev <- sum(dev.resids(G$y, wtdmu, weights))
 
   n.ok <- nobs - sum(weights == 0)
   nulldf <- n.ok - as.integer(intercept)
   
   edf <- object$edf ## effective degrees of freedom
   trA <- sum(edf)
   edf1 <- 2*edf - rowSums(t(object$F)*object$F) ## alternative form (needed for reference degrees of freedom used in computing test statistic and the p-values, but since the null distributions are non-standard the reference df is not very interpretable????)
    
   ## calculating the approximate covariance matrices 
   ## (dealing with the expected Hessian of the log likelihood) ...
   ## get the inverse of the expected Hessian...
   if (!scale.known) sig2 <-  dev/(nobs-trA)  # scale estimate
   Vb <- tcrossprod(object$P) * sig2 
          ## P%*%t(P)*sig2 # Bayesian posterior covariance matrix for the parameters 
   Ve <- crossprod(object$K%*%t(object$P)) *sig2
        #PKt%*%t(PKt)*sig2 # covariance matrix of the parameter estimators 
   ## Delta method to get covariance matrix for the reparametrized parameters...
   df.p <- rep(1,q)
   df.p[object$iv] <- object$beta.t[object$iv]
   Vb.t <- t(df.p*t(df.p*Vb))
   Ve.t <- t(df.p*t(df.p*Ve))
 
   residuals <- rep.int(NA, nobs)
   g.deriv <- 1/object$family$mu.eta(eta) # diag(G)
   residuals <- (G$y-mu)*g.deriv # the working residuals for the fitted model
      
   aic.model <- object$family$aic(y, nobs, mu, weights, dev) + 2*trA
   if (G$AR1.rho!=0) { ## correct aic for AR1 transform
        df <- 1 ## if (getARs) sum(G$AR.start) else 1
        aic.model <- aic.model - 2*(nobs-df)*log(1/sqrt(1-G$AR1.rho^2)) 
   }
    
   list (null.dev=null.dev, df.null=nulldf,Vb=Vb,Vb.t=Vb.t,Ve=Ve,Ve.t=Ve.t,rank=object$rank,
        sig2=sig2,edf=edf,edf1=edf1, trA=trA, deviance=dev,residuals=residuals, aic=aic.model, mu=mu, eta=eta)
} ## end of scam.fit.post1





