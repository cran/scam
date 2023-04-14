## (c) Natalya Pya (2012-2023). Released under GPL2.
## efsudr.scam based on (c) Simon N Wood (efsudr(mgcv))


#########################################################
# Function to return gcv/ubre ...                      ##
#########################################################

gcv.ubre <- function(rho,G,env,control) 
{  ## function to get GCV.UBRE value for optim()...
   if (length(rho)!= length(G$off)) stop (paste("length of rho and n.terms has to be the same"))
   sp <- exp(rho)
   b <- scam.fit(G=G, sp=sp, env=env,control=scam.control()) 
   if (G$scale.known) #  value of Mallow's Cp/UBRE/AIC ....
      {  n <- nrow(G$X)
         gcv.ubre <- b$dev/n - G$sig2 +2*G$gamma*b$trA*G$sig2/n
      }  else   # value of GCV ...
           gcv.ubre <- b$gcv
   return(gcv.ubre)
}

#########################################################
## function to get the gradient of the gcv/ubre.....   ##
#########################################################

gcv.ubre.derivative <- function(rho,G, env, control)  
{  ## function to return derivative of GCV or UBRE for optim...
   gcv.ubre_grad(rho, G, env, control=control)$gcv.ubre.rho
}

#############################################################################
## for nlm() function to get the gcv/ubre and gradient of the gcv/ubre.....##
#############################################################################

dgcv.ubre.nlm <- function(rho,G, env,  control) 
{  ## GCV UBRE objective function for nlm
   gg <- gcv.ubre_grad(rho, G, env, control=control) 
   attr(gg$gcv.ubre,"gradient") <- gg$gcv.ubre.rho
   gg$gcv.ubre
}

#######################################################
#### estimate.scam()....                             ##
#######################################################

estimate.scam <- function(G,optimizer,optim.method,rho, env,  control)
                   ##  check.analytical, del, devtol.fit, steptol.fit)
## function to select smoothing parameter by minimizing GCV/UBRE criterion...

{  ## set 'newton' method for coeff estimation if not specified 
   ## (it can happen if 'optimizer' was supplied with one element, specifying 
   ## the sp optimization method to use as, e.g., 'optimizer="efs")
   
   if (is.na(optimizer[2])) 
         optimizer[2] <- "newton"

   if (!(optimizer[2] %in% c("newton", "bfgs")) )
          stop("unknown optimization method to use for the model coefficient estimation")
   if (optimizer[2]== "bfgs")
        if (optimizer[1] !="efs") {
           warning("`bfgs` method for the coefficient estimation works only together with `efs` method; `efs' was used")
           optimizer[1] <- "efs"
        }

   if (!(optimizer[1] %in% c("bfgs", "nlm", "optim","nlm.fd","efs")) )
          stop("unknown optimization method to use to optimize the smoothing parameter estimation criterion")

   if (length(rho)==0) { ## no sp estimation to do -- run a fit instead
     optimizer[1] <- "no.sps" ## will cause scam.fit/scam.fit1 to be called, below
   }

   if (optimizer[1] == "bfgs") {## minimize GCV/UBRE by BFGS...
         b <- bfgs_gcv.ubre(gcv.ubre_grad,rho=rho, G=G,env=env, control=control) ## check.analytical=check.analytical, del=del,devtol.fit=devtol.fit, steptol.fit=steptol.fit) 
         sp <- exp(b$rho)
         object <- b$object
         object$gcv.ubre <- b$gcv.ubre
         object$dgcv.ubre <- b$dgcv.ubre
         object$termcode <- b$termcode
         object$check.grad <- b$check.grad
         object$dgcv.ubre.check <- b$dgcv.ubre.check
         object$conv.bfgs <- b$conv.bfgs
         object$iterations <- b$iterations
         object$score.hist <- b$score.hist
   }
   else if (optimizer[1]=="optim"){  ## gr=gcv.ubre.derivative
              if (!(optim.method[1] %in% c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN")))
                 {   warning("unknown optim() method, `L-BFGS-B' were used")
                     optim.method[1] <- "L-BFGS-B"
                 }
              if (is.na(optim.method[2])) 
                 {   warning("the second parameter of optim.method argument is not supplied, 
                      finite-difference approximation of the gradient were used")
                     grr <- NULL
                 }  else if (!(optim.method[2] %in% c("fd","grad")))
                      {   warning("only `fd' and `grad' options are possible, finite-difference 
                              approximation of the gradient were used")
                          grr <- NULL
                      }  else if (optim.method[2] == "grad")
                                 grr <- gcv.ubre.derivative
                              else 
                                 grr <- NULL
              b <- optim(par=rho,fn=gcv.ubre, gr=grr, method=optim.method[1],G=G, control=
         list(factr=control$optim$factr,lmm=min(5,length(rho))), env=env) 
              sp <- exp(b$par)
              gcv.ubre <- b$value
              dgcv.ubre <- NULL
              iterations <- b$counts
              termcode <- b$convergence
              if (termcode == 0)
                     conv <- "Successful completion"
              else if (termcode == 1)  
                      conv <- "The iteration limit `maxit' had been reached"
              else if (termcode == 10)  
                      conv <- "Degeneracy of the Nelder-Mead simplex"
              else if (termcode == 51)  
                      conv <- "A warning from the `L-BFGS-B' method; see help for `optim' for further details"
              else if (termcode == 52)  
                      conv <- "An error from the `L-BFGS-B' method; see help for `optim' for further details"
   }
   else if (optimizer[1]=="nlm.fd") {## nlm() with finite difference derivatives...
             b <- nlm(f=gcv.ubre, p=rho,typsize=rho, stepmax = control$nlm$stepmax, ndigit = control$nlm$ndigit,
	               gradtol = control$nlm$gradtol, steptol = control$nlm$steptol, 
                       iterlim = control$nlm$iterlim,  G=G, env=env, control=control) 
   }
   else if (optimizer[1]=="nlm"){  ## nlm() with analytical derivatives...
            b <- nlm(f=dgcv.ubre.nlm, p=rho,typsize=rho, stepmax = control$nlm$stepmax, ndigit = control$nlm$ndigit, gradtol = control$nlm$gradtol, steptol = control$nlm$steptol, iterlim = control$nlm$iterlim, G=G,env=env, control=control) 
   }
   else if (optimizer[1]=="efs"){  ## Extended Fellner-Schall method
          ## if bfgs method is used for model coeff. estimation, the inner fit function is scam.fit1(),
          ## if newton method then the fit function is scam.fit()...
          fit.fn <- if (optimizer[2]== "bfgs") scam.fit1 else scam.fit
          b <- efsudr.scam2(fit.fn=fit.fn,G=G,lsp=rho,env=env, control=control)
          sp <- b$sp
          object <- b
          object$gcv.ubre <- b$gcv.ubre
          object$outer.info <- b$outer.info
          object$inner.info <- b$inner.info
         ## if (optimizer[2]== "bfgs") { ## bfgs method for the model coeff. estimation  
         ##     b <- efsudr.scam2(G=G,lsp=rho,env=env, control=control)
         ##     sp <- b$sp
         ##     object <- b
         ##   #  object$iterations <- b$niter
         ##   #  object$conv <- b$outer.info$conv
         ##   #  object$score.hist <- b$outer.info$score.hist
         ##     object$gcv.ubre <- b$gcv.ubre
         ##     object$outer.info <- b$outer.info
         ##     object$inner.info <- b$inner.info
         ## } else { ## newton method for the model coeff. estimation
         ##     b <- efsudr.scam(G=G,lsp=rho,env=env, control=control)
         ##   ##  b <- efsudr.scam1(G=G,lsp=rho,env=env, control=control)
         ##     sp <- b$sp
         ##     object <- b
         ##     object$gcv.ubre <- b$gcv.ubre
         ##     object$outer.info <- b$outer.info
         ##     object$inner.info <- b$inner.info
         ## }
   } 

   ## get some convergence info from the optimization method used...
   if (optimizer[1]== "nlm.fd" || optimizer[1]== "nlm"){ 
          sp <- exp(b$estimate)
          gcv.ubre <- b$minimum
          dgcv.ubre <- b$gradient
          iterations <- b$iterations
          termcode <- b$code
          if (termcode == 1)
                 conv <- "Relative gradient is close to zero, current iterate is probably solution"
          else if (termcode == 2)  
                 conv <- "Successive iterates within tolerance, current iterate is probably solution"
          else if (termcode == 3)  
                 conv <- "Last global step failed to locate a point lower than `estimate'. Either 
                       `estimate' is an approximate local minimum of the function or 
                        `steptol' is too small"
          else if (termcode == 4)  
                 conv <- "Iteration limit exceeded"
          else if (termcode == 5)  
                 conv <- "Maximum step size `stepmax' exceeded five consecutive 
                         times. Either the function is unbounded below, becomes asymptotic 
                         to a finite value from above in some direction or stepmax is too small"   
   }
   ## fit the model using the optimal sp from "optim" or "nlm"...
   if (optimizer[1]== "nlm.fd" || optimizer[1]== "nlm" || optimizer[1]== "optim"){
         object <- scam.fit(G=G, sp=sp,env=env,control=control) 
         object$conv <- object$conv
         object$gcv.ubre <- gcv.ubre
         object$dgcv.ubre <- dgcv.ubre 
         object$outer.info <- list(termcode=termcode,conv=conv,iterations=iterations)
   }
   else if (optimizer[1]=="no.sps"){
        # object <- scam.fit(G=G, sp=exp(rho),env=env,control=control)      
         sp <- G$sp    
         object <- if (optimizer[2] == "newton")  
                      scam.fit(G=G, sp=sp,env=env, control=control) 
                   else scam.fit1(G=G, sp=sp,env=env, control=control) ## BFGS optimization 
        # object$optimizer[1] <- "NA"      
   }

   if (optimizer[1]=="optim"){
         object$optim.method <- rep(NA,2)
         object$optim.method[1] <- optim.method[1]
         if (!is.null(grr))
              object$optim.method[2] <- "grad"
   }
   
   object$sp <- sp
   object$q.f <- G$q.f
   object$p.ident <- G$p.ident
   object$S <- G$S
   object$optimizer <- optimizer
   if (optimizer[1]=="no.sps") object$optimizer[1] <- "NA" 
   object
} ## estimate.scam



efsudr.scam2 <- function(fit.fn,G,lsp,env, control)
## Extended Fellner-Schall method for regular families  as in efsudr(mgcv), however,
## rather than minimizing REML (as in gam(mgcv)), the function minimizes GCV/UBRE criterion;
## also there is an alternative BFGS method in place of the full Newton.
## The dependence of H on lambda is neglected. 

## 'fit.fn' - the routine to estimate the model coefficients: 
##       fit.fn = scam.fit1 if a quasi-Newton/BFGS method is used for the coefficient estimation, 
##       fit.fn = scam.fit in case of the full Newton (PIRLS).
## This function is the same as efsudr.scam1 with the only difference that scam.fit() is either  
## scam.fit1() or scam.fit()...

{  ## deriv.dev.edf <- function(G,fit,sp){ ## the expected approx Hessian is not correct here, not +ve def!!
   ##  ## function to calculate derivatives of deviance and tr(A) w.r.t. log(sp)
   ## ## Hessian is replaced by a BFGS approximation from the scam.fit1,
   ## ## the dependence of H on sp/lambda is neglected
   ##  b <- fit
   ##  nsp <- length(G$S) 
   ##  H <- b$hess ## expected approximate Hessian 
   ##  Vb <- b$inv.B ## inverse of the approximate Hessian
   ##  trA.rho<- rep(0,nsp)
   ##  beta.rho <- matrix(0,ncol(G$X),nsp) ## matrix of the parameters derivatives wrt log(sp)
   ##  for (j in 1:nsp){
   ##        VbS <- Vb%*%b$S[[j]] 
   ##       trA.rho[j] <- -sp[j]*sum((VbS%*%Vb)*t(H)) ## -lambda*tr(Vb*Sj*Vb*H)
   ##        beta.rho[,j] <- -sp[j]*(VbS%*%b$beta) 
   ##     }
   ##  y.mu <- drop(b$y)-b$mu
   ##  c <- -2*y.mu/(b$Var*b$dlink.mu)
   ##  D.beta <- t(b$X1)%*%c  ## derivative of the deviance w.r.t. beta   
   ##  D.rho <- t(D.beta)%*%beta.rho  ## derivative of the deviance w.r.t. log(sp)
   ##  list(trA.rho=trA.rho,D.rho=D.rho)
   ## } ## deriv.dev.edf()

  deriv.dev.edf <- function(G,fit,sp){
  ## function to calculate derivatives of deviance and tr(A) w.r.t. log(sp),
  ## using expected Hessian and neglecting its dependence on sp
     b <- fit
     nsp <- length(G$S) 
     H <- crossprod(b$wX1) ## expected Hessian = CtXtWXC 
     Vb <- tcrossprod(b$P) ## P%*%t(P)*sig2 Bayesian posterior covariance matrix for the parameters
                 ## but without ? scale parameter here (inverse of the Hessian)
     trA.rho<- rep(0,nsp)
     beta.rho <- matrix(0,ncol(G$X),nsp) ## matrix of the parameters derivatives wrt log(sp)
    ## VbS <- matrix(0,dim(Vb)[1],dim(Vb)[2])
     for (j in 1:nsp){
           VbS <- Vb%*%b$S[[j]] 
           trA.rho[j] <- -sp[j]*sum((VbS%*%Vb)*t(H)) ## -lambda*tr(Vb*Sj*Vb*H)
           beta.rho[,j] <- -sp[j]*(VbS%*%b$beta) 
        }
     y.mu <- drop(b$y)-b$mu
     c <- -2*y.mu/(b$Var*b$dlink.mu)
     D.beta <- t(b$X1)%*%c  ## derivative of the deviance w.r.t. beta   
     D.rho <- t(D.beta)%*%beta.rho  ## derivative of the deviance w.r.t. log(sp)
     list(trA.rho=trA.rho,D.rho=D.rho)
  } ## deriv.dev.edf()


  crit.gcv <- function(fit,G,n){
  ## function to get the GCV criterion value
  ## gcv=dev*nobs/(nobs-gamma*trA)^2
    fit$gcv
  } 

  crit.ubre <- function(fit,G,n){
  ## function to get the UBRE criterion value
  ## ubre = dev/n - sig2 +2*gamma*trA*sig2/n
    ## cr <- fit$deviance +2*G$gamma*fit$trA*G$sig2
    cr <- fit$deviance/n +2*G$gamma*fit$trA*G$sig2/n - G$sig2
    cr
  } 

  if (G$scale.known) ## use UBRE...
         crit <- crit.ubre
  else crit <- crit.gcv ## use GCV
 
  nsp <- length(G$S) # number of smoothing parameters
  spind <- 1:nsp ## index of smoothing params in lsp
  lsp[spind] <- lsp[spind] + 2.5 
  mult <- 1
  fit <- fit.fn(G=G,sp=exp(lsp), env=env,control=control)
  q <- ncol(G$X)
  n <- nrow(G$X)
  score.hist <- rep(0,200)
   
  D.rho <- tau.rho <- rep(0,nsp) ## initialize derivatives of deviance and edf (tau/trA) wrt log(sp)
  for (iter in 1:200) { ## loop for GCV/UBRE minimization...
     ## calculation of the derivatives of deviance and edf wrt rho...
     der <-  deriv.dev.edf(G=G, fit=fit,sp=exp(lsp))
     D.rho <- der$D.rho
     tau.rho <- der$trA.rho
     
     if (G$scale.known)
        a <- pmax(0,-2*G$gamma*G$sig2*tau.rho)
     else a <- pmax(0,-2*G$gamma*fit$deviance*tau.rho/(n-fit$trA))
     r <- a/pmax(0,D.rho)
    # r[a==0 & D.rho==0] <- 1 
     r[a==0 | D.rho==0] <- 1 
     r[!is.finite(r)] <- 1e6
     lsp1 <- lsp
     lsp1 <- pmin(lsp + log(r)*mult,control$efs.lspmax)
     max.step <- max(abs(lsp1-lsp))
     old.gcv <- crit(fit, G,n) ## fit$gcv
     fit <- fit.fn(G=G,sp=exp(lsp1), env=env,control=control)
      ## some step length control...
     fit.gcv <- crit(fit, G,n) ## undated gcv/ubre value
   
     if (fit.gcv<=old.gcv) { ## improvement
       if (max.step<.05) { ## consider step extension (near optimum)
         lsp2 <- lsp
         lsp2 <- pmin(lsp + log(r)*mult*2,control$efs.lspmax) ## try extending step...
         fit2 <- fit.fn(G=G,sp=exp(lsp2), env=env,control=control)
         fit2.gcv <- crit(fit2,G,n)
         if (fit2.gcv < fit.gcv) { ## improvement - accept extension
           fit <- fit2;lsp <- lsp2
           fit.gcv <-  fit2.gcv 
	   mult <- mult * 2
          } else { ## accept old step
            lsp <- lsp1
           }
         } else lsp <- lsp1
     } else { ## no improvement 
      # while (fit$gcv > old.gcv&&mult>1) { ## don't contract below 1 as update doesn't have to improve GCV
      #     mult <- mult/2 ## contract step
      #	   lsp1 <- lsp
      #     lsp1 <- pmin(lsp + log(r)*mult,control$efs.lspmax)
      #	   fit <- fit.fn(G=G,sp=exp(lsp1), gamma=gamma, env=env,control=control)
      #  }
        gcv.thresh <- 10*(.1 +abs(old.gcv))*.Machine$double.eps^.5
        ii <- 1
        maxHalf.fit <- 15
        while (is.na(fit.gcv) || (fit.gcv-old.gcv) > gcv.thresh) { # 'step reduction' approach
             if (ii > maxHalf.fit) 
                 break ## stop ("step reduction failed")
             ii <- ii+1
             mult <- mult/2 ## contract step
             lsp1 <- lsp
             lsp1 <- pmin(lsp + log(r)*mult,control$efs.lspmax)
             fit <- fit.fn(G=G,sp=exp(lsp1), env=env,control=control)
             fit.gcv <- crit(fit,G,n)
         }
        lsp <- lsp1
        if (mult<1) mult <- 1
      }
     score.hist[iter] <- fit.gcv
     ## break if EFS step small and GCV change negligible over last 3 steps.
     if (iter>3 && max.step<.05 && max(abs(diff(score.hist[(iter-3):iter])))<control$efs.tol) break
     ## or break if deviance not changing...
     if (iter==1) old.dev <- fit$deviance else {
       if (abs(old.dev-fit$deviance) < 100*control$devtol.fit*abs(fit$deviance)) break
       old.dev <- fit$deviance
     }
  } # end of 'for' loop for GCV/UBRE minimization
 
  fit$sp <- exp(lsp)
  fit$niter <- iter
  fit$outer.info <- list(iter = iter,score.hist=score.hist[1:iter])
  fit$outer.info$conv <- if (iter==200) "iteration limit reached" else "full convergence"
  fit$gcv.ubre <- fit.gcv
  fit$inner.info <- list(iter=fit$iter,conv=fit$conv,pdev.hist=fit$pdev.hist) ## convergence info from BFGS for pen. deviance minimization 
  fit
} ## efsudr.scam2




##########################################################################
## Extended Fellner-Schall method for regular families as in efsudr(mgcv) (c) Simon N Wood
## with PIRLS performed by scam.fit. Rather than minimizing REML (as in gam(mgcv)), 
## the function minimizes GCV/UBRE criterion...
##########################################################################

efsudr.scam <- function(G,lsp,env, control){##maxit=200, devtol.fit=1e-7, steptol.fit=1e-7,
                                            # gamma=1, start=NULL, etastart=NULL, mustart=NULL, env=env){
## Extended Fellner-Schall method for regular families,
## with full-Newton/PIRLS method used for the coefficient estimation, performed by scam.fit().
## In this implementation I use the third order derivatives of the log likelihood when calculating
## the derivatives of the tau/edf. So I did not neglect the dependence of H on lambda 
## as it is in Wood and Fasiolo (2017)  

  deriv.dev.edf <- function(G,fit,sp){
  ## function to calculate derivatives of deviance and tr(A) w.r.t. log(sp)
     b <- fit
     nsp <- length(G$S) 
     q <- ncol(G$X)
     n <- nrow(G$X) 
     y.mu <- drop(b$y)-b$mu
     c <- -2*y.mu/(b$Var*b$dlink.mu)
     D.beta <- t(b$X1)%*%c  # derivative of the deviance w.r.t. beta 
     
     D.rho <- rep(0,nsp)    # define derivative of the deviance wrt rho
     D.rho <- t(D.beta)%*%b$dbeta.rho
     ## calculating the derivative of tr(A) w.r.t. log(sp)
     d2link.dlink <- b$d2link.mu/b$dlink.mu
     a1 <- as.numeric(y.mu*d2link.dlink)    # a constant for updating the derivative of z
     a2 <- as.numeric(b$w^2*(b$dvar.mu*b$dlink.mu+2*b$Var*b$d2link.mu)) # a constant for updating the derivative of w
     eta.rho <- matrix(0,n,nsp) # define derivatives of the linear predictor
     N_rho <- matrix(0,q,nsp) # define diagonal elements of N_j
     w.rho <- matrix(0,n,nsp) # define derivatives of the diagonal elements of W
     alpha.rho <- matrix(0,n,nsp) # define derivatives of alpha
     T_rho <- matrix(0,n,nsp) # define diagonal elements of T_j
     ## a constant for the derivative of the alpha
     dvar.var <- b$dvar.mu/b$Var
     alpha1 <- as.numeric(-(dvar.var+d2link.dlink)/b$dlink.mu -
                        y.mu*(dvar.var^2 + d2link.dlink^2 - b$d2var.mu/b$Var -
                                  b$d3link.mu/b$dlink.mu)/b$dlink.mu)
     eta.rho <- b$X1%*%b$dbeta.rho
     N_rho[b$iv,] <- b$dbeta.rho[b$iv,]
     alpha.rho <- alpha1*eta.rho
     w.rho <- -a2*b$alpha*eta.rho + b$w1*alpha.rho
     T_rho <- w.rho / b$abs.w
     z2 <- b$dlink.mu*y.mu    # term for the derivative of E
     w1.rho <- matrix(0,n,nsp)  # define derivatives of the diagonal elements of W1
     T1_rho <- matrix(0,n,nsp)  # define diagonal elements of T1_j
     Q_rho <- matrix(0,q,nsp)   # derivative of E diagonal wrt log(sp[j]) 
     w1.rho <- -a2*eta.rho
     T1_rho <- w1.rho/b$w1
     term <- T1_rho*z2 + a1*eta.rho- eta.rho
     Q_rho <- N_rho*drop(b$C2diag*crossprod(b$X,b$w1*z2)) +  
           b$C1diag*crossprod(b$X,b$w1*term)    
     ## efficient version of derivative of trA...
     KtIL <- t((b$L*b$I.plus)*b$K)
     KtILK <- KtIL%*%b$K
     KKtILK <- b$K%*%KtILK
     trA.rho<-rep(0,nsp)
     for (j in 1:nsp){
           trA.rho[j] <- - 2*sum(KtILK*(b$KtIQ1R%*%(N_rho[,j]*b$P))) -
                  sum((T_rho[,j]*KKtILK)*b$K) -  
                  sp[j]*sum((t(b$P)%*%b$S[[j]]%*%b$P)*t(KtILK) )+
                  sum( (t(b$P)%*%(c(Q_rho[,j])*b$P))*t(KtILK) ) +
                  2*sum( b$KtILQ1R*t(N_rho[,j]*b$P) ) +
                  sum(KtIL*t(T1_rho[,j]*b$K))
        }
    list(trA.rho=trA.rho,D.rho=D.rho)
  } ## deriv.dev.edf()

  crit.gcv <- function(fit,G,n){
  ## function to get the GCV criterion value
  ## gcv=dev*nobs/(nobs-gamma*trA)^2
    fit$gcv
  } 

  crit.ubre <- function(fit,G,n){
  ## function to get the UBRE criterion value
  ## ubre = dev/n - sig2 +2*gamma*trA*sig2/n
    ## cr <- fit$deviance +2*G$gamma*fit$trA*G$sig2
    cr <- fit$deviance/n +2*G$gamma*fit$trA*G$sig2/n - G$sig2
    cr
  } 

  if (G$scale.known) ## use UBRE...
         crit <- crit.ubre
  else crit <- crit.gcv ## use GCV
 
  nsp <- length(G$S) # number of smoothing parameters
  spind <- 1:nsp ## index of smoothing params in lsp
  lsp[spind] <- lsp[spind] + 2.5 
  mult <- 1
  fit <- scam.fit(G=G,sp=exp(lsp), env=env,control=control)
  q <- ncol(G$X)
  n <- nrow(G$X)
  score.hist <- rep(0,200)
   
  D.rho <- tau.rho <- rep(0,nsp) ## initialize derivatives of deviance and edf (tau/trA) wrt log(sp)
  for (iter in 1:200) { ## loop for GCV/UBRE minimization...
     ## calculation of the derivatives of deviance and edf wrt rho...
     der <-  deriv.dev.edf(G=G, fit=fit,sp=exp(lsp))
     D.rho <- der$D.rho
     tau.rho <- der$trA.rho
     
     if (G$scale.known)
        a <- pmax(0,-2*G$gamma*G$sig2*tau.rho)
     else a <- pmax(0,-2*G$gamma*fit$deviance*tau.rho/(n-fit$trA))
     r <- a/pmax(0,D.rho)
    # r[a==0 & D.rho==0] <- 1 
     r[a==0 | D.rho==0] <- 1 
     r[!is.finite(r)] <- 1e6
     lsp1 <- lsp
     lsp1 <- pmin(lsp + log(r)*mult,control$efs.lspmax)
     max.step <- max(abs(lsp1-lsp))
     old.gcv <- crit(fit, G,n) ## fit$gcv
     fit <- scam.fit(G=G,sp=exp(lsp1), env=env,control=control)
      ## some step length control...
     fit.gcv <- crit(fit, G,n) ## undated gcv/ubre value
   
     if (fit.gcv<=old.gcv) { ## improvement
       if (max.step<.05) { ## consider step extension (near optimum)
         lsp2 <- lsp
         lsp2 <- pmin(lsp + log(r)*mult*2,control$efs.lspmax) ## try extending step...
         fit2 <- scam.fit(G=G,sp=exp(lsp2), env=env,control=control)
         fit2.gcv <- crit(fit2,G,n)
         if (fit2.gcv < fit.gcv) { ## improvement - accept extension
           fit <- fit2;lsp <- lsp2
           fit.gcv <-  fit2.gcv 
	   mult <- mult * 2
          } else { ## accept old step
            lsp <- lsp1
           }
         } else lsp <- lsp1
     } else { ## no improvement 
      # while (fit$gcv > old.gcv&&mult>1) { ## don't contract below 1 as update doesn't have to improve GCV
      #     mult <- mult/2 ## contract step
      #	   lsp1 <- lsp
      #     lsp1 <- pmin(lsp + log(r)*mult,control$efs.lspmax)
      #	   fit <- scam.fit(G=G,sp=exp(lsp1), gamma=gamma, env=env,control=control)
      #  }
        gcv.thresh <- 10*(.1 +abs(old.gcv))*.Machine$double.eps^.5
        ii <- 1
        maxHalf.fit <- 15
        while (is.na(fit.gcv) || (fit.gcv-old.gcv) > gcv.thresh) { # 'step reduction' approach
             if (ii > maxHalf.fit) 
                 break ## stop ("step reduction failed")
             ii <- ii+1
             mult <- mult/2 ## contract step
             lsp1 <- lsp
             lsp1 <- pmin(lsp + log(r)*mult,control$efs.lspmax)
             fit <- scam.fit(G=G,sp=exp(lsp1), env=env,control=control)
             fit.gcv <- crit(fit,G,n)
         }
        lsp <- lsp1
        if (mult<1) mult <- 1
      }
     score.hist[iter] <- fit.gcv
     ## break if EFS step small and GCV change negligible over last 3 steps.
     if (iter>3 && max.step<.05 && max(abs(diff(score.hist[(iter-3):iter])))<control$efs.tol) break
     ## or break if deviance not changing...
     if (iter==1) old.dev <- fit$deviance else {
       if (abs(old.dev-fit$deviance) < 100*control$devtol.fit*abs(fit$deviance)) break
       old.dev <- fit$deviance
     }
  } # end of 'for' loop for GCV/UBRE minimization
 
  fit$sp <- exp(lsp)
  fit$niter <- iter
  fit$outer.info <- list(iter = iter,score.hist=score.hist[1:iter])
  fit$outer.info$conv <- if (iter==200) "iteration limit reached" else "full convergence"
  fit$gcv.ubre <- fit.gcv
  fit$inner.info <- list(iter=fit$iter,conv=fit$conv,pdev.hist=fit$pdev.hist) ## convergence info from the Newton optimization
  fit
} ## efsudr.scam



efsudr.scam1 <- function(G,lsp,env, control)
## Extended Fellner-Schall method for regular families,
## with PIRLS/full-Newton method used for the coefficient estimation, performed by scam.fit().
## Here only I use the first and the second order derivatives of the log likelihood when calculating
## the derivatives of the tau/edf. The dependence of H on lambda is neglected,
## as it is in Wood and Fasiolo (2017)  

{  deriv.dev.edf <- function(G,fit,sp){
  ## function to calculate derivatives of deviance and tr(A) w.r.t. log(sp),
  ## using expected Hessian and neglecting its dependence on sp
     b <- fit
     nsp <- length(G$S) 
     H <- crossprod(b$wX1) ## expected Hessian = CtXtWXC 
     Vb <- tcrossprod(b$P) ## P%*%t(P)*sig2 Bayesian posterior covariance matrix for the parameters
                 ## but without ? scale parameter here (inverse of the Hessian)
     trA.rho<- rep(0,nsp)
     beta.rho <- matrix(0,ncol(G$X),nsp) ## matrix of the parameters derivatives wrt log(sp)
    ## VbS <- matrix(0,dim(Vb)[1],dim(Vb)[2])
     for (j in 1:nsp){
           VbS <- Vb%*%b$S[[j]] 
           trA.rho[j] <- -sp[j]*sum((VbS%*%Vb)*t(H)) ## -lambda*tr(Vb*Sj*Vb*H)
           beta.rho[,j] <- -sp[j]*(VbS%*%b$beta) 
        }
     y.mu <- drop(b$y)-b$mu
     c <- -2*y.mu/(b$Var*b$dlink.mu)
     D.beta <- t(b$X1)%*%c  ## derivative of the deviance w.r.t. beta   
     D.rho <- t(D.beta)%*%beta.rho  ## derivative of the deviance w.r.t. log(sp)
     list(trA.rho=trA.rho,D.rho=D.rho)
  } ## deriv.dev.edf()

  crit.gcv <- function(fit,G,n){
  ## function to get the GCV criterion value
  ## gcv=dev*nobs/(nobs-gamma*trA)^2
    fit$gcv
  } 

  crit.ubre <- function(fit,G,n){
  ## function to get the UBRE criterion value
  ## ubre = dev/n - sig2 +2*gamma*trA*sig2/n
    ## cr <- fit$deviance +2*G$gamma*fit$trA*G$sig2
    cr <- fit$deviance/n +2*G$gamma*fit$trA*G$sig2/n - G$sig2
    cr
  } 

  if (G$scale.known) ## use UBRE...
         crit <- crit.ubre
  else crit <- crit.gcv ## use GCV
 
  nsp <- length(G$S) # number of smoothing parameters
  spind <- 1:nsp ## index of smoothing params in lsp
  lsp[spind] <- lsp[spind] + 2.5 
  mult <- 1
  fit <- scam.fit(G=G,sp=exp(lsp), env=env,control=control)
  q <- ncol(G$X)
  n <- nrow(G$X)
  score.hist <- rep(0,200)
   
  D.rho <- tau.rho <- rep(0,nsp) ## initialize derivatives of deviance and edf (tau/trA) wrt log(sp)
  for (iter in 1:200) { ## loop for GCV/UBRE minimization...
     ## calculation of the derivatives of deviance and edf wrt rho...
     der <-  deriv.dev.edf(G=G, fit=fit,sp=exp(lsp))
     D.rho <- der$D.rho
     tau.rho <- der$trA.rho
     
     if (G$scale.known)
        a <- pmax(0,-2*G$gamma*G$sig2*tau.rho)
     else a <- pmax(0,-2*G$gamma*fit$deviance*tau.rho/(n-fit$trA))
     r <- a/pmax(0,D.rho)
    # r[a==0 & D.rho==0] <- 1 
     r[a==0 | D.rho==0] <- 1 
     r[!is.finite(r)] <- 1e6
     lsp1 <- lsp
     lsp1 <- pmin(lsp + log(r)*mult,control$efs.lspmax)
     max.step <- max(abs(lsp1-lsp))
     old.gcv <- crit(fit, G,n) ## fit$gcv
     fit <- scam.fit(G=G,sp=exp(lsp1), env=env,control=control)
      ## some step length control...
     fit.gcv <- crit(fit, G,n) ## undated gcv/ubre value
   
     if (fit.gcv<=old.gcv) { ## improvement
       if (max.step<.05) { ## consider step extension (near optimum)
         lsp2 <- lsp
         lsp2 <- pmin(lsp + log(r)*mult*2,control$efs.lspmax) ## try extending step...
         fit2 <- scam.fit(G=G,sp=exp(lsp2), env=env,control=control)
         fit2.gcv <- crit(fit2,G,n)
         if (fit2.gcv < fit.gcv) { ## improvement - accept extension
           fit <- fit2;lsp <- lsp2
           fit.gcv <-  fit2.gcv 
	   mult <- mult * 2
          } else { ## accept old step
            lsp <- lsp1
           }
         } else lsp <- lsp1
     } else { ## no improvement 
      # while (fit$gcv > old.gcv&&mult>1) { ## don't contract below 1 as update doesn't have to improve GCV
      #     mult <- mult/2 ## contract step
      #	   lsp1 <- lsp
      #     lsp1 <- pmin(lsp + log(r)*mult,control$efs.lspmax)
      #	   fit <- scam.fit(G=G,sp=exp(lsp1), gamma=gamma, env=env,control=control)
      #  }
        gcv.thresh <- 10*(.1 +abs(old.gcv))*.Machine$double.eps^.5
        ii <- 1
        maxHalf.fit <- 15
        while (is.na(fit.gcv) || (fit.gcv-old.gcv) > gcv.thresh) { # 'step reduction' approach
             if (ii > maxHalf.fit) 
                 break ## stop ("step reduction failed")
             ii <- ii+1
             mult <- mult/2 ## contract step
             lsp1 <- lsp
             lsp1 <- pmin(lsp + log(r)*mult,control$efs.lspmax)
             fit <- scam.fit(G=G,sp=exp(lsp1), env=env,control=control)
             fit.gcv <- crit(fit,G,n)
         }
        lsp <- lsp1
        if (mult<1) mult <- 1
      }
     score.hist[iter] <- fit.gcv
     ## break if EFS step small and GCV change negligible over last 3 steps.
     if (iter>3 && max.step<.05 && max(abs(diff(score.hist[(iter-3):iter])))<control$efs.tol) break
     ## or break if deviance not changing...
     if (iter==1) old.dev <- fit$deviance else {
       if (abs(old.dev-fit$deviance) < 100*control$devtol.fit*abs(fit$deviance)) break
       old.dev <- fit$deviance
     }
  } # end of 'for' loop for GCV/UBRE minimization
 
  fit$sp <- exp(lsp)
  fit$niter <- iter
  fit$outer.info <- list(iter = iter,score.hist=score.hist[1:iter])
  fit$outer.info$conv <- if (iter==200) "iteration limit reached" else "full convergence"
  fit$gcv.ubre <- fit.gcv
  fit$inner.info <- list(iter=fit$iter,conv=fit$conv,pdev.hist=fit$pdev.hist) ## convergence info from the Newton optimization
  fit
} ## efsudr.scam1








