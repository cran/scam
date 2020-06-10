
#########################################################
# Function to return gcv/ubre ...                      ##
#########################################################

gcv.ubre <- function(rho,G,gamma,env,control) 
{  ## function to get GCV.UBRE value for optim()...
   if (length(rho)!= length(G$off)) stop (paste("length of rho and n.terms has to be the same"))
   sp <- exp(rho)
   b <- scam.fit(G=G, sp=sp, env=env,gamma=gamma, control=scam.control()) 
   if (G$scale.known) #  value of Mallow's Cp/UBRE/AIC ....
      {  n <- nrow(G$X)
         gcv.ubre <- b$dev/n - G$sig2 +2*gamma*b$trA*G$sig2/n
      }  else   # value of GCV ...
           gcv.ubre <- b$gcv
   return(gcv.ubre)
}

#########################################################
## function to get the gradient of the gcv/ubre.....   ##
#########################################################

gcv.ubre.derivative <- function(rho,G, gamma,env, control)  
{  ## function to return derivative of GCV or UBRE for optim...
   gcv.ubre_grad(rho, G, gamma,env, control=control)$gcv.ubre.rho
}


#############################################################################
## for nlm() function to get the gcv/ubre and gradient of the gcv/ubre.....##
#############################################################################

dgcv.ubre.nlm <- function(rho,G, gamma,env,  control) 
{  ## GCV UBRE objective function for nlm
   gg <- gcv.ubre_grad(rho, G, gamma,env, control=control) 
   attr(gg$gcv.ubre,"gradient") <- gg$gcv.ubre.rho
   gg$gcv.ubre
}



#######################################################
#### estimate.scam()....                             ##
#######################################################


estimate.scam <- function(G,optimizer,optim.method,rho, gamma,env,  control)
                   ##  check.analytical, del, devtol.fit, steptol.fit)
{  ## function to select smoothing parameter...
   if (!(optimizer %in% c("bfgs", "nlm", "optim","nlm.fd","efs")) )
          stop("unknown outer optimization method")

   if (length(rho)==0) { ## no sp estimation to do -- run a fit instead
     optimizer <- "no.sps" ## will cause scam.fit to be called, below
   }

   if (optimizer == "bfgs") {## minimize GCV/UBRE by BFGS...
         b <- bfgs_gcv.ubre(gcv.ubre_grad,rho=rho, G=G,gamma=gamma,env=env, control=control) ## check.analytical=check.analytical, del=del,devtol.fit=devtol.fit, steptol.fit=steptol.fit) 
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
   else if (optimizer=="optim"){  ## gr=gcv.ubre.derivative
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
         list(factr=control$optim$factr,lmm=min(5,length(rho))), gamma=gamma,env=env) 
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
   else if (optimizer=="nlm.fd") {## nlm() with finite difference derivatives...
             b <- nlm(f=gcv.ubre, p=rho,typsize=rho, stepmax = control$nlm$stepmax, ndigit = control$nlm$ndigit,
	               gradtol = control$nlm$gradtol, steptol = control$nlm$steptol, 
                       iterlim = control$nlm$iterlim,  G=G, gamma=gamma,env=env, control=control) 
   }
   else if (optimizer=="nlm"){  ## nlm() with analytical derivatives...
            b <- nlm(f=dgcv.ubre.nlm, p=rho,typsize=rho, stepmax = control$nlm$stepmax, ndigit = control$nlm$ndigit, gradtol = control$nlm$gradtol, steptol = control$nlm$steptol, iterlim = control$nlm$iterlim, G=G,gamma=gamma,env=env, control=control) 
   }
   else if (optimizer=="efs"){  ## Extended Fellner-Schall method
             b <- efsudr.scam(G=G,lsp=rho,gamma=gamma,env=env, control=control)
             sp <- b$sp
             object <- b
             object$iterations <- b$niter
             object$conv <- b$outer.info$conv
             object$score.hist <- b$outer.info$score.hist
             object$gcv.ubre <- b$gcv
  }


   if (optimizer== "nlm.fd" || optimizer== "nlm"){ 
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
   if (optimizer== "nlm.fd" || optimizer== "nlm" || optimizer== "optim"){
         object <- scam.fit(G=G, sp=sp,env=env,control=control) 
         object$gcv.ubre <- gcv.ubre
         object$dgcv.ubre <- dgcv.ubre 
         object$termcode <- termcode
         object$conv <- conv
         object$iterations <- iterations 
   }
   else if (optimizer=="no.sps"){
         object <- scam.fit(G=G, sp=exp(rho),env=env,control=control)      
         sp <- rep(0,0)    
   }
   if (optimizer=="optim"){
         object$optim.method <- rep(NA,2)
         object$optim.method[1] <- optim.method[1]
         if (!is.null(grr))
              object$optim.method[2] <- "grad"
   }
   
   object$sp <- sp
   object$q.f <- G$q.f
   object$p.ident <- G$p.ident
   object$S <- G$S
   object
} ## estimate.scam

##########################################################################
## Extended Fellner-Schall method for regular families as in gam(mgcv)
## with PIRLS performed by scam.fit. Rather than minimizing REML (as in gam(mgcv)), 
## the function minimizes GCV criterion...
##########################################################################


efsudr.scam <- function(G,lsp,gamma,env, control){##maxit=200, devtol.fit=1e-7, steptol.fit=1e-7,
                                            # gamma=1, start=NULL, etastart=NULL, mustart=NULL, env=env){
## Extended Fellner-Schall method for regular families,
## with PIRLS performed by scam.fit

  deriv.dev.edf <- function(G,fit,sp){
  ## function to calculate derivatives of deviance and tr(A) w.r.t. log(sp)
     b <- fit
     nsp <- length(G$S) 
     q <- ncol(G$X)
     n <- nrow(G$X) 
     y.mu <- b$y-b$mu
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
  } # deriv.dev.edf
 
  nsp <- length(G$S) # number of smoothing parameters
  spind <- 1:nsp ## index of smoothing params in lsp
  lsp[spind] <- lsp[spind] + 2.5 
  mult <- 1
  fit <- scam.fit(G=G,sp=exp(lsp), gamma=gamma, env=env,control=control)
  q <- ncol(G$X)
  n <- nrow(G$X)
  score.hist <- rep(0,200)


  D.rho <- tau.rho <- rep(0,nsp) ## initialize derivatives of deviance and edf (tau/trA) wrt log(sp)
  for (iter in 1:200) {
     ## calculation of the derivatives of deviance and edf wrt rho...
     der <-  deriv.dev.edf(G=G, fit=fit,sp=exp(lsp))
     D.rho <- der$D.rho
     tau.rho <- der$trA.rho
     a <- pmax(0,-2*fit$deviance*tau.rho/(n-fit$trA))
     r <- a/pmax(0,D.rho)
    # r[a==0 & D.rho==0] <- 1 
     r[a==0 | D.rho==0] <- 1 
     r[!is.finite(r)] <- 1e6
     lsp1 <- lsp
     lsp1 <- pmin(lsp + log(r)*mult,control$efs.lspmax)
     max.step <- max(abs(lsp1-lsp))
     old.gcv <- fit$gcv
     fit <- scam.fit(G=G,sp=exp(lsp1), gamma=gamma, env=env,control=control)
      ## some step length control...
   
     if (fit$gcv<=old.gcv) { ## improvement
       if (max.step<.05) { ## consider step extension (near optimum)
         lsp2 <- lsp
         lsp2 <- pmin(lsp + log(r)*mult*2,control$efs.lspmax) ## try extending step...
         fit2 <- scam.fit(G=G,sp=exp(lsp2), gamma=gamma, env=env,control=control)
         if (fit2$gcv < fit$gcv) { ## improvement - accept extension
           fit <- fit2;lsp <- lsp2
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
        while (is.na(fit$gcv) || (fit$gcv-old.gcv) > gcv.thresh) { # 'step reduction' approach
             if (ii > maxHalf.fit) 
                 break ## stop ("step reduction failed")
             ii <- ii+1
             mult <- mult/2 ## contract step
             lsp1 <- lsp
             lsp1 <- pmin(lsp + log(r)*mult,control$efs.lspmax)
             fit <- scam.fit(G=G,sp=exp(lsp1), gamma=gamma, env=env,control=control)
         }
        lsp <- lsp1
        if (mult<1) mult <- 1
      }
     score.hist[iter] <- fit$gcv
     ## break if EFS step small and GCV change negligible over last 3 steps.
     if (iter>3 && max.step<.05 && max(abs(diff(score.hist[(iter-3):iter])))<control$efs.tol) break
     ## or break if deviance not changing...
     if (iter==1) old.dev <- fit$deviance else {
       if (abs(old.dev-fit$deviance) < 1000*control$devtol.fit*abs(fit$deviance)) break
       old.dev <- fit$deviance
     }
  } # end of 'for' loop
  fit$sp <- exp(lsp)
  fit$niter <- iter
  fit$outer.info <- list(iter = iter,score.hist=score.hist[1:iter])
  fit$outer.info$conv <- if (iter==200) "iteration limit reached" else "full convergence"
  fit
} ## efsudr.scam











