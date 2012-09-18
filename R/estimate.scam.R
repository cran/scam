
#########################################################
# Function to return gcv/ubre ...                      ##
#########################################################

gcv.ubre <- function(rho,G,gamma,ee,eb,esp,SVD){
   ## function to get GCV.UBRE value for optim()...
    if (length(rho)!= length(G$off)) stop (paste("length of rho and n.terms has to be the same"))
    sp <- exp(rho)
    b <- scam.fit(G=G, sp=sp, SVD=SVD,ee=ee,eb=eb,esp=esp)

    if (G$scale.known){ #  value of Mallow's Cp/UBRE/AIC ....
           n <- nrow(G$X)
           gcv.ubre <- b$dev/n - G$sig2 +2*gamma*b$trA*G$sig2/n
    }
    else      # value of GCV ...
           gcv.ubre <- b$gcv
    return(gcv.ubre)
}

#########################################################
## function to get the gradient of the gcv/ubre.....   ##
#########################################################

gcv.ubre.derivative <- function(rho,G, gamma,ee,eb,esp,SVD,
                            check.analytical=FALSE, del)
{# gradient of the gcv function 
 # G - object from scam setup 
 # rho - values of the logarithms of the smoothing parameters
 # X - model matrix
 # S - total penalty matrix
 # q0 - number of the parameters of the strictly parametric model
 # p.ident - vector of 0's & 1's for the model parameters identification: 
 # 1's stand for the parameters which will be exponentiated, 0's - otherwise
 # n.terms - number of the smooth terms in the mono-GAM
 # gamma - an ad hoc parametrer of the GCV
 # check.analytical - logical whether the analytical gradient of GCV/UBRE should be checked
 # del - increment for finite differences when checking analytical gradients
 ### ------------------------------------------------------------------
  y <- G$y;  X <- G$X;  S <- G$S
  q0 <- G$q0; q.f <- G$q.f
  p.ident <- G$p.ident; n.terms <- G$n.terms
  family <- G$family; intercept <- G$intercept; offset <- G$offset;
  weights <- G$weights;  scale.known <- G$scale.known; sig2 <- G$sig2
  n.pen <- length(S) # number of penalties
  if (length(rho)!=n.pen) stop (paste("length of rho and penalties has to be the same"))
  sp <- exp(rho)
 # fit the model with the given values of the smoothing parameters ----------
  b <- scam.fit(G=G, sp=sp, SVD=SVD,ee=ee,eb=eb,esp=esp)

  n <- nrow(b$X)
  q <- ncol(b$X)
  c <- -2*(y-b$mu)/(b$Var*b$dlink.mu) 
  diag.C<-rep(1,q)   # diagonal elements of matrix C with estimates of beta at convergence
  diag.C[b$iv]<-b$beta.t[b$iv]
  D.beta<-(diag.C*t(b$X))%*%c  # derivative of the deviance w.r.t. beta
  PPt<-b$P%*%t(b$P)  # inverse matrix of the Hessian
# ----------------------------------------------------------------------------
# calculation of the derivatives of beta by the Implicit Function Theorem starts here
# together with the derivative of the deviance wrt rho
  beta.rho<-matrix(0,q,n.pen) # define matrix of the parameters derivatives
  D.rho<-rep(0,n.pen)         # define derivative of the deviance wrt rho
  for (j in 1:n.pen)
     {beta.rho[,j]<--sp[j]*PPt%*%(b$S[[j]]%*%b$beta) # derivative of beta wrt rho[j]
      D.rho[j]<-t(D.beta)%*%beta.rho[,j]  # derivative of the deviance wrt rho[j] 
      }
# end of calculating the parameters derivatives and derivatives of  the deviance 
# ----------------------------------------------------------------------------
# calculating the derivative of tr(A) w.r.t. log(sp)
  wX1<-b$wX1                     # weighted model matrix
  a1<-(y-b$mu)*b$d2link.mu/b$dlink.mu    # a constant for updating the derivative of z
  a2<-(b$w)^2*(b$dvar.mu*b$dlink.mu+2*b$Var*b$d2link.mu) # a constant for updating the derivative of w
  eta.rho<-matrix(0,n,n.pen) # define derivatives of the linear predictor
  N_rho<-matrix(0,q,n.pen) # define diagonal elements of N_j
  w.rho<-matrix(0,n,n.pen) # define derivatives of the diagonal elements of W
  alpha.rho<-matrix(0,n,n.pen) # define derivatives of alpha
  T_rho<-matrix(0,n,n.pen) # define diagonal elements of T_j
     # a constant for the derivative of the alpha
  alpha1<-{-1/b$dlink.mu*(b$dvar.mu/b$Var+b$d2link.mu/b$dlink.mu)-
       1/b$dlink.mu*(y-b$mu)*((b$dvar.mu/b$Var)^2+
      (b$d2link.mu/b$dlink.mu)^2-b$d2var.mu/b$Var-b$d3link.mu/b$dlink.mu)}
  for (j in 1:n.pen)
    {eta.rho[,j]<-b$X1%*%beta.rho[,j]
     N_rho[b$iv,j]<-beta.rho[,j][b$iv]
     alpha.rho[,j]<-eta.rho[,j]*alpha1
     w.rho[,j]<--a2*b$alpha*eta.rho[,j]+b$w1*alpha.rho[,j]
     T_rho[,j]<-c(w.rho[,j]/b$abs.w)  
      }
  X2<-t((c(b$w1)*b$X)%*%b$C1) # term for the derivative of E
  z2<-b$dlink.mu*(y-b$mu)     # term for the derivative of E
  w1.rho<-matrix(0,n,n.pen) # define derivatives of the diagonal elements of W1
  T1_rho<-matrix(0,n,n.pen) # define diagonal elements of T1_j
# vector Q_rho[j] - derivative of E diagonal wrt log(sp[j])  
  Q_rho<-matrix(0,q,n.pen)
  for (j in 1:n.pen)
    {w1.rho[,j]<--a2*eta.rho[,j]
     T1_rho[,j]<-c(w1.rho[,j]/b$w1) 
     Q_rho[,j]<-{(N_rho[,j]*X2)%*%z2+X2%*%(T1_rho[,j]*z2)+
             X2%*%(a1*eta.rho[,j])-X2%*%eta.rho[,j]}}
#  KPt<-b$K%*%t(b$P)
#  KKt<-b$K%*%t(b$K)
   KtILK<-t(b$K)%*%(b$L*b$I.plus*b$K)
    # derivative of the trace of A
  trA.rho<-rep(0,n.pen)
  for (j in 1:n.pen)
    {# trA.rho[j]<-{-sum(diag((b$L*b$I.plus*KPt)%*%(N_rho[,j]*t(wX1))%*%(b$I.plus*KKt)))-
   # sum(diag((b$L*b$I.plus*KKt)%*%(T_rho[,j]*KKt)))-
   # sum(diag((b$L*b$I.plus*KKt)%*%(b$I.plus*wX1)%*%(N_rho[,j]*t(KPt))))-
  #  sp[j]*sum(diag((b$L*b$I.plus*KPt)%*%b$S[[j]]%*%t(KPt)))+
  #  sum(diag((b$L*b$I.plus*KPt)%*%(c(Q_rho[,j])*t(KPt))))+
  #  sum(diag((b$L*b$I.plus*KPt)%*%(N_rho[,j]*t(wX1))))+
  #  sum(diag(b$L*b$I.plus*T1_rho[,j]*KKt))+
  #  sum(diag((N_rho[,j]*t(KPt))%*%(b$L*b$I.plus*wX1)))}
  trA.rho[j]<-{-sum((t(b$P)%*%(N_rho[,j]*t(wX1)))*t((b$I.plus*b$K)%*%KtILK))-
    sum((t(b$K)%*%(T_rho[,j]*b$K))*t(KtILK))-
    sum((t(b$K)%*%(b$I.plus*wX1))*t((N_rho[,j]*b$P)%*%KtILK))-   
    sp[j]*sum((t(b$P)%*%b$S[[j]]%*%b$P)*t(KtILK))+
    sum((t(b$P)%*%(c(Q_rho[,j])*b$P))*t(KtILK))+
    sum(((N_rho[,j]*t(wX1))%*%((b$L*b$I.plus*b$K)))*(b$P))+
    sum((b$L*b$I.plus*T1_rho[,j]*b$K)*b$K)+
    sum((N_rho[,j]*b$P)*t(t(b$K)%*%(b$L*b$I.plus*wX1)))}
      } 
# Calculating the derivatives of the trA is completed here -----------------
# --------------------------------------------------------------------------
 if (scale.known){ #  derivative of Mallow's Cp/UBRE/AIC wrt log(sp) ....
       ubre.rho <- rep(0,n.pen) 
       for (j in 1:n.pen)
          {ubre.rho[j] <- D.rho[j]/n +2*gamma*trA.rho[j]*sig2/n
       }  
      gcv.ubre.rho <- ubre.rho
      gcv.ubre <- b$dev/n - sig2 +2*gamma*b$trA*sig2/n
  }
  else { # derivative of GCV wrt log(sp) ...
       gcv.rho<-rep(0,n.pen) # define the derivatives of the gcv
       for (j in 1:n.pen)
           {gcv.rho[j]<-n*(D.rho[j]*(n-gamma*b$trA)+2*gamma*b$dev*trA.rho[j])/(n-gamma*b$trA)^3
       }  
      gcv.ubre <- b$gcv
      gcv.ubre.rho <- gcv.rho
  }

# checking the derivatives by finite differences----------------------------
   dgcv.ubre.check <- NULL
   if (check.analytical){
      del <- 1e-4
      sp1 <- rep(0,n.pen)
      dbeta.check <- matrix(0,q,n.pen)
      dtrA.check <- rep(0,n.pen)
      dgcv.ubre.check <- rep(0,n.pen)
      for (j in 1:n.pen){
          sp1<-sp; sp1[j]<-sp[j]*exp(del)
          b1 <- scam.fit(G=G,sp=sp1,SVD,ee,eb,esp)
       ## calculating the derivatives of beta estimates by finite differences...
          dbeta.check[,j]<-1/del*(b1$beta-b$beta)
       ## calculating the derivatives of the trA by finite differences...
          dtrA.check[j]<-(1/del)*(b1$trA-b$trA)
       ## calculating derivatives of GCV/UBRE by finite differences...
          if (scale.known) gcv.ubre1 <- b1$dev/n - sig2 +2*gamma*b1$trA*sig2/n
          else gcv.ubre1 <- b1$gcv
          dgcv.ubre.check[j] <- (1/del)*(gcv.ubre1-gcv.ubre)
       }
  }
  check.grad <- NULL
  if (check.analytical){
          check.grad <- 100*(gcv.ubre.rho-dgcv.ubre.check)/dgcv.ubre.check
  }      

# end of checking the derivatives -------------------------------------------
  return(gcv.ubre.rho)
}



#############################################################################
## for nlm() function to get the gcv/ubre and gradient of the gcv/ubre.....##
#############################################################################

dgcv.ubre.nlm <- function(rho,G, gamma,ee,eb,esp,SVD,
                            check.analytical=FALSE, del){## gradient of the gcv/ubre function
 # G - object from scam setup 
 # rho - values of the logarithms of the smoothing parameters
 # X - model matrix
 # S - total penalty matrix
 # q0 - number of the parameters of the strictly parametric model
 # p.ident - vector of 0's & 1's for the model parameters identification: 
 # 1's stand for the parameters which will be exponentiated, 0's - otherwise
 # n.terms - number of the smooth terms in the mono-GAM
 # gamma - an ad hoc parametrer of the GCV
 # check.analytical - logical whether the analytical gradient of GCV/UBRE should be checked
 # del - increment for finite differences when checking analytical gradients
 ### ------------------------------------------------------------------
  y <- G$y;  X <- G$X;  S <- G$S
  q0 <- G$q0; q.f <- G$q.f
  p.ident <- G$p.ident; n.terms <- G$n.terms
  family <- G$family; intercept <- G$intercept; offset <- G$offset;
  weights <- G$weights;  scale.known <- G$scale.known; sig2 <- G$sig2
  n.pen <- length(S) # number of penalties
  if (length(rho)!=n.pen) stop (paste("length of rho and # penalties has to be the same"))
  sp<-exp(rho)
 # fit the model with the given values of the smoothing parameters ----------
  b <- scam.fit(G=G,sp=sp,SVD,ee,eb,esp)
                  
  n<-nrow(b$X)
  q<-ncol(b$X)
  c<--2*(y-b$mu)/(b$Var*b$dlink.mu) 
  diag.C<-rep(1,q)   # diagonal elements of matrix C with estimates of beta at convergence
  diag.C[b$iv]<-b$beta.t[b$iv]
  D.beta<-(diag.C*t(b$X))%*%c  # derivative of the deviance w.r.t. beta
  PPt<-b$P%*%t(b$P)  # inverse matrix of the Hessian
# ----------------------------------------------------------------------------
# calculation of the derivatives of beta by the Implicit Function Theorem starts here
# together with the derivative of the deviance wrt rho
  beta.rho<-matrix(0,q,n.pen) # define matrix of the parameters derivatives
  D.rho<-rep(0,n.pen)         # define derivative of the deviance wrt rho
  for (j in 1:n.pen)
     {beta.rho[,j]<--sp[j]*PPt%*%(b$S[[j]]%*%b$beta) # derivative of beta wrt rho[j]
      D.rho[j]<-t(D.beta)%*%beta.rho[,j]  # derivative of the deviance wrt rho[j] 
      }
# end of calculating the parameters derivatives and derivatives of  the deviance 
# ----------------------------------------------------------------------------
# calculating the derivative of tr(A) w.r.t. log(sp)
  wX1<-b$wX1                     # weighted model matrix
  a1<-(y-b$mu)*b$d2link.mu/b$dlink.mu    # a constant for updating the derivative of z
  a2<-(b$w)^2*(b$dvar.mu*b$dlink.mu+2*b$Var*b$d2link.mu) # a constant for updating the derivative of w
  eta.rho<-matrix(0,n,n.pen) # define derivatives of the linear predictor
  N_rho<-matrix(0,q,n.pen) # define diagonal elements of N_j
  w.rho<-matrix(0,n,n.pen) # define derivatives of the diagonal elements of W
  alpha.rho<-matrix(0,n,n.pen) # define derivatives of alpha
  T_rho<-matrix(0,n,n.pen) # define diagonal elements of T_j
     # a constant for the derivative of the alpha
  alpha1<-{-1/b$dlink.mu*(b$dvar.mu/b$Var+b$d2link.mu/b$dlink.mu)-
       1/b$dlink.mu*(y-b$mu)*((b$dvar.mu/b$Var)^2+
      (b$d2link.mu/b$dlink.mu)^2-b$d2var.mu/b$Var-b$d3link.mu/b$dlink.mu)}
  for (j in 1:n.pen)
    {eta.rho[,j]<-b$X1%*%beta.rho[,j]
     N_rho[b$iv,j]<-beta.rho[,j][b$iv]
     alpha.rho[,j]<-eta.rho[,j]*alpha1
     w.rho[,j]<--a2*b$alpha*eta.rho[,j]+b$w1*alpha.rho[,j]
     T_rho[,j]<-c(w.rho[,j]/b$abs.w)  
      }
  X2<-t((c(b$w1)*b$X)%*%b$C1) # term for the derivative of E
  z2<-b$dlink.mu*(y-b$mu)     # term for the derivative of E
  w1.rho<-matrix(0,n,n.pen) # define derivatives of the diagonal elements of W1
  T1_rho<-matrix(0,n,n.pen) # define diagonal elements of T1_j
# vector Q_rho[j] - derivative of E diagonal wrt log(sp[j])  
  Q_rho<-matrix(0,q,n.pen)
  for (j in 1:n.pen)
    {w1.rho[,j]<--a2*eta.rho[,j]
     T1_rho[,j]<-c(w1.rho[,j]/b$w1) 
     Q_rho[,j]<-{(N_rho[,j]*X2)%*%z2+X2%*%(T1_rho[,j]*z2)+
             X2%*%(a1*eta.rho[,j])-X2%*%eta.rho[,j]}}
#  KPt<-b$K%*%t(b$P)
#  KKt<-b$K%*%t(b$K)
   KtILK<-t(b$K)%*%(b$L*b$I.plus*b$K)
    # derivative of the trace of A
  trA.rho<-rep(0,n.pen)
  for (j in 1:n.pen)
    {# trA.rho[j]<-{-sum(diag((b$L*b$I.plus*KPt)%*%(N_rho[,j]*t(wX1))%*%(b$I.plus*KKt)))-
   # sum(diag((b$L*b$I.plus*KKt)%*%(T_rho[,j]*KKt)))-
   # sum(diag((b$L*b$I.plus*KKt)%*%(b$I.plus*wX1)%*%(N_rho[,j]*t(KPt))))-
  #  sp[j]*sum(diag((b$L*b$I.plus*KPt)%*%b$S[[j]]%*%t(KPt)))+
  #  sum(diag((b$L*b$I.plus*KPt)%*%(c(Q_rho[,j])*t(KPt))))+
  #  sum(diag((b$L*b$I.plus*KPt)%*%(N_rho[,j]*t(wX1))))+
  #  sum(diag(b$L*b$I.plus*T1_rho[,j]*KKt))+
  #  sum(diag((N_rho[,j]*t(KPt))%*%(b$L*b$I.plus*wX1)))}
  trA.rho[j]<-{-sum((t(b$P)%*%(N_rho[,j]*t(wX1)))*t((b$I.plus*b$K)%*%KtILK))-
    sum((t(b$K)%*%(T_rho[,j]*b$K))*t(KtILK))-
    sum((t(b$K)%*%(b$I.plus*wX1))*t((N_rho[,j]*b$P)%*%KtILK))-   
    sp[j]*sum((t(b$P)%*%b$S[[j]]%*%b$P)*t(KtILK))+
    sum((t(b$P)%*%(c(Q_rho[,j])*b$P))*t(KtILK))+
    sum(((N_rho[,j]*t(wX1))%*%((b$L*b$I.plus*b$K)))*(b$P))+
    sum((b$L*b$I.plus*T1_rho[,j]*b$K)*b$K)+
    sum((N_rho[,j]*b$P)*t(t(b$K)%*%(b$L*b$I.plus*wX1)))}
      } 
# Calculating the derivatives of the trA is completed here -----------------
# --------------------------------------------------------------------------
 if (scale.known){ #  derivative of Mallow's Cp/UBRE/AIC wrt log(sp) ....
       ubre.rho <- rep(0,n.pen) 
       for (j in 1:n.pen)
          {ubre.rho[j] <- D.rho[j]/n +2*gamma*trA.rho[j]*sig2/n
       }  
      gcv.ubre.rho <- ubre.rho
      gcv.ubre <- b$dev/n - sig2 +2*gamma*b$trA*sig2/n
  }
  else { # derivative of GCV wrt log(sp) ...
       gcv.rho<-rep(0,n.pen) # define the derivatives of the gcv
       for (j in 1:n.pen)
           {gcv.rho[j]<-n*(D.rho[j]*(n-gamma*b$trA)+2*gamma*b$dev*trA.rho[j])/(n-gamma*b$trA)^3
       }  
      gcv.ubre <- b$gcv
      gcv.ubre.rho <- gcv.rho
  }
  # checking the derivatives by finite differences----------------------------
   dgcv.ubre.check <- NULL
   if (check.analytical){
      del <- 1e-4
      sp1 <- rep(0,n.pen)
      dbeta.check <- matrix(0,q,n.pen)
      dtrA.check <- rep(0,n.pen)
      dgcv.ubre.check <- rep(0,n.pen)
      for (j in 1:n.pen){
          sp1<-sp; sp1[j]<-sp[j]*exp(del)
          b1 <- scam.fit(G=G,sp=sp1,SVD,ee,eb,esp)
       ## calculating the derivatives of beta estimates by finite differences...
          dbeta.check[,j]<-1/del*(b1$beta-b$beta)
       ## calculating the derivatives of the trA by finite differences...
          dtrA.check[j]<-(1/del)*(b1$trA-b$trA)
       ## calculating derivatives of GCV/UBRE by finite differences...
          if (scale.known) gcv.ubre1 <- b1$dev/n - sig2 +2*gamma*b1$trA*sig2/n
          else gcv.ubre1 <- b1$gcv
          dgcv.ubre.check[j] <- (1/del)*(gcv.ubre1-gcv.ubre)
       }
  }
  check.grad <- NULL
  if (check.analytical){
          check.grad <- 100*(gcv.ubre.rho-dgcv.ubre.check)/dgcv.ubre.check
  }      

# end of checking the derivatives ------------------------------
# -------------------------------------------------------
  attr(gcv.ubre,"gradient") <- gcv.ubre.rho
  gcv.ubre
}


#######################################################
#### estimate.scam()....                             ##
#######################################################


estimate.scam <- function(G,optimizer,optim.method,rho, gamma=1,
                    ee,eb,esp,SVD=TRUE, 
                    check.analytical, del,epsilon=epsilon){
## function to select smoothing parameter...
   if (!(optimizer %in% c("bfgs", "nlm", "optim","nlm.fd")) )
        stop("unknown outer optimization method.")
    
   if (optimizer == "bfgs"){ ## minimize GCV/UBRE by BFGS...
        b <- bfgs_gcv.ubre(gcv.ubre_grad,rho=rho, G=G,gamma=gamma,ee=ee,eb=eb,esp=esp,
                         check.analytical=check.analytical, del=del)
        sp <- exp(b$rho)
        object <- b$object
        object$gcv.ubre <- b$gcv.ubre
        object$dgcv.ubre <- b$dgcv.ubre
        object$termcode <- b$termcode
        object$check.grad <- b$check.grad
        object$dgcv.ubre.check <- b$dgcv.ubre.check
        object$conv.bfgs <- b$conv.bfgs
        object$iterations <- b$iterations
   }
   else if (optimizer=="optim"){  ## gr=gcv.ubre.derivative
          if (!(optim.method[1] %in% c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"))) {
              warning("unknown optim() method, `L-BFGS-B' were used")
              optim.method[1] <- "L-BFGS-B"
          }
          if (is.na(optim.method[2])) {
               warning("the second parameter of optim.method argument is not supplied, 
                   finite-difference approximation of the gradient were used")
               grr <- NULL
          }
          else if (!(optim.method[2] %in% c("fd","grad"))) {
               warning("only `fd' and `grad' options are possible, finite-difference 
                        approximation of the gradient were used")
               grr <- NULL
          }        
          else if (optim.method[2] == "grad")
                grr <- gcv.ubre.derivative
          else 
                grr <- NULL

         b <- optim(par=rho,fn=gcv.ubre, gr=grr, method=optim.method[1], 
                   G=G, gamma=gamma,ee=ee,eb=eb,esp=esp,SVD=TRUE) ##,
               ##     check.analytical=check.analytical, del=del)
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
   else if (optimizer=="nlm.fd"){## nlm() with finite difference derivatives...
          b <- nlm(f=gcv.ubre, p=rho,iterlim=100, G=G, gamma=gamma,ee=ee,eb=eb,esp=esp,
                    SVD=TRUE)
   }
   else if (optimizer=="nlm"){## nlm() with analytical derivatives...
          b <- nlm(f=dgcv.ubre.nlm, p=rho,iterlim=100,G=G,gamma=1,
                    ee=ee,eb=eb,esp=esp, SVD=TRUE,
                    check.analytical=check.analytical, del=del)
   }
   if (optimizer== "nlm.fd" || optimizer== "nlm") {
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
        object <- scam.fit(G=G, sp=sp, SVD=TRUE,
                      ee=ee,eb=eb,esp=esp,  epsilon=epsilon)
        object$gcv.ubre <- gcv.ubre
        object$dgcv.ubre <- dgcv.ubre 
        object$termcode <- termcode
        object$conv <- conv
        object$iterations <- iterations 
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
}


