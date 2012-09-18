

## BFGS for GCV minimization 
## with modifications on convergence control, scaling, etc...


## function to get the gcv value and the gradient of the gcv

gcv.ubre_grad <- function(rho, G, gamma,ee,eb,esp,SVD=TRUE,
                         check.analytical, del){
# gradient of the gcv function....
  # G - object from scam setup  
 # y - values of the response variable
 # rho - values of the logarithms of the smoothing parameters
 # X - model matrix
 # S - total penalty matrix
 # q0 - number of the parameters of the strictly parametric model
 # p.ident - vector of 0's & 1's for the model parameters identification: 
 # 1's stand for the parameters which will be exponentiated, 0's - otherwise
 # n.terms - number of the smooth terms in the SCAM
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
  b<- scam.fit(G=G,sp=sp,SVD=SVD,ee=ee,eb=eb,esp=esp)
                  
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
   KtILK<-t(b$K)%*%(b$L*b$I.plus*b$K)
    # derivative of the trace of A
  trA.rho<-rep(0,n.pen)
  for (j in 1:n.pen)
    { trA.rho[j]<-{-sum((t(b$P)%*%(N_rho[,j]*t(wX1)))*t((b$I.plus*b$K)%*%KtILK))-
    sum((t(b$K)%*%(T_rho[,j]*b$K))*t(KtILK))-
    sum((t(b$K)%*%(b$I.plus*wX1))*t((N_rho[,j]*b$P)%*%KtILK))-   
    sp[j]*sum((t(b$P)%*%b$S[[j]]%*%b$P)*t(KtILK))+
    sum((t(b$P)%*%(c(Q_rho[,j])*b$P))*t(KtILK))+
    sum(((N_rho[,j]*t(wX1))%*%((b$L*b$I.plus*b$K)))*(b$P))+
    sum((b$L*b$I.plus*T1_rho[,j]*b$K)*b$K)+
    sum((N_rho[,j]*b$P)*t(t(b$K)%*%(b$L*b$I.plus*wX1)))}
      } 
# Calculating the derivatives of the trA is completed here -----
# -------------------------------------------------------------
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


# checking the derivatives by finite differences----------------
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

# end of checking the derivatives -------------------------------
 list(dgcv.ubre=gcv.ubre.rho,gcv.ubre=gcv.ubre, scale.est=b$dev/(n-b$trA), 
       check.grad=check.grad, dgcv.ubre.check=dgcv.ubre.check, object = b)
}



#####################################################
## BFGS for gcv miminization...                    ##
#####################################################

bfgs_gcv.ubre <- function(fn=gcv.ubre_grad,rho,ini.fd=TRUE,G,gamma=1,
           ee,eb,esp,SVD=TRUE, n.pen=length(rho),
          typx=rep(1,n.pen), typf=1, steptol= 1e-7, ## 3.66*1e-11, 
          gradtol = 6.0554*1e-06, maxNstep = 5, maxHalf = 30, 
          check.analytical, del) 
{
## fn - GCV/UBRE Function which returs the GCV/UBRE value and its derivative wrt log(sp)
## rho - log of the initial values of the smoothing parameters
## ini.fd - if TRUE, a finite difference to the Hessian is used to find the initial inverse Hessian
## typx - vector whose component is a positive scalar specifying the typical magnitude of sp
## typf - a positive scalar estimating the magnitude of the gcv near the minimum
## gradtol - a scalar giving a tolerance at which the gradient is considered
##            to be close enougth to 0 to terminate the algorithm
## steptol - a positive scalar giving the tolerance at which the scaled distance between
##          two successive iterates is considered close enough to zero to terminate the algorithm 
## maxNstep - a positive scalar which gives the maximum allowable step length
## maxHalf - a positive scalar which gives the maximum number of step halving in "backtracking"
   
 #  n.pen <- length(rho) # smoothing parameter dimension
   Sp <- 1/typx     # diagonal of the scaling matrix 
   ## storage for solution track
   rho1 <- rho
   old.rho <- rho
   
   b <- fn(rho,G,gamma=1,ee,eb,esp,SVD=TRUE, check.analytical=FALSE,del)
   old.score <- score <- b$gcv.ubre
  score.plot<-rep(0,200) ## for plotting the gcv
  score.plot[1] <- score

   grad <-  b$dgcv.ubre
   ## The initial inverse Hessian ...
   if (ini.fd) {
      B <- matrix(0,n.pen,n.pen)
      for (j in 1:n.pen) {
          rho2 <- rho; rho2[j] <- rho[j] + 1e-6
          b2 <- fn(rho2,G,gamma=1,ee,eb,esp,SVD=TRUE,check.analytical=FALSE, del)
          B[,j] <- (b2$dgcv.ubre - grad)/1e-6   
         }
      B <- B + t(B)
      eh <- eigen(B)
      ind <- eh$values < 0
      eh$values[ind] <- -eh$values[ind] ## condition number??????
      B <- eh$vectors%*%(t(eh$vectors)/eh$values)
     } else B <- diag(n.pen)*100
  
   score.scale <- b$scale.est + score  ## ??? for UBRE
   unconv.ind <- abs(grad) > score.scale*gradtol # checking the gradient is within the tolerance
   if (!sum(unconv.ind)) ## if at least one is false
       unconv.ind <- unconv.ind | TRUE

   consecmax <- 0
   ## Quasi-Newton algorithm to minimize GCV...
   for (i in 1:200) {  
      
      ## compute a BFGS search direction ...
      Nstep <-0*grad  # initialize the quasi-Newton step
      Nstep[unconv.ind] <- -drop(B[unconv.ind, unconv.ind]%*%grad[unconv.ind])
    #  Nstep <- -B%*%grad
      Dp <- Sp*Nstep
      Newtlen <- (sum(Dp^2))^0.5  ## Euclidean norm to get the Newton step
      if (Newtlen > maxNstep){ ## reduce if ms is greater than the max allowable
            Nstep <- maxNstep*Nstep/Newtlen
            Newtlen <- maxNstep
      }

      maxtaken <- FALSE
      retcode <- 2
      initslope <- sum(grad * Nstep) # initial slope
      rellength <- max(abs(Nstep)/max(abs(rho),1/Sp)) # relative length of rho for the stopping criteria
      minalpha <- steptol/rellength
      c1 <- 1e-4 # constant for the sufficient decrease condition 
      alpha <- 1 ## initialize step length
      ii <- 0 ## initialize the number of "step halving"
      step <- alpha*Nstep 
        
      ## step length selection ...
      
      curv.condition <- TRUE

      repeat {
          rho1 <- rho + alpha*Nstep
          b <- fn(rho=rho1,G,gamma=1,ee,eb,esp,SVD=TRUE,check.analytical=FALSE, del)
          score1 <- b$gcv.ubre
          if (score1 <= score+c1*alpha*initslope) { 
             grad1 <- b$dgcv.ubre
             newslope <- sum(grad1 * Nstep)
             curv.condition <- TRUE
             if (newslope < 0.9*initslope){ # the curvature condition is not satisfied
                 if (alpha == 1 && Newtlen<maxNstep) {
                      maxalpha <- maxNstep/Newtlen
                      repeat {
                           old.alpha <- alpha
                           old.score1 <- score1
                           alpha <- min(2*alpha, maxalpha)
                           rho1 <- rho + alpha*Nstep
                           b <- fn(rho=rho1,G,gamma=1,ee,eb,esp,SVD=TRUE,
                                     check.analytical=FALSE, del)
                           score1 <- b$gcv.ubre
                           if (score1 <= score+c1*alpha*initslope){
                                grad1 <- b$dgcv.ubre
                                newslope <- sum(grad1*Nstep)
                           }
                           if (score1 > score+c1*alpha*initslope)
                               break
                           if (newslope >= 0.9*initslope)
                               break
                           if (alpha >= maxalpha)
                               break
                      } 
                 }
                 if ((alpha < 1) || (alpha>1 && (score1>score+c1*alpha*initslope))){
                       alpha.lo <- min(alpha, old.alpha)
                       alpha.diff <- abs(old.alpha - alpha)
                       if (alpha < old.alpha){
                            sc.lo <- score1
                            sc.hi <- old.score1
                       }    
                       else{
                            sc.lo <- old.score1
                            sc.hi <- score1 
                       }
                       repeat{
                            alpha.incr <- -newslope*alpha.diff^2/(2*(sc.hi-(sc.lo+newslope*alpha.diff)))
                            if (alpha.incr < 0.2*alpha.diff) 
                                   alpha.incr <- 0.2*alpha.diff
                            alpha <- alpha.lo+alpha.incr
                            rho1 <- rho + alpha*Nstep
                            b <- fn(rho=rho1,G,gamma=1,ee,eb,esp,SVD=TRUE,
                                  check.analytical=FALSE, del)
                            score1 <- b$gcv.ubre
                            if (score1 > score+c1*alpha*initslope){
                                 alpha.diff <- alpha.incr
                                 sc.hi <- score1
                            }
                            else {
                                 grad1 <- b$dgcv.ubre
                                 newslope <- sum(grad1*Nstep)
                                 if (newslope < 0.9*initslope){
                                       alpha.lo <- alpha
                                       alpha.diff <- alpha.diff-alpha.incr
                                       sc.lo <- score1
                                 }
                            }
                            if (newslope >= 0.9*initslope)
                                 break
                            if (alpha.diff < minalpha)
                                 break 
                       }
                       if (newslope < 0.9*initslope){ ## couldn't satisfy curvature condition
                            curv.condition <- FALSE
                            score1 <- sc.lo
                            rho1 <- rho + alpha.lo*Nstep
                            b <- fn(rho=rho1,G,gamma=1,ee,eb,esp,SVD=TRUE,
                                   check.analytical=FALSE, del)
                       } 
                  }
               }
               retcode <- 0
               if (newslope < 0.9*initslope) ## couldn't satisfy curvature condition
                            curv.condition <- FALSE
               if (alpha*Newtlen > 0.99*maxNstep)  
                    maxtaken <- TRUE
          } ## end of (if (score1 <= ...))
          else if (alpha < minalpha){ ## no satisfactory rho+ can be found suff-ly distinct from previous rho
                   retcode <- 1
                   rho1 <- rho
                   b <- fn(rho=rho1,G,gamma=1,ee,eb,esp,SVD=TRUE,
                            check.analytical=FALSE,del)
               }
          else {  ## backtracking to satisfy the sufficient decrease condition...
             ii <- ii+1
             if (alpha ==1){ ## first backtrack, quadratic fit
                 alpha.temp <- -initslope/(2*(score1-score-initslope))
             }
             else { ## all subsequent backtracts, cubic fit 
                 A1 <- matrix(0,2,2)
                 bb1 <-rep(0,2)
                 ab <- rep(0,2)
                 A1[1,1] <- 1/alpha^2
                 A1[1,2] <- -1/old.alpha^2
                 A1[2,1] <- -old.alpha/alpha^2
                 A1[2,2] <- alpha/old.alpha^2
                 bb1[1] <- score1-score-alpha*initslope
                 bb1[2] <- old.score1 -score-old.alpha*initslope
                 ab <- 1/(alpha-old.alpha)*A1%*%bb1
                 disc <- ab[2]^2-3*ab[1]*initslope
                 if (ab[1] ==0) ## cubic is a quadratic
                     alpha.temp <- -initslope/(2*ab[2])
                 else     ## legitimate cubic
                     alpha.temp <- (-ab[2]+disc^0.5)/(3*ab[1])
                 if (alpha.temp > 0.5*alpha)
                      alpha.temp <- 0.5*alpha   
              }
              old.alpha <- alpha
              old.score1 <- score1
              if (alpha.temp <= 0.1*alpha) alpha <- 0.1*alpha
              else alpha <- alpha.temp 
          }
          if (ii == maxHalf)
                break
          if (retcode < 2)
                break 
       } ## end of REPEAT for the step length selection

      ## rho1 is now new point. 

      b <- fn(rho=rho1,G, gamma=1,ee,eb,esp,SVD=TRUE,check.analytical,del)
  
      step <- alpha*Nstep
      old.score <-score
      old.rho <- rho
      rho <- rho1
      old.grad <- grad
      score <- score1
      grad <- b$dgcv.ubre 

  score.plot[i+1] <- score
 # rho.plot[i+1] <- rho  

      ## update B...
      yg <- grad - old.grad
      rr <- 1/sum(yg * step)
      skipupdate <- TRUE 

      ## skip update if `step' is sufficiently close to B%*%yg ...
      for (i in 1:n.pen){
           closeness <- step[i]-B[i,]%*%yg  
           if (abs(closeness) >= gradtol*max(abs(grad[i]),abs(old.grad[i])))   
               skipupdate <- FALSE
        }
      ## skip update if curvature condition is not satisfied...
      if (!curv.condition)
               skipupdate <- TRUE

      if (!skipupdate) {
               B <- B - rr * step %*% (t(yg)%*%B) 
               B <- B - rr*(B %*% yg) %*% t(step) + rr * step %*% t(step)
           }
 
      
      ## check the termination condition ...
      termcode <- 0
      if (retcode ==1) 
               termcode <- 3
      else if (max(abs(grad)*max(abs(rho),1/Sp)/max(abs(score),typf))<= gradtol)
               termcode <- 1
      else if (max(abs(rho-old.rho)/max(abs(rho),1/Sp))<= steptol)
               termcode <- 2
      else if (i==200) 
               termcode <- 4 
      else if (maxtaken){ ## step of length maxNstep was taken
                consecmax <- consecmax +1
                if (consecmax ==5)
                     termcode <- 5 # limit of 5 maxNsteps was reached
            }
      else consecmax <- 0 
      

     #if (ii == maxHalf)
     #      converged <- TRUE
      if (termcode > 0)
           break
      else { ## if not converged...
            converged <- TRUE
            score.scale <- b$scale.est + score
            unconv.ind <- abs(grad) > score.scale * gradtol
            if (sum(unconv.ind))
                  converged <- FALSE
            if (abs(old.score - score) > score.scale*gradtol) {
                if (converged)
                unconv.ind <- unconv.ind | TRUE
                converged <- FALSE
             }
      } # end of ELSE

    } ## end of the Quasi-Newton algorithm 

  ## printing why the algorithm terminated... 
  if (termcode == 1)
       ct <- "Full convergence"
  else if (termcode == 2)
       ct <- "Successive iterates within tolerance, current iterate is probably solution"
  else if (termcode == 3)
       ct <- "Last step failed to locate a lower point than old.rho"
  else if (termcode == 4)
       ct <- "Iteration limit reached"
  else if (termcode ==5)
       ct <- "Five conseqcutive steps of length maxNstep have been taken" 
 #  else if (termcode == 6)
 #      ct <- "Curvature condition couldn't be satisfied"

 list (gcv.ubre=score,rho=rho,dgcv.ubre=grad,iterations=i, B=B,conv.bfgs = ct,object = b$object,
      score.plot=score.plot[1:(i+1)], termcode=termcode, check.grad= b$check.grad,
       dgcv.ubre.check=b$dgcv.ubre.check) 
}


