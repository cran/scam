
## BFGS for GCV/UBRE minimization 
## with modifications on convergence control, scaling, etc...


##############################################################
## function to get the gcv/ubre value and their gradient    ##
##############################################################

gcv.ubre_grad <- function(rho, G, env, control)
 { ## G - object from scam setup via gam(..., fit=FALSE) 
   ## rho - logarithm of the smoothing parameters
   ## gamma - an ad hoc parametrer of the GCV
   ## check.analytical - logical whether the analytical gradient of GCV/UBRE should be checked
   ## del - increment for finite differences when checking analytical gradients
   y <- drop(G$y);  X <- G$X; gamma <- G$gamma
   S <- G$S;
   not.exp <- G$not.exp
   q0 <- G$q0; q.f <- G$q.f
   p.ident <- G$p.ident; n.terms <- G$n.terms
   family <- G$family; intercept <- G$intercept; offset <- G$offset;
   weights <- G$weights;  scale.known <- G$scale.known; sig2 <- G$sig2
   n.pen <- length(S) # number of penalties
   if (length(rho)!=n.pen) stop (paste("length of rho and # penalties has to be the same"))
   sp <- exp(rho)
   ## fit the model with the given values of the smoothing parameters...
   b <- scam.fit(G=G,sp=sp, env=env, control=control) 
   n <- nobs <- nrow(G$X)
   q <- ncol(G$X)
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
        X <- rwMatrix(end,row,weight.r,G$X)
        y <- rwMatrix(end,row,weight.r,G$y)    
      }   


   y.mu <- y - b$mu##  b$y-b$mu
   c <- -2*y.mu/(b$Var*b$dlink.mu) 
  ## diag.C <- rep(1,q)   # diagonal elements of matrix C with estimates of beta at convergence
  ## diag.C[b$iv] <- b$beta.t[b$iv]
  ## diag.C <- b$Cdiag
  ## D.beta <- (diag.C*t(b$X))%*%c  # derivative of the deviance w.r.t. beta 
   D.beta <- t(b$X1)%*%c  # derivative of the deviance w.r.t. beta 
   ## -----------------------------------------------------------
   ## calculation of the deviance derivative wrt rho
   D.rho <- rep(0,n.pen)         # define derivative of the deviance wrt rho
   D.rho <- t(D.beta)%*%b$dbeta.rho
   ## ------------------------------------------------------------
   ## calculating the derivative of tr(A) w.r.t. log(sp)
   d2link.dlink <- b$d2link.mu/b$dlink.mu
   a1 <- as.numeric(y.mu*d2link.dlink)    # a constant for updating the derivative of z
   a2 <- as.numeric(b$w^2*(b$dvar.mu*b$dlink.mu+2*b$Var*b$d2link.mu)) # a constant for updating the derivative of w
   eta.rho <- matrix(0,n,n.pen) # define derivatives of the linear predictor
   N_rho <- matrix(0,q,n.pen) # define diagonal elements of N_j
   w.rho <- matrix(0,n,n.pen) # define derivatives of the diagonal elements of W
   alpha.rho <- matrix(0,n,n.pen) # define derivatives of alpha
   T_rho <- matrix(0,n,n.pen) # define diagonal elements of T_j
   ## a constant for the derivative of the alpha
   dvar.var <- b$dvar.mu/b$Var
   alpha1 <- as.numeric(-(dvar.var+d2link.dlink)/b$dlink.mu -
                        y.mu*(dvar.var^2 + d2link.dlink^2 - b$d2var.mu/b$Var -
                                  b$d3link.mu/b$dlink.mu)/b$dlink.mu)
   eta.rho <- b$X1%*%b$dbeta.rho
   N_rho[b$iv,] <- b$dbeta.rho[b$iv,]
   alpha.rho <- alpha1*eta.rho
   w.rho <- -a2*b$alpha*eta.rho + b$w1*alpha.rho
   T_rho <- w.rho /b$abs.w
   z2 <- b$dlink.mu*y.mu    # term for the derivative of E
   w1.rho <- matrix(0,n,n.pen)  # define derivatives of the diagonal elements of W1
   T1_rho <- matrix(0,n,n.pen)  # define diagonal elements of T1_j
   Q_rho <- matrix(0,q,n.pen)   # derivative of E diagonal wrt log(sp[j]) 
   w1.rho <- -a2*eta.rho
   T1_rho <- w1.rho/b$w1
   term <- T1_rho*z2 + a1*eta.rho- eta.rho
   Q_rho <- N_rho*drop(b$C2diag*crossprod(X,b$w1*z2)) +  
           b$C1diag*crossprod(X,b$w1*term)    
   ## Q_rho <- N_rho*drop(b$C1diag*(t(b$X)%*%(b$w1*z2))) +
      ##     b$C1diag*(t(b$X)%*%(b$w1*(T1_rho*z2 + a1*eta.rho- eta.rho)))
   ## efficient version of derivative of trA...
   KtIL <- t((b$L*b$I.plus)*b$K)
   KtILK <- KtIL%*%b$K
   KKtILK <- b$K%*%KtILK
 ##  KtIQ1R <-  if (!not.exp) crossprod(b$I.plus*b$K,b$wX1) else crossprod(b$I.plus*b$K,b$wXC1)  
 ##  KtILQ1R<-if (!not.exp) crossprod(b$L*b$I.plus*b$K,b$wX1) else crossprod(b$L*b$I.plus*b$K,b$wXC1)
   trA.rho<-rep(0,n.pen)
   #if (n.pen>0) 
                for (j in 1:n.pen){
                   trA.rho[j] <- - 2*sum(KtILK*(b$KtIQ1R%*%(N_rho[,j]*b$P))) -
                     sum((T_rho[,j]*KKtILK)*b$K) -  
                     sp[j]*sum((t(b$P)%*%b$S[[j]]%*%b$P)*t(KtILK) )+
                     sum( (t(b$P)%*%(c(Q_rho[,j])*b$P))*t(KtILK) ) +
                     2*sum( b$KtILQ1R*t(N_rho[,j]*b$P) ) +
                     #2*sum(KtILQ1R*t(N_rho[,j]*b$P) ) +
                     sum(KtIL*t(T1_rho[,j]*b$K))
                 }
   ## Calculating the derivatives of the trA is completed here -----
   if (scale.known) ##  derivative of Mallow's Cp/UBRE/AIC wrt log(sp) 
      { ubre.rho <- rep(0,n.pen) 
        #if (n.pen>0) 
                     for (j in 1:n.pen){
                        ubre.rho[j] <- D.rho[j]/n +2*gamma*trA.rho[j]*sig2/n
                     }  
        gcv.ubre.rho <- ubre.rho
        gcv.ubre <- b$dev/n - sig2 +2*gamma*b$trA*sig2/n
      }
   else { # derivative of GCV wrt log(sp) ...
          gcv.rho<-rep(0,n.pen) # define the derivatives of the gcv
          #if (n.pen>0) 
                      for (j in 1:n.pen){
                          gcv.rho[j]<-n*(D.rho[j]*(n-gamma*b$trA)+2*gamma*b$dev*trA.rho[j])/(n-gamma*b$trA)^3
                       }  
          gcv.ubre <- b$gcv
          gcv.ubre.rho <- gcv.rho
        }

   ## checking the derivatives by finite differences-----------
   dgcv.ubre.check <- NULL
   if (control$bfgs$check.analytical){
        ##del <- del #1e-4
         sp1 <- rep(0,n.pen)
         dbeta.check <- matrix(0,q,n.pen)
         dtrA.check <- rep(0,n.pen)
         dgcv.ubre.check <- rep(0,n.pen)
         #if (n.pen>0) 
           for (j in 1:n.pen){
               sp1 <- sp; sp1[j] <- sp[j]*exp(control$bfgs$del)
               b1 <- scam.fit(G=G,sp=sp1, env=env,control=control) 
               ## calculating the derivatives of beta estimates by finite differences...
               dbeta.check[,j] <- (b1$beta - b$beta)/control$bfgs$del
               ## calculating the derivatives of the trA by finite differences...
               dtrA.check[j] <- (b1$trA-b$trA)/control$bfgs$del
               ## calculating derivatives of GCV/UBRE by finite differences...
               if (scale.known) gcv.ubre1 <- b1$dev/n - sig2 + 2*gamma*b1$trA*sig2/n
               else gcv.ubre1 <- b1$gcv
               dgcv.ubre.check[j] <- (gcv.ubre1-gcv.ubre)/control$bfgs$del
           }
   }
   check.grad <- NULL
   if (control$bfgs$check.analytical)
          check.grad <- 100*(gcv.ubre.rho-dgcv.ubre.check)/dgcv.ubre.check
   ## end of checking the derivatives ----------------------------
   list(dgcv.ubre=gcv.ubre.rho,gcv.ubre=gcv.ubre, scale.est=b$dev/(n-b$trA), 
         check.grad=check.grad, dgcv.ubre.check=dgcv.ubre.check, object = b, trA.rho=trA.rho,D.rho=D.rho)
 } ## end of gcv.ubre_grad



#####################################################
## BFGS for gcv/ubre miminization..                ##
#####################################################

bfgs_gcv.ubre <- function(fn=gcv.ubre_grad, rho, ini.fd=TRUE, G, env,
              n.pen=length(rho), typx=rep(1,n.pen), typf=1, control)
             ## steptol.bfgs= 1e-7, gradtol.bfgs = 1e-06, maxNstep = 5, maxHalf = 30, check.analytical, del, devtol.fit, steptol.fit)  
{  ## fn - GCV/UBRE Function which returns the GCV/UBRE value and its derivative wrt log(sp)
   ## rho - log of the initial values of the smoothing parameters
   ## ini.fd - if TRUE, a finite difference to the Hessian is used to find the initial inverse Hessian
   ## typx - vector whose component is a positive scalar specifying the typical magnitude of sp
   ## typf - a positive scalar estimating the magnitude of the gcv near the minimum
   ## gradtol.bfgs - a scalar giving a tolerance at which the gradient is considered
   ##            to be close enougth to 0 to terminate the BFGS algorithm
   ## steptol.bfgs - a positive scalar giving the tolerance at which the scaled distance between
   ##          two successive iterates is considered close enough to zero to terminate the BFGS algorithm 
   ## maxNstep - a positive scalar which gives the maximum allowable step length
   ## maxHalf - a positive scalar which gives the maximum number of step halving in "backtracking"
   ## check.analytical - logical whether the analytical gradient of GCV/UBRE should be checked
   ## del - increment for finite differences when checking analytical gradients
   ## devtol.fit, steptol.fit - scalars giving the tolerance for the full Newton methods to estimate model coefficients
   
   Sp <- 1/typx     # diagonal of the scaling matrix 
   maxNstep <- control$bfgs$maxNstep ## the maximum allowable step length
   ## storage for solution track
   rho1 <- rho
   old.rho <- rho
   not.exp <- G$not.exp
   b <- fn(rho,G, env, control=control) 
   old.score <- score <- b$gcv.ubre

   max.step <- 200
   score.hist <- rep(NA,max.step+1) ## storing scores for history, for plotting the gcv
   score.hist[1] <- score
   grad <- b$dgcv.ubre
   scale.est <- b$scale.est
   rm(b)

   ## The initial inverse Hessian ...
   if (ini.fd) 
      {  B <- matrix(0,n.pen,n.pen)
         feps <- 1e-4
         #if (n.pen>0) 
                    for (j in 1:n.pen){ 
                        rho2 <- rho; rho2[j] <- rho[j] + feps
                        b2 <- fn(rho2,G,env, control=control) 
                        B[,j] <- (b2$dgcv.ubre - grad)/feps  
                        rm(b2) 
                      }
         ## force initial Hessian to +ve def and invert... 
         B <- (B+t(B))/2 ## B + t(B)
         eh <- eigen(B,symmetric=TRUE)
         eh$values <- abs(eh$values) ## eh$values[eh$values < 0] <- -eh$values[eh$values < 0] 
         ##ind <- eh$values > max(eh$values)*.Machine$double.eps^75 ## index of non zero eigenvalues
         ##eh$values[ind] <- 1/eh$values[ind]; eh$values[!ind] <- 0 
         ##B <- eh$vectors%*%(eh$values*t(eh$vectors)) 
         thresh <- max(eh$values) * 1e-4
         eh$values[eh$values<thresh] <- thresh  
         B <- eh$vectors%*%(t(eh$vectors)/eh$values)       
      } 
   else B <- diag(n.pen)*100

   uconv.ind <- rep(TRUE,ncol(B))
   ## Wolfe condition constants...
   c1 <- 1e-4   ## Wolfe 1, constant for the sufficient decrease condition 
   c2 <- .9   ## Wolfe 2, curvature condition  
 
   score.scale <- abs(scale.est) + abs(score)  ## ??? for UBRE
   unconv.ind <- abs(grad) > score.scale*control$bfgs$gradtol.bfgs # checking the gradient is within the tolerance
   if (!sum(unconv.ind))  ## if at least one is false
         unconv.ind <- unconv.ind | TRUE
   consecmax <- 0

   ## Quasi-Newton algorithm to minimize GCV...
   for (i in 1:max.step) {
        ## compute a BFGS search direction ...
         Nstep <- 0*grad   ## initialize the quasi-Newton step
         Nstep[unconv.ind] <- -drop(B[unconv.ind, unconv.ind]%*%grad[unconv.ind])
 
         if (sum(Nstep*grad)>=0) { ## step not descending!
           ## Following would really be in the positive definite space... 
           ##step[uconv.ind] <- -solve(chol2inv(chol(B))[uconv.ind,uconv.ind],initial$grad[uconv.ind])
           Nstep <- -diag(B)*grad ## simple scaled steepest descent 
           Nstep[!uconv.ind] <- 0 ## don't move if apparently converged 
         }

         Dp <- Sp*Nstep
         Newtlen <- (sum(Dp^2))^.5  ## Euclidean norm to get the Newton step
         if (Newtlen > maxNstep) { ## reduce if max step is greater than the max allowable
              Nstep <- maxNstep*Nstep/Newtlen
              Newtlen <- maxNstep
         }
         maxtaken <- FALSE
         retcode <- 2

         initslope <- sum(Nstep* grad) ## initial slope, directional derivative
         rellength <- max(abs(Nstep)/max(abs(rho),1/Sp)) ## relative length of rho for the stopping criteria
         alpha.min <- control$bfgs$steptol.bfgs/rellength
         
         alpha.max <- maxNstep/Newtlen

         ms <- max(abs(Nstep))
         if (ms > maxNstep) { ## initialize step length, alpha
            alpha <- maxNstep/ms  
            alpha.max <- alpha*1.05
         } else {
              alpha <- 1 
              alpha.max <- min(2,maxNstep/ms) ## 1*maxNstep/ms
           }
    
         ##alpha <- 1   ## initialize step length
         ii <- 0 ## initialize the number of "step halving"
         step <- alpha*Nstep 
         ## step length selection ...
         curv.condition <- TRUE
         repeat 
             {  rho1 <- rho + alpha*Nstep 
                b <- fn(rho=rho1,G,env, control=control) 
                score1 <- b$gcv.ubre
                if (score1 <= score+c1*alpha*initslope) {## check the first Wolfe condition (sufficient decrease)...
                       grad1 <- b$dgcv.ubre ## Wolfe 1 met
                       newslope <- sum(grad1 * Nstep) ## directional derivative
                       curv.condition <- TRUE
                       if (newslope < c2*initslope) {# the curvature condition (Wolfe 2) is not satisfied
                              if (alpha == 1 && Newtlen < maxNstep)
                                 {  ## alpha.max <- maxNstep/Newtlen
                                     repeat 
                                        {  old.alpha <- alpha
                                           old.score1 <- score1
                                           alpha <- min(2*alpha, alpha.max)
                                           rho1 <- rho + alpha*Nstep
                                           b <- fn(rho=rho1,G, env, control=control) 
                                           score1 <- b$gcv.ubre
                                           if (score1 <= score+c1*alpha*initslope)
                                              {   grad1 <- b$dgcv.ubre
                                                  newslope <- sum(grad1*Nstep)
                                              }
                                           if (score1 > score+c1*alpha*initslope)
                                               break
                                           if (newslope >= c2*initslope)
                                               break
                                           if (alpha >= alpha.max)
                                               break
                                         } 
                                  }
                               if ((alpha < 1) || (alpha>1 && (score1>score+c1*alpha*initslope)))
                                  {   alpha.lo <- min(alpha, old.alpha)
                                      alpha.diff <- abs(old.alpha - alpha)
                                      if (alpha < old.alpha)
                                         {  sc.lo <- score1
                                            sc.hi <- old.score1
                                          }    
                                      else
                                          {  sc.lo <- old.score1
                                             sc.hi <- score1 
                                          }
                                      repeat
                                         {  alpha.incr <- -newslope*alpha.diff^2/(2*(sc.hi-(sc.lo+newslope*alpha.diff)))
                                            if (alpha.incr < .2*alpha.diff) 
                                                  alpha.incr <- .2*alpha.diff
                                            alpha <- alpha.lo+alpha.incr
                                            rho1 <- rho + alpha*Nstep
                                            b <- fn(rho=rho1,G, env, control=control) 
                                            score1 <- b$gcv.ubre
                                            if (score1 > score+c1*alpha*initslope)
                                               {  alpha.diff <- alpha.incr
                                                  sc.hi <- score1
                                                }
                                            else
                                                {   grad1 <- b$dgcv.ubre
                                                    newslope <- sum(grad1*Nstep)
                                                    if (newslope < c2 *initslope)
                                                       {  alpha.lo <- alpha
                                                          alpha.diff <- alpha.diff-alpha.incr
                                                          sc.lo <- score1
                                                        }
                                                 }
                                            if (abs(newslope) <= -c2*initslope) ## met Wolfe 2, curvature condition
                                                    break
                                            if (alpha.diff < alpha.min)
                                                    break 
                                          }
                                      if (newslope < c2*initslope)   ## couldn't satisfy curvature condition
                                         {   curv.condition <- FALSE
                                             score1 <- sc.lo
                                             rho1 <- rho + alpha.lo*Nstep
                                             b <- fn(rho=rho1,G,env, control=control)
                                          } 
                                  } ## end of "if ((alpha < 1) || (alpha>1 && ..."
                           }  ## end of "if (newslope < c2*initslope) ...", if Wolfe 2 not met
                       retcode <- 0
                       if (newslope < c2*initslope) ## couldn't satisfy curvature condition, failed Wolfe 2
                             curv.condition <- FALSE
                       if (alpha*Newtlen > .99*maxNstep)  
                             maxtaken <- TRUE
                   } ## end of "if (score1 <= ...) ...", checking the first Wolfe condition (sufficient decrease)
                else if (alpha < alpha.min) {## no satisfactory rho+ can be found suff-ly distinct from previous rho
                        retcode <- 1
                        rho1 <- rho
                        b <- fn(rho=rho1,G,env, control=control)
                     }
                else   ## backtracking to satisfy the sufficient decrease condition...
                    {   ii <- ii+1
                        if (alpha == 1) {## first backtrack, quadratic fit
                              alpha.temp <- -initslope/(score1-score-initslope)/2
                        }
                        else {  ## all subsequent backtracts, cubic fit 
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
                                if (ab[1] == 0) ## cubic is a quadratic
                                       alpha.temp <- -initslope/ab[2]/2
                                else     ## legitimate cubic
                                       alpha.temp <- (-ab[2]+disc^.5)/(3*ab[1])
                                if (alpha.temp > .5*alpha)
                                       alpha.temp <- .5*alpha   
                             }
                        old.alpha <- alpha
                        old.score1 <- score1
                        if (alpha.temp <= .1*alpha) alpha <- .1*alpha
                        else alpha <- alpha.temp 
                    }
                if (ii == control$bfgs$maxHalf)
                    break
                if (retcode < 2)
                    break 
             } ## end of REPEAT for the step length selection
         ## rho1 is now new point. 
         step <- alpha*Nstep
         old.score <-score
         old.rho <- rho
         rho <- rho1
         old.grad <- grad
         score <- score1
         grad <- b$dgcv.ubre 
         score.hist[i+1] <- score
         ## update B...
         yg <- grad - old.grad
         
         skipupdate <- TRUE 
         ## skip update if `step' is sufficiently close to B%*%yg ...
         #if (n.pen>0) 
                     for (j in 1:n.pen){
                          closeness <- step[j]-B[j,]%*%yg  
                          if (abs(closeness) >= control$bfgs$gradtol.bfgs*max(abs(grad[j]),abs(old.grad[j])))   
                     skipupdate <- FALSE
                       }
         ## skip update if curvature condition is not satisfied...
         if (!curv.condition)
               skipupdate <- TRUE
         if (!skipupdate) {
               if (i==1) { ## initial step --- adjust Hessian as p143 of N&W
                  B <- B * alpha ## this is Simon's version 
                  ## B <- B * sum(yg*step)/sum(yg*yg) ## this is N&W
               }
               rr <- 1/sum(yg * step)
               B <- B - rr * step %*% crossprod(yg,B) # (t(yg)%*%B) 
               B <- B - rr*tcrossprod((B %*% yg),step) + rr *tcrossprod(step)  # B - rr*(B %*% yg) %*% t(step) + rr * step %*% t(step)
         }
 
         ## check the termination condition ...
         termcode <- 0
         if (retcode ==1) {
               if (max(abs(grad)*max(abs(rho),1/Sp)/max(abs(score),typf))<= control$bfgs$gradtol.bfgs*6.0554)
                    termcode <- 1
               else termcode <- 3
         }
         else if (max(abs(grad)*max(abs(rho),1/Sp)/max(abs(score),typf))<= control$bfgs$gradtol.bfgs*6.0554)
               termcode <- 1
         else if (max(abs(rho-old.rho)/max(abs(rho),1/Sp))<= control$bfgs$steptol.bfgs)
               termcode <- 2
         else if (i==max.step) 
               termcode <- 4 
         else if (maxtaken) ## step of length maxNstep was taken
               {  consecmax <- consecmax +1
                  if (consecmax ==5)
                      termcode <- 5 # limit of 5 maxNsteps was reached
               }
         else consecmax <- 0 
         ##---------------------
         if (termcode > 0)
               break
         else  ## if not converged...
             {   converged <- TRUE
                 score.scale <- abs(b$scale.est) + abs(score)
                 unconv.ind <- abs(grad) > score.scale * control$bfgs$gradtol.bfgs ##*.1
                 if (sum(unconv.ind))
                        converged <- FALSE
                 if (abs(old.score - score) > score.scale*control$bfgs$gradtol.bfgs) 
                     {  if (converged)
                            unconv.ind <- unconv.ind | TRUE
                        converged <- FALSE
                     }
             } # end of ELSE
      } ## end of the Quasi-Newton algorithm 

   ## final fit...
  ## b <- fn(rho=rho,G,gamma=gamma,env, control=control)  

   ## printing why the algorithm terminated... 
   if (termcode == 1)
         ct <- "Full convergence"
   else if (termcode == 2)
         ct <- "Successive iterates within tolerance, current iterate is probably solution"
   else if (termcode == 3)
         ct <- "Last step failed to locate a lower point than 'rho'"
   else if (termcode == 4)
         ct <- "Iteration limit reached"
   else if (termcode ==5)
         ct <- "Five consecutive steps of length maxNstep have been taken" 
   list (gcv.ubre=score, rho=rho, dgcv.ubre=grad, iterations=i, B=B, conv.bfgs = ct, object=b$object, score.hist=score.hist[!is.na(score.hist)], termcode = termcode, check.grad= b$check.grad,
       dgcv.ubre.check = b$dgcv.ubre.check) 
} ## end bfgs_gcv.ubre


