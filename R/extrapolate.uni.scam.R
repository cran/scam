## an attempt to predict future data for scam...


##############################################################
### function to predict future values of the response variable in case of a single univariate shape constrained term ....
################################################################

extrapolate.uni.scam <- function(object, newdata){
## object is a scam object,
## newdata is a dataframe to be used in extrapolation
## Warning: extrapolation for "micv" and "mdcx" will flatten the smooth, the very short forecasting is recommended (parameters extensions differ from those for other smooths due to the Sigma matrix construction) 

  if (!inherits(object,"scam")) stop("extrapolate.scam can only be used to extrapolate from scam objects")
  if (!inherits(object$smooth[[1]],c("mpi.smooth","cv.smooth", "cx.smooth","micx.smooth","micv.smooth","mpd.smooth","mdcx.smooth","mdcv.smooth"))) stop("extrapolate.scam can only be used to extrapolate univariate shape constrained smooths only")

  n <- length(object$data[,1])
  m <- object$smooth[[1]]$m
  knots <- object$smooth[[1]]$knots   ## evenly spaced knots
  n.kn <- length(knots)   ## number of knots
  h <- knots[2]-knots[1]  ## length between two knots
## getting additional knots...  
  if (max(newdata[[1]]) > max(object$data[,2:ncol(object$data)])) {
      knots.add <- knots[n.kn]+h
      k <- 1
      while (max(newdata[[1]]) > max(knots.add)) {
      knots.add <- c(knots.add,knots.add[k]+h)
      k <- k+1
           }
      for (j in 1:(m+1))
            knots.add <- c(knots.add,knots.add[k]+h*j)
      n.p <- length(knots.add) ## number of additional coeff
      ## constracting model matrix...
      Xp <- splineDesign(knots=c(knots,knots.add),newdata[[1]],ord=m+2)
     }
   else { n.p <- 0;
          Xp <- splineDesign(knots=knots,newdata[[1]],ord=m+2)
     }

 ## getting spline coefficients....
  beta.t <- object$coefficients.t
  q <- length(beta.t)
  qp <- q+n.p
  beta.p <- rep(beta.t[q],n.p) ## additional predicted coeffs 
  Sig <- matrix(0,q+n.p,q+n.p)   # Define Matrix Sigma
  if (inherits(object$smooth[[1]],"mpi.smooth")){
         for (i in 1:qp)  Sig[i,1:i]<-1
  } else if (inherits(object$smooth[[1]],"micv.smooth")){
         Sig[,1] <- rep(1,qp)
         Sig[2:q,2:q] <- object$smooth[[1]]$Sig
         for (i in (q+1):qp)  Sig[i,1:q] <- Sig[q,1:q];
     #    for (i in (q+1):qp) Sig[i,(q+1)]<-1
                                 # Sig[i,(q+1):i]<-1
         }
    else if (inherits(object$smooth[[1]],"micx.smooth")){
         Sig1 <- matrix(0,(qp-1),(qp-1))   
         for (i in 1:(qp-1)) Sig1[i:(qp-1),i] <- c(1:(qp-i))
         Sig[,1] <- rep(1,qp)
         Sig [2:qp,2:qp] <- Sig1
         }
    else if (inherits(object$smooth[[1]],"mpd.smooth")){
         Sig[,1] <- 1
         for (i in 2:qp)  Sig[i,2:i]<- -1
         }
    else if (inherits(object$smooth[[1]],"mdcv.smooth")){
         Sig1 <- matrix(0,(qp-1),(qp-1))   
         for (i in 1:(qp-1)) Sig1[i:(qp-1),i] <- -c(1:(qp-i))
         Sig[,1] <- rep(1,qp)
         Sig [2:qp,2:qp] <- Sig1 
         }
    else if (inherits(object$smooth[[1]],"mdcx.smooth")){
         Sig[,1] <- rep(1,qp)
         Sig[2:q,2:q] <- object$smooth[[1]]$Sig
         for (i in (q+1):qp)  Sig[i,1:q] <- Sig[q,1:q];
        # for (i in (q+1):qp)  Sig[i,(q+1):i]<- -1
         }
    else if (inherits(object$smooth[[1]],"cv.smooth")){
         Sig1 <- matrix(0,(qp-1),(qp-1)) 
         Sig1[1:(qp-1),1]<- c(1:(qp-1))  
         for (i in 2:(qp-1)) Sig1[i:(qp-1),i]<--c(1:(qp-i)) 
         Sig[,1] <- rep(1,qp)
         Sig [2:qp,2:qp] <- Sig1
         }
    else if (inherits(object$smooth[[1]],"cx.smooth")){
         Sig1 <- matrix(0,(qp-1),(qp-1)) 
         Sig1[1:(qp-1),1]<- -c(1:(qp-1))  
         for (i in 2:(qp-1)) Sig1[i:(qp-1),i]<- c(1:(qp-i)) 
         Sig[,1] <- rep(1,qp)
         Sig [2:qp,2:qp] <- Sig1
         }

 
  beta.p <- c(beta.t,beta.p)
  XpSig <- Xp%*%Sig
##  gamma <- Sig%*%beta.p ## re-parametrized coefficients
  f.p <- XpSig%*%beta.p  ## extrapolated values

## get se....
  first<-object$smooth[[1]]$first.para;
  last<-object$smooth[[1]]$last.para
  Vp <- matrix(0,qp,qp)
  Vp[1:q,1:q] <- object$Vp.t[c(1,first:last),c(1,first:last)] 
  Vp[,1] <- rep(0,nrow(Vp))
  Vp[1,] <- rep(0,ncol(Vp))
  if (qp >q){
     Vp[q:qp,q:qp] <- Vp[q,q]
     for (j in 1:(q-1)) Vp[j,(q+1):qp] <- Vp[(q+1):qp,j] <- rep(Vp[j,q],(qp-q))
    }

  if (nrow(Xp)==1) # prediction vector if prediction is made for only one value of covariates
        X1 <- c(1,t(XpSig[,first:qp]))
  else 
        X1 <- cbind(rep(1,nrow(XpSig)),XpSig[,first:qp]  ) # prediction matrix

   # X0 - model matrix of the original data....
   X0 <- cbind(rep(1,nrow(object$X)),object$X[,first:last]) 
   onet <- matrix(rep(1,nrow(X0)),1,nrow(X0))
   A <- onet%*%X0
   qrX <- qr(X0)
   R <- qr.R(qrX) # matrix R from QR decomposition
   qrA <- qr(t(A))
   R <- R[-1,]  ## (q-1) by q matrix
   RZa <- t(qr.qty(qrA,t(R)))[,2:q] ## (q-1) by (q-1) matrix
   RZa.inv <- solve(RZa)
   RZaR <- RZa.inv%*%R ## (q-1) by q matrix
   Za <- qr.Q(qrA,complete=T)[,-1] 
   ZaRZaR <- Za%*%RZaR 
  ## ?? ZaRZaR <- qr.qy(qrA,t(RZaR))

    Gp <- matrix(0,qp,qp)
   Gp[1:q,1:q] <- ZaRZaR
   if (qp >q){
      for (j in (q+1):qp) Gp[j,1:q] <- ZaRZaR[q,]
     }
   XSigGp <- X1%*%Gp
   se <- sqrt(rowSums((XSigGp%*%Vp)*XSigGp))*sqrt(1+1/n)

  list(fit=f.p,se=se)
}












