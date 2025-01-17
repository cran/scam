## (c) Natalya Pya (2024). Provided under GPL 2.
## routines for univariate local SCOP-spline construction, LSOP-spline, 
## in collaboration with Jens Lichter and Thomas Kneib, University of Gottingen


#####################################################
### Local monotone increasing LSCOP-spline construction, increasing up until the change point, xc, 
######################################################
## building B-spline bases functions over the whole range of x
## and making m knots at xc change point ...

smooth.construct.lmpi.smooth.spec<- function(object, data, knots)
## construction of the monotone increasing smooth
## xc specifies a change point
## scop-spline up until xc and unconstrained p-spline from xc, 
## using ceiling(q/2) basis functions for both parts
{ 
  m <- object$p.order[1]
  if (is.na(m)) m <- 2 ## default for cubic spline
  if (m<1) stop("silly m supplied")
  if (object$bs.dim<0) object$bs.dim <- 10 ## default
  q <- object$bs.dim 
 ## q1 <- q2 <- ceiling(q/2) ## basis dimention for each part 
  xc <- object$xt$xc ## a change point
  rg <- range(x)
  share <- (xc-rg[1])/(rg[2]-rg[1]) ## share of the x values up untill the change point
  q1 <- max(ceiling(q*share),5) ## basis dimention for the constrained part 
  q2 <- max((q-q1),5)   
  n <- length(x)
  nk <- q1+q2+1 ## number of knots 
  if (nk<=0) stop("k too small for m")
  x <- data[[object$term]]  ## the data
  xk <- knots[[object$term]] ## will be NULL if none supplied
  if (!is.null(xk)) ## if not NULL 
   stop(paste("'lmpi' smooth currenly does not work with user-supplied knots"))
  if (is.null(xc))  
   stop(paste("a change point 'xc' is not supplied"))

  ## getting equally spaced knots from both sides of the change point...
  xk <- rep(NA,q1+q2+1+1)
  xk[(m+2):(q1+1)] <- seq(rg[1],xc,length=q1-m) ## inner knots for the  constrained part (to the left from the change point)
  xk[(q1+1):(q1+q2-m+1)] <- seq(xc,rg[2],length=q2-m+1) ## inner knots for the unconstrained part (to the right from the change point), add one knot more than for the constrained part, as otherwise there 5 basis functions for the constrained part and only 3 for the unconstrained 
  for (i in 1:(m+1)) ## outer knots on the left-hand-side of the x range
       xk[i] <- xk[m+2]-(m+2-i)*(xk[m+3]-xk[m+2])
  for (i in (q1+q2-m+2):(q1+q2+2)) ## outer knots on the right-hand-side of the x range
       xk[i] <- xk[q1+q2-m]+(i-q1-q2+m)*(xk[q1+2]-xk[q1+1])
  xk <- c(xk[which(xk < xc)],rep(xc,m),xk[which(xk > xc)]) ## including xc change point m times

  if (!is.null(object$point.con[[1]])) ## a point constraint is supplied?
    stop(paste("'lmpi' smooth does not work with a point constraint; use 'miso' for a start-at-zero constraint, or 'mifo' for a finish-at-zero constraint"))

  ## get model matrix...
  X <- splineDesign(xk,x,ord=m+2) 
  Sig1 <- matrix(1,q1,q1)  ## coef summation matrix
  Sig1[upper.tri(Sig1)] <- 0
  q.t <- ncol(X) ## number of coefficients of the final lscop-spline, which is (q-2)
  Sig <- matrix(0,q.t,q.t)
  Sig[1:q1,1:q1] <- Sig1
  Sig[(q1+1):q.t,(q1+1):q.t] <- diag(1,q.t-q1)
  X <- X%*%Sig

  ## applying sum-to-zero (centering) constraint...
  cmx <- colMeans(X)
  X <- sweep(X,2,cmx) ## subtract cmx from columns 
  object$X <- X # the final model matrix
  object$cmX <- cmx

  object$P <- list()
  object$S <- list()
  object$Sigma <- Sig

  if (!object$fixed) 
  {
   # P <- diff(diag(q.t-1),difference=1)
   # P <- rbind(rep(0,q.t-1),P) ## adding 1st row of zeros
   # P <- cbind(rep(0,q.t-1),P) ## adding first column of zeros
    ## making 1st order difference penalty for the constrained part and 2nd order diff-s for the unconstrained
    P <- matrix(0,q.t-3,q.t) 
    d1 <- diff(diag(q1),difference=1)
    P[2:(q1),2:(q1+1)] <- d1
    d1 <- diff(diag(q.t-q1-1),difference=2)
    P[(q1+1):(q.t-3),(q1+2):q.t] <- d1

    object$P[[1]] <- P
    object$S[[1]] <- crossprod(P)
    ## create a block diagonal penalty matrix, penalizing separately constrained and unconstrained parts
    ## block diagonal penalty showed slightly worse mse when compared on two simul. examples ...
     ## penalty for the constrained first part
  #  P <- matrix(0,q.t-2,q.t) # matrix(0,q.t-1,q.t) 
  #  d1 <- diff(diag(q1-1),difference=1)
  #  P[2:(q1-1),2:q1] <- d1
  #   # P[q1:(q.t-1),(q1+1):q.t] <-diag(1,(q.t-q1))
  #     P[q1:(q.t-2),(q1+1):(q.t-1)] <-diag(1,(q.t-q1-1))
  #  object$P[[1]] <- P
  #  object$S[[1]] <- crossprod(P)
  #   ## block for the unconstrained second part...
  #  P <- matrix(0,q.t-2,q.t) # matrix(0,q.t-1,q.t) 
  #  d1 <- diff(diag(q.t-q1),difference=2)
  #  P[(q1+1):(q.t-2),(q1+1):q.t] <- d1
  #  object$P[[2]] <- P
  #  object$S[[2]] <- crossprod(P)
  }
  object$p.ident <- c(FALSE,rep(TRUE,q1-1),rep(FALSE,q.t-q1)) ## p.ident is an indicator of which coefficients must be positive (exponentiated)

  object$rank <- ncol(object$X)  # penalty rank
  object$null.space.dim <- 2 ##  ##m+1  # dim. of unpenalized space, 2 as the basis of a straight line is two-dimensional
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  object$knots <- xk;
  object$m <- m;
  object$df<-ncol(object$X)     # maximum DoF 
  object$q1 <- q1 

  ## get model matrix for 1st and 2nd derivatives of the smooth...
  h <- xk[m+3]-xk[m+2] ## distance between two adjacent knots
  object$Xdf1 <- splineDesign(xk,x,ord=m+1)[,1:(q.t-1)]/h ## ord is by one less for the 1st derivative
  object$Xdf2 <- splineDesign(xk,x,ord=m)[,1:(q.t-2)]/h^2 ## ord is by two less for the 2nd derivative

  class(object)<-"lmpi.smooth"  # Give object a class
  object
}


## Prediction matrix for the `lmpi` smooth class... 
Predict.matrix.lmpi.smooth<-function(object,data)
## prediction method function for the `mpi' smooth class
{ m <- object$m # spline order, m+1=3 default for cubic spline
  q <- object$df 
  q1 <- object$q1 ## basis dimention for the constrained part 
  x <- data[[object$term]]
 
  Sig1 <- matrix(1,q1,q1)  ## coef summation matrix
  Sig1[upper.tri(Sig1)] <-0
  Sig <- matrix(0,q,q)
  Sig[1:q1,1:q1] <- Sig1
  Sig[(q1+1):q,(q1+1):q] <- diag(1,q-q1)
 
  ## find spline basis inner knot range...
  ll <- object$knots[m+2];ul <- object$knots[length(object$knots)-m-1]
  n <- length(x)
  ind <- x<=ul & x>=ll ## data in range
  if (sum(ind)==n)  ## all in range
     X <- spline.des(object$knots,x,m+2)$design     
   else { ## some extrapolation needed 
     ## matrix mapping coefs to value and slope at end points...
     D <- spline.des(object$knots,c(ll,ll,ul,ul),m+2,c(0,1,0,1))$design
     X <- matrix(0,n,ncol(D)) ## full predict matrix
     if (sum(ind)> 0)  X[ind,] <- spline.des(object$knots,x[ind],m+2)$design ## interior rows
     ## Now add rows for linear extrapolation...
     ind <- x < ll 
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ll)%*%D[1:2,]
     ind <- x > ul
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ul)%*%D[3:4,]
  }
  X <- X%*%Sig 
  X <- sweep(X,2,object$cmX)
  X
}


#####################################################
### Local monotone increasing LSCOP-spline construction:
### increasing up until the change point, xc, and reaching a plateau from xc
######################################################


smooth.construct.lipl.smooth.spec<- function(object, data, knots)
## construction of the monotone increasing smooth
## xc specifies a change point
## scop-spline up until xc and a plateau from xc, 
## using ceiling(q/2) basis functions for both parts
{ 
  m <- object$p.order[1]
  if (is.na(m)) m <- 2 ## default for cubic spline
  if (m<1) stop("silly m supplied")
  if (object$bs.dim<0) object$bs.dim <- 10 ## default
  q <- object$bs.dim 
 ## q1 <- q2 <- ceiling(q/2) ## basis dimention for each part 
  xc <- object$xt$xc ## a change point
  rg <- range(x)
  share <- (xc-rg[1])/(rg[2]-rg[1]) ## share of the x values up untill the change point
  q1 <- max(ceiling(q*share),5) ## basis dimention for the constrained part 
  q2 <- max((q-q1),5)   
  n <- length(x)
  nk <- q1+q2+1 ## number of knots 
  if (nk<=0) stop("k too small for m")
  x <- data[[object$term]]  ## the data
  xk <- knots[[object$term]] ## will be NULL if none supplied
  if (!is.null(xk)) ## if not NULL 
   stop(paste("'lipl' smooth currenly does not work with user-supplied knots"))
  if (is.null(xc))  
   stop(paste("a change point 'xc' is not supplied"))

  ## getting equally spaced knots from both sides of the change point...
  xk <- rep(NA,q1+q2+1+1)
  xk[(m+2):(q1+1)] <- seq(rg[1],xc,length=q1-m) ## inner knots for the  constrained part (to the left from the change point)
  xk[(q1+1):(q1+q2-m+1)] <- seq(xc,rg[2],length=q2-m+1) ## inner knots for the unconstrained part (to the right from the change point), add one knot more than for the constrained part, as otherwise there 5 basis functions for the constrained part and only 3 for the unconstrained 
  for (i in 1:(m+1)) ## outer knots on the left-hand-side of the x range
       xk[i] <- xk[m+2]-(m+2-i)*(xk[m+3]-xk[m+2])
  for (i in (q1+q2-m+2):(q1+q2+2)) ## outer knots on the right-hand-side of the x range
       xk[i] <- xk[q1+q2-m]+(i-q1-q2+m)*(xk[q1+2]-xk[q1+1])
  xk <- c(xk[which(xk < xc)],rep(xc,m),xk[which(xk > xc)]) ## including xc change point m times

  if (!is.null(object$point.con[[1]])) ## a point constraint is supplied?
    stop(paste("'lipl' smooth does not work with a point constraint; use 'miso' for a start-at-zero constraint, or 'mifo' for a finish-at-zero constraint"))

  ## get model matrix...
  X <- splineDesign(xk,x,ord=m+2) 
  Sig1 <- matrix(1,q1,q1)  ## coef summation matrix
  Sig1[upper.tri(Sig1)] <- 0
  q.t <- ncol(X) ## number of coefficients of the final lscop-spline, which is (q-2)
  Sig <- matrix(0,q.t,q.t)
  Sig[1:q1,1:q1] <- Sig1
 ## Sig[(q1+1):q.t,(q1+1):q.t] <- diag(1,q.t-q1)
  X <- X%*%Sig
  X <- X[,-c((q1+1):q.t)]
    
  ## applying sum-to-zero (centering) constraint...
  cmx <- colMeans(X)
  X <- sweep(X,2,cmx) ## subtract cmx from columns 
  object$X <- X # the final model matrix
  object$cmX <- c(cmx, rep(0,q.t-q1))

  object$P <- list()
  object$S <- list()
  object$Sigma <- Sig

  if (!object$fixed) 
  {
    ## making 1st order difference penalty for the constrained part and...
  ##  P <- matrix(0,q.t-3,q.t) 
  ##  d1 <- diff(diag(q1),difference=1)
  ##  P[2:(q1),2:(q1+1)] <- d1
  ##  d1 <- diff(diag(q.t-q1-1),difference=2)
  ##  P[(q1+1):(q.t-3),(q1+2):q.t] <- d1
    P <- matrix(0,q1-1,q1) 
    d1 <- diff(diag(q1-1),difference=1)
    P[2:(q1-1),2:(q1)] <- d1
    object$P[[1]] <- P
    object$S[[1]] <- crossprod(P)
  }
  object$p.ident <- c(FALSE,rep(TRUE,q1-1)) ## p.ident is an indicator of which coefficients must be positive (exponentiated)
  object$n.zero.col <- q.t-q1 ## number of zeroed coeff/columns removed, needed to be added in predict and plot functions
  object$rank <- ncol(object$X)  # penalty rank
  object$null.space.dim <- 2 ##  ##m+1  # dim. of unpenalized space, 2 as the basis of a straight line is two-dimensional
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  object$knots <- xk;
  object$m <- m;
  object$df<-ncol(object$X)     # maximum DoF 
  object$q1 <- q1 

  class(object)<-"lipl.smooth"  # Give object a class
  object
}


## Prediction matrix for the `lipl` smooth class... 
Predict.matrix.lipl.smooth<-function(object,data)
## prediction method function for the `mpi' smooth class
{ m <- object$m # spline order, m+1=3 default for cubic spline
  q <- object$df 
  q1 <- object$q1 ## basis dimention for the constrained part 
  x <- data[[object$term]] 
  
  ## find spline basis inner knot range...
  ll <- object$knots[m+2];ul <- object$knots[length(object$knots)-m-1]
  n <- length(x)
  ind <- x<=ul & x>=ll ## data in range
  if (sum(ind)==n)  ## all in range
     X <- spline.des(object$knots,x,m+2)$design     
   else { ## some extrapolation needed 
     ## matrix mapping coefs to value and slope at end points...
     D <- spline.des(object$knots,c(ll,ll,ul,ul),m+2,c(0,1,0,1))$design
     X <- matrix(0,n,ncol(D)) ## full predict matrix
     if (sum(ind)> 0)  X[ind,] <- spline.des(object$knots,x[ind],m+2)$design ## interior rows
     ## Now add rows for linear extrapolation...
     ind <- x < ll 
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ll)%*%D[1:2,]
     ind <- x > ul
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ul)%*%D[3:4,]
  }
  Sig1 <- matrix(1,q1,q1)  ## coef summation matrix
  Sig1[upper.tri(Sig1)] <- 0
  q.t <- ncol(X) ## number of coefficients of the final lscop-spline, which is (q-2)
  Sig <- matrix(0,q.t,q.t)
  Sig[1:q1,1:q1] <- Sig1
  X <- X%*%Sig 
  X <- sweep(X,2,object$cmX)
  X
}















