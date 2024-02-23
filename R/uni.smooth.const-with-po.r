## (c) Natalya Pya (2012-2024). Provided under GPL 2.
## routines for univariate SCOP-spline construction 
## with sum-to-zero identifiability constraints (2023)



#####################################################
### Adding Monotone increasing SCOP-spline construction 
######################################################

smooth.construct.mpi.smooth.spec<- function(object, data, knots)
## construction of the monotone increasing smooth
{ 
   #require(splines)
  m <- object$p.order[1]
  if (is.na(m)) m <- 2 ## default for cubic spline
  if (m<1) stop("silly m supplied")
  if (object$bs.dim<0) object$bs.dim <- 10 ## default
  q <- object$bs.dim 
  nk <- q+m+2 ## number of knots
  if (nk<=0) stop("k too small for m")
  x <- data[[object$term]]  ## the data
  xk <- knots[[object$term]] ## will be NULL if none supplied
  if (is.null(xk)) # space knots through data
     { n<-length(x)
       xk<-rep(0,q+m+2)
       xk[(m+2):(q+1)]<-seq(min(x),max(x),length=q-m)
       for (i in 1:(m+1)) {xk[i]<-xk[m+2]-(m+2-i)*(xk[m+3]-xk[m+2])}
       for (i in (q+2):(q+m+2)) {xk[i]<-xk[q+1]+(i-q-1)*(xk[m+3]-xk[m+2])}
     }
  if (length(xk)!=nk) # right number of knots?
   stop(paste("there should be ", nk," supplied knots"))
  if (!is.null(object$point.con[[1]])) ## a point constraint is supplied?
    stop(paste("'mpi' smooth does not work with a point constraint; use 'miso' for a start-at-zero constraint, or 'mifo' for a finish-at-zero constraint"))

  #######################################################################
  ## indentifiability constraint by first dropping the first columns of XSigma (setting beta_1=0)
  ## and then applying sum-to-zero (centering) constraint...
  #  get model matrix-------------
  X1 <- splineDesign(xk,x,ord=m+2)
  # get matrix Sigma and remove the first column for the model matrix XSig
 ##  Sig <- matrix(as.numeric(rep(1:q,q)>=rep(1:q,each=q)),q,q) ## coef summation matrix
  Sig <- matrix(1,q,q)  ## coef summation matrix
  Sig[upper.tri(Sig)] <-0
  X <- X1%*%Sig
  X <- X[,-1]
  
  ## applying sum-to-zero (centering) constraint...
  cmx <- colMeans(X)
  X <- sweep(X,2,cmx) ## subtract cmx from columns 
  object$X <- X # the final model matrix
  object$cmX <- cmx

  object$P <- list()
  object$S <- list()
  object$Sigma <- Sig[-1,-1]

  if (!object$fixed) # create the penalty matrix
  { P <- diff(diag(q-1),difference=1)
    object$P[[1]] <- P
    object$S[[1]] <- crossprod(P)
  }
  object$p.ident <- rep(TRUE,q-1) ## p.ident is an indicator of which coefficients must be positive (exponentiated)

  object$rank <- ncol(object$X)-1  # penalty rank
  object$null.space.dim <- 2 ##  ##m+1  # dim. of unpenalized space, 2 as the basis of a straight line is two-dimensional
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  object$knots <- xk;
  object$m <- m;
  object$df<-ncol(object$X)     # maximum DoF 

  ## get model matrix for 1st and 2nd derivatives of the smooth...
  h <- (max(x)-min(x))/(q-m-1) ## distance between two adjacent knots
  object$Xdf1 <- splineDesign(xk,x,ord=m+1)[,2:(q-1)]/h ## ord is by one less for the 1st derivative
  object$Xdf2 <- splineDesign(xk,x,ord=m)[,2:(q-2)]/h^2 ## ord is by two less for the 2nd derivative

  class(object)<-"mpi.smooth"  # Give object a class
  object
}


## Prediction matrix for the `mpi` smooth class... 

Predict.matrix.mpi.smooth<-function(object,data)
## prediction method function for the `mpi' smooth class
{ m <- object$m # spline order, m+1=3 default for cubic spline
  q <- object$df +1
 # Sig <- matrix(as.numeric(rep(1:q,q)>=rep(1:q,each=q)),q,q) ## coef summation matrix
  Sig <- matrix(1,q,q)  ## coef summation matrix
  Sig[upper.tri(Sig)] <-0
  ## find spline basis inner knot range...
  ll <- object$knots[m+2];ul <- object$knots[length(object$knots)-m-1]
  x <- data[[object$term]]
  n <- length(x)
  ind <- x<=ul & x>=ll ## data in range
  if (sum(ind)==n) { ## all in range
     X <- spline.des(object$knots,x,m+2)$design
     X <- X%*%Sig 
  } else { ## some extrapolation needed 
     ## matrix mapping coefs to value and slope at end points...
     D <- spline.des(object$knots,c(ll,ll,ul,ul),m+2,c(0,1,0,1))$design
     X <- matrix(0,n,ncol(D)) ## full predict matrix
     if (sum(ind)> 0)  X[ind,] <- spline.des(object$knots,x[ind],m+2)$design ## interior rows
     ## Now add rows for linear extrapolation...
     ind <- x < ll 
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ll)%*%D[1:2,]
     ind <- x > ul
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ul)%*%D[3:4,]
     X <- X%*%Sig 
  ## X <- sweep(X,2,c(0,object$cmX))
  }
  X <- sweep(X,2,c(0,object$cmX))
  X
}

####################################################################################################
### Adding monotone increasing SCOP-spline construction without applying identifiability constraints
### to be used with numeric 'by' variable...
####################################################################################################
## when 'by' variable takes more than one value, the smooth terms are identifiable without a
## 'zero intercept' constraint, so they are left unconstrained.  

smooth.construct.mpiBy.smooth.spec<- function(object, data, knots)
## construction of the monotone increasing smooth
{ 
  m <- object$p.order[1]
  if (is.na(m)) m <- 2 ## default for cubic spline
  if (m<1) stop("silly m supplied")
  if (object$bs.dim<0) object$bs.dim <- 10 ## default
  q <- object$bs.dim 
  nk <- q+m+2 ## number of knots
  if (nk<=0) stop("k too small for m")
  x <- data[[object$term]]  ## the data
  xk <- knots[[object$term]] ## will be NULL if none supplied
  if (is.null(xk)) # space knots through data
     { n<-length(x)
       xk<-rep(0,q+m+2)
       xk[(m+2):(q+1)]<-seq(min(x),max(x),length=q-m)
       for (i in 1:(m+1)) {xk[i]<-xk[m+2]-(m+2-i)*(xk[m+3]-xk[m+2])}
       for (i in (q+2):(q+m+2)) {xk[i]<-xk[q+1]+(i-q-1)*(xk[m+3]-xk[m+2])}
     }
  if (length(xk)!=nk) # right number of knots?
   stop(paste("there should be ", nk," supplied knots"))
  if (!is.null(object$point.con[[1]])) ## a point constraint is supplied?
    stop(paste("'mpi' smooth does not work with a point constraint; use 'miso' for a start-at-zero constraint, or 'mifo' for a finish-at-zero constraint"))

  #  get model matrix-------------
  X1 <- splineDesign(xk,x,ord=m+2)
  # get matrix Sigma...
  Sig <- matrix(1,q,q)  ## coef summation matrix
  Sig[upper.tri(Sig)] <-0
  X <- X1%*%Sig 
  object$X <- X # the final model matrix
  object$P <- list()
  object$S <- list()
  object$Sigma <- Sig

  if (!object$fixed) # create the penalty matrix
  { P <- diff(diag(q-1),difference=1)
    P <- rbind(rep(0,q-1),P) ## adding 1st row of zeros
    P <- cbind(rep(0,q-1),P) ## adding first column of zeros
    object$P[[1]] <- P
    object$S[[1]] <- crossprod(P)
  }
  object$p.ident <- c(FALSE,rep(TRUE,q-1)) ## p.ident is an indicator of which coefficients must be positive (exponentiated)

  object$rank <- ncol(object$X)  # penalty rank
  object$null.space.dim <- 2 ## m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  object$knots <- xk;
  object$m <- m;
  object$df<-ncol(object$X)     # maximum DoF 

  ## get model matrix for 1st and 2nd derivatives of the smooth...
  h <- (max(x)-min(x))/(q-m-1) ## distance between two adjacent knots
  object$Xdf1 <- splineDesign(xk,x,ord=m+1)[,1:(q-1)]/h ## ord is by one less for the 1st derivative
  object$Xdf2 <- splineDesign(xk,x,ord=m)[,1:(q-2)]/h^2 ## ord is by two less for the 2nd derivative

  class(object)<-"mpiBy.smooth"  # Give object a class
  object
}


## Prediction matrix for the `mpiBy` smooth class... 
Predict.matrix.mpiBy.smooth<-function(object,data)
## prediction method function for the `mpiBy' smooth class
{ m <- object$m # spline order, m+1=3 default for cubic spline
  q <- object$df 
  Sig <- matrix(1,q,q)  ## coef summation matrix
  Sig[upper.tri(Sig)] <-0
  ## find spline basis inner knot range...
  ll <- object$knots[m+2];ul <- object$knots[length(object$knots)-m-1]
  x <- data[[object$term]]
  n <- length(x)
  ind <- x<=ul & x>=ll ## data in range
  if (sum(ind)==n) { ## all in range
     X <- spline.des(object$knots,x,m+2)$design
     X <- X%*%Sig 
  } else { ## some extrapolation needed 
     ## matrix mapping coefs to value and slope at end points...
     D <- spline.des(object$knots,c(ll,ll,ul,ul),m+2,c(0,1,0,1))$design
     X <- matrix(0,n,ncol(D)) ## full predict matrix
     if (sum(ind)> 0)  X[ind,] <- spline.des(object$knots,x[ind],m+2)$design ## interior rows
     ## Now add rows for linear extrapolation...
     ind <- x < ll 
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ll)%*%D[1:2,]
     ind <- x > ul
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ul)%*%D[3:4,]
     X <- X%*%Sig 
  }
  X
}


##############################################################################################
### Adding Monotone increasing SCOP-spline construction with a 'start at zero' constraint,
### a SCOP-spline additionally constrained to be zero at a start, at the left-end point of the covariate range...
#############################################################################################

smooth.construct.miso.smooth.spec<- function(object, data, knots)
## construction of the monotone increasing smooth with a 'start-at-zero' constraint;
## achieved simply by setting the first (m+1) spline coefficients to zero...
{ 
  m <- object$p.order[1]
  if (is.na(m)) m <- 2 ## default for cubic spline
  if (m<1) stop("silly m supplied")
  if (object$bs.dim<0) object$bs.dim <- 10 ## default
  q <- object$bs.dim 
  nk <- q+m+2 ## number of knots
  if (nk<=0) stop("k too small for m")
  x <- data[[object$term]]  ## the data
  xk <- knots[[object$term]] ## will be NULL if none supplied
  if (is.null(xk)) { # space knots evenly through data
       n<-length(x)
       xk<-rep(0,q+m+2)
       xk[(m+2):(q+1)]<-seq(min(x),max(x),length=q-m)
       for (i in 1:(m+1)) {xk[i]<-xk[m+2]-(m+2-i)*(xk[m+3]-xk[m+2])}
       for (i in (q+2):(q+m+2)) {xk[i]<-xk[q+1]+(i-q-1)*(xk[m+3]-xk[m+2])}
     }
  ## move m knots, xk[m+3],...,xk[m+3+m-1] (2nd,...,(m+1) inner knots), to the first inner knot xk[m+2],
  ## join these m knots with the constrained xk[m+2] knot to avoid a plateau start...
  xk[(m+3):(m+2+m)] <- xk[m+2]
  if (length(xk)!=nk) # right number of knots?
   stop(paste("there should be ", nk," supplied knots"))
  if (!is.null(object$point.con[[1]])) ## a point constraint is supplied?
    stop(paste("'miso' smooth works only with a 'start at zero' constraint; 'pc' argument does not work here"))
  ## get model matrix...
  X1 <- splineDesign(xk,x,ord=m+2)
  ## get unconstrained matrix Sigma and remove the first m+1 columns and rows...
  Sig <- matrix(1,q,q)  ## coef summation matrix
  Sig[upper.tri(Sig)] <-0
  ind <- 1:(m+1) ## 1:3 if m=2;
  ## Sig[,ind] <- 0; Sig[ind,] <- 0
  Sig <- Sig[-ind,-ind]
  X <- X1[,-ind]%*%Sig # drop (m+1) start terms, model submatrix for the scop-term
  object$X <- X # the final model matrix
  object$P <- list()
  object$S <- list()
  object$Sigma <- Sig
 
  if (!object$fixed) # create the penalty matrix
  { P <- diff(diag(q-length(ind)),difference=1)
    object$P[[1]] <- P
    object$S[[1]] <- crossprod(P)
  }
  object$p.ident <- rep(TRUE,q-length(ind)) ## p.ident is an indicator of which coefficients must be positive (exponentiated)
  object$n.zero <- ind
  object$rank <- ncol(object$X)-1  # penalty rank
  object$null.space.dim <- 2 ## m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  object$knots <- xk;
  object$m <- m;
  object$df<-ncol(object$X)     # maximum DoF 

  ## get model matrix for 1st and 2nd derivatives of the smooth...
  h <- xk[q]-xk[q-1] ##(max(x)-min(x))/(q-m-1) ## distance between two adjacent knots (both expressions give the same distance)
  ## Xdf1, Xdf2 need to be checked...
  object$Xdf1 <- splineDesign(xk,x,ord=m+1)[,(m+2):(q-1)]/h ## ord is by one less for the 1st derivative
  object$Xdf2 <- splineDesign(xk,x,ord=m)[,(m+2):(q-2)]/h^2 ## ord is by two less for the 2nd derivative

  class(object)<-"miso.smooth"  # Give object a class
  object
}


## Prediction matrix for the `miso` smooth class... 
Predict.matrix.miso.smooth<-function(object,data)
## prediction method function for the `miso' smooth class
{ m <- object$m # spline order, m+1=3 default for cubic spline
  q <- object$bs.dim ## object$df +1 ## ?????
 # Sig <- matrix(as.numeric(rep(1:q,q)>=rep(1:q,each=q)),q,q) ## coef summation matrix
  Sig <- matrix(1,q,q)  ## coef summation matrix
  Sig[upper.tri(Sig)] <-0
  ind <- object$n.zero ## 1:(m+1)
  Sig[,ind] <- 0; Sig[ind,] <- 0
  ## find spline basis inner knot range...
  ll <- object$knots[m+2];ul <- object$knots[length(object$knots)-m-1]
  x <- data[[object$term]]
  n <- length(x)
  ind <- x<=ul & x>=ll ## data in range
  if (sum(ind)==n) { ## all in range
     X <- spline.des(object$knots,x,m+2)$design
     X <- X%*%Sig 
  } else { ## some extrapolation needed 
     ## matrix mapping coefs to value and slope at end points...
     D <- spline.des(object$knots,c(ll,ll,ul,ul),m+2,c(0,1,0,1))$design
     X <- matrix(0,n,ncol(D)) ## full predict matrix
     if (sum(ind)> 0)  X[ind,] <- spline.des(object$knots,x[ind],m+2)$design ## interior rows
     ## Now add rows for linear extrapolation...
     ind <- x < ll 
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ll)%*%D[1:2,]
     ind <- x > ul
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ul)%*%D[3:4,]
     X <- X%*%Sig 
  }
  X
}


##############################################################################################
### Adding Monotone increasing SCOP-spline construction with an 'finish at zero' constraint,
### a SCOP-spline additionally constrained to be zero at the end, at the right-end point of 
### the covariate range...
#############################################################################################

smooth.construct.mifo.smooth.spec<- function(object, data, knots)
## construction of the monotone increasing smooth with a 'finish-at-zero' constraint;
## achieved simply by setting the last (m+1) spline coefficients to zero...
{ 
  if (!is.null(object$point.con[[1]])) ## a point constraint is supplied?
    stop(paste("'mifo' smooth works only with a 'finish at zero' constraint; 'pc' argument does not work here"))
  m <- object$p.order[1]
  if (is.na(m)) m <- 2 ## default for cubic spline
  if (m<1) stop("silly m supplied")
  if (object$bs.dim<0) object$bs.dim <- 10 ## default
  q <- object$bs.dim 
  nk <- q+m+2 ## number of knots
  if (nk<=0) stop("k too small for m")
  x <- data[[object$term]]  ## the data
  xk <- knots[[object$term]] ## will be NULL if none supplied
  if (is.null(xk)) { # space knots evenly through data
       n<-length(x)
       xk<-rep(0,q+m+2)
       xk[(m+2):(q+1)]<-seq(min(x),max(x),length=q-m)
       for (i in 1:(m+1)) {xk[i]<-xk[m+2]-(m+2-i)*(xk[m+3]-xk[m+2])}
       for (i in (q+2):(q+m+2)) {xk[i]<-xk[q+1]+(i-q-1)*(xk[m+3]-xk[m+2])}
     }
  ## move m knots, xk[q-m+1],...,xk[q] (the m pre-last inner knots), to the last inner knot xk[q+1],
  ## join these m knots with the constrained xk[q+1] knot to avoid a plateau end...
  xk[(q-m+1):q] <- xk[q+1]
  if (length(xk)!=nk) # right number of knots?
   stop(paste("there should be ", nk," supplied knots"))
  ##  get model matrix...
  X1 <- splineDesign(xk,x,ord=m+2)
  ## get unconstrained matrix Sigma and remove the last m+1 columns and rows...
  Sig <- matrix(1,q,q)  ## coef summation matrix
  Sig[upper.tri(Sig)] <-0
  ind <- (q-m):q
 ## Sig[,ind] <- 0; Sig[ind,] <- 0
  Sig <- Sig[-ind,-ind]
  X <- X1[,-ind]%*%Sig # drop (m+1) end terms, model submatrix for the scop-term
  object$X <- X # the final model matrix
  object$P <- list()
  object$S <- list()
  object$Sigma <- Sig

  if (!object$fixed) # create the penalty matrix
  { P <- diff(diag(q-length(ind)),difference=1)
    object$P[[1]] <- P
    object$S[[1]] <- crossprod(P)
  }
  object$p.ident <- c(FALSE,rep(TRUE,q-length(ind)-1)) ## p.ident is an indicator of which coefficients must be positive (exponentiated)
  object$n.zero <- ind
  object$rank <- ncol(object$X)-1  # penalty rank
  object$null.space.dim <- 2 ##m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  object$knots <- xk;
  object$m <- m;
  object$df<-ncol(object$X)     # maximum DoF 

  ## get model matrix for 1st and 2nd derivatives of the smooth...
  h <- (max(x)-min(x))/(q-m-1) ## distance between two adjacent knots
  ## Xdf1, Xdf2 need to be checked...
  object$Xdf1 <- splineDesign(xk,x,ord=m+1)[,1:(q-1-(m+1))]/h ## ord is by one less for the 1st derivative
  object$Xdf2 <- splineDesign(xk,x,ord=m)[,1:(q-2-(m+1))]/h^2 ## ord is by two less for the 2nd derivative

  class(object)<-"mifo.smooth"  # Give object a class
  object
}


## Prediction matrix for the `mifo` smooth class... 
Predict.matrix.mifo.smooth<-function(object,data)
## prediction method function for the `mifo' smooth class
{ m <- object$m # spline order, m+1=3 default for cubic spline
  q <- object$bs.dim ## object$df +1 ## ?????
  Sig <- matrix(1,q,q)  ## coef summation matrix
  Sig[upper.tri(Sig)] <-0
  ind <- object$n.zero 
  Sig[,ind] <- 0; Sig[ind,] <- 0
  ## find spline basis inner knot range...
  ll <- object$knots[m+2];ul <- object$knots[length(object$knots)-m-1]
  x <- data[[object$term]]
  n <- length(x)
  ind <- x<=ul & x>=ll ## data in range
  if (sum(ind)==n) { ## all in range
     X <- spline.des(object$knots,x,m+2)$design
     X <- X%*%Sig 
  } else { ## some extrapolation needed 
     ## matrix mapping coefs to value and slope at end points...
     D <- spline.des(object$knots,c(ll,ll,ul,ul),m+2,c(0,1,0,1))$design
     X <- matrix(0,n,ncol(D)) ## full predict matrix
     if (sum(ind)> 0)  X[ind,] <- spline.des(object$knots,x[ind],m+2)$design ## interior rows
     ## Now add rows for linear extrapolation...
     ind <- x < ll 
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ll)%*%D[1:2,]
     ind <- x > ul
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ul)%*%D[3:4,]
     X <- X%*%Sig 
  }
  X
}



########################################################
### Adding Monotone decreasing SCOP-spline construction 
########################################################

smooth.construct.mpd.smooth.spec<- function(object, data, knots)
## construction of the monotone decreasing smooth
{ 
  # require(splines)
  m <- object$p.order[1]
  if (is.na(m)) m <- 2 ## default for cubic splines
  if (m<1) stop("silly m supplied")
  if (object$bs.dim<0) object$bs.dim <- 10 ## default
  q <- object$bs.dim 
  nk <- q+m+2 ## number of knots
  if (nk<=0) stop("k too small for m")
  x <- data[[object$term]]  ## the data
  xk <- knots[[object$term]] ## will be NULL if none supplied
  if (is.null(xk)) # space knots through data
     { n <- length(x)
       xk <- rep(0,q+m+2)
       xk[(m+2):(q+1)] <- seq(min(x),max(x),length=q-m)
       for (i in 1:(m+1)) {xk[i] <- xk[m+2]-(m+2-i)*(xk[m+3]-xk[m+2])}
       for (i in (q+2):(q+m+2)) {xk[i] <- xk[q+1]+(i-q-1)*(xk[m+3]-xk[m+2])}
  }
  if (length(xk)!=nk) # right number of knots?
  stop(paste("there should be ",nk," supplied knots"))
  #  get model matrix-------------
  X1 <- splineDesign(xk,x,ord=m+2)
  # get matrix Sigma and remove the first column for the monotone smooth (column and row for Sig)
 # Sig <- matrix(as.numeric(rep(1:q,q)>=rep(1:q,each=q)),q,q) ## coef summation matrix
  Sig <- matrix(-1,q,q)  ## coef summation matrix
  Sig[upper.tri(Sig)] <-0
  Sig[,1] <- -Sig[,1] ## monotone decrease case
  X <- X1[,-1]%*%Sig[-1,-1] # drop intercept term
  
  ## applying sum-to-zero (centering) constraint...
  cmx <- colMeans(X)
  X <- sweep(X,2,cmx) ## subtract cmx from columns 
  object$X <- X # the final model matrix
  object$cmX <- cmx  

  object$Sigma <- Sig[-1,-1]
  object$P <- list()
  object$S <- list()

  if (!object$fixed) # create the penalty matrix
  { P <- diff(diag(q-1),difference=1)
    object$P[[1]] <- P
    object$S[[1]] <- crossprod(P)
  }
  object$p.ident <- rep(TRUE,q-1)  ## p.ident is an indicator of which coefficients must be positive (exponentiated) 
  object$rank <- ncol(object$X)-1  # penalty rank
  object$null.space.dim <-2 ## m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  
  ## get model matrix for 1st and 2nd derivatives of the smooth...
  h <- (max(x)-min(x))/(q-m-1) ## distance between two adjacent knots
  object$Xdf1 <- splineDesign(xk,x,ord=m+1)[,2:(q-1)]/h ## ord is by one less for the 1st derivative
  object$Xdf2 <- splineDesign(xk,x,ord=m)[,2:(q-2)]/h^2 ## ord is by two less for the 2nd derivative

  object$knots <- xk;
  object$m <- m;
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  class(object)<-"mpd.smooth"  # Give object a class
  object
}


## Prediction matrix for the `mpd` smooth class 

Predict.matrix.mpd.smooth<-function(object,data)
## prediction method function for the `mpd' smooth class
{ m <- object$m+1; # spline order, m+1=3 default for cubic spline
  q <- object$df +1
  # elements of matrix Sigma for decreasing smooth...
 # Sig <- matrix(as.numeric(rep(1:q,q)>=rep(1:q,each=q)),q,q) ## coef summation matrix
  Sig <- matrix(-1,q,q)  ## coef summation matrix
  Sig[upper.tri(Sig)] <-0
  Sig[,1] <- -Sig[,1] ## monotone decrease case
  ## find spline basis inner knot range...
  ll <- object$knots[m+1];ul <- object$knots[length(object$knots)-m]
  m <- m + 1
  x <- data[[object$term]]
  n <- length(x)
  ind <- x<=ul & x>=ll ## data in range
  if (sum(ind)==n) { ## all in range
     X <- spline.des(object$knots,x,m)$design
     X <- X%*%Sig 
  } else { ## some extrapolation needed 
     ## matrix mapping coefs to value and slope at end points...
     D <- spline.des(object$knots,c(ll,ll,ul,ul),m,c(0,1,0,1))$design
     X <- matrix(0,n,ncol(D)) ## full predict matrix
     if (sum(ind)> 0)  X[ind,] <- spline.des(object$knots,x[ind],m)$design ## interior rows
     ## Now add rows for linear extrapolation...
     ind <- x < ll 
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ll)%*%D[1:2,]
     ind <- x > ul
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ul)%*%D[3:4,]
     X <- X%*%Sig 
    ## X <- sweep(X,2,c(0,object$cmX)) 
  }
  X <- sweep(X,2,c(0,object$cmX))
  X
}

#######################################################
### Adding Monotone decreasing SCOP-spline construction without applying identifiability constraints
### to be used with numeric 'by' variable...
########################################################

smooth.construct.mpdBy.smooth.spec<- function(object, data, knots)
## construction of the monotone decreasing smooth
{ 
  # require(splines)
  m <- object$p.order[1]
  if (is.na(m)) m <- 2 ## default for cubic splines
  if (m<1) stop("silly m supplied")
  if (object$bs.dim<0) object$bs.dim <- 10 ## default
  q <- object$bs.dim 
  nk <- q+m+2 ## number of knots
  if (nk<=0) stop("k too small for m")
  x <- data[[object$term]]  ## the data
  xk <- knots[[object$term]] ## will be NULL if none supplied
  if (is.null(xk)) # space knots through data
     { n <- length(x)
       xk <- rep(0,q+m+2)
       xk[(m+2):(q+1)] <- seq(min(x),max(x),length=q-m)
       for (i in 1:(m+1)) {xk[i] <- xk[m+2]-(m+2-i)*(xk[m+3]-xk[m+2])}
       for (i in (q+2):(q+m+2)) {xk[i] <- xk[q+1]+(i-q-1)*(xk[m+3]-xk[m+2])}
  }
  if (length(xk)!=nk) # right number of knots?
  stop(paste("there should be ",nk," supplied knots"))
  #  get model matrix-------------
  X1 <- splineDesign(xk,x,ord=m+2)
  # get matrix Sigma and remove the first column for the monotone smooth (column and row for Sig)
 # Sig <- matrix(as.numeric(rep(1:q,q)>=rep(1:q,each=q)),q,q) ## coef summation matrix
  Sig <- matrix(-1,q,q)  ## coef summation matrix
  Sig[upper.tri(Sig)] <-0
  Sig[,1] <- -Sig[,1] ## monotone decrease case

  X <- X1%*%Sig # no drop intercept term
  object$X<-X # the final model matrix
  object$Sigma <- Sig
  object$P <- list()
  object$S <- list()

  if (!object$fixed) # create the penalty matrix
  { P <- diff(diag(q-1),difference=1)
    P <- rbind(rep(0,q-1),P) ## adding 1st row of zeros
    P <- cbind(rep(0,q-1),P) ## adding first column of zeros
    object$P[[1]] <- P
    object$S[[1]] <- crossprod(P)
  }
  object$p.ident <- c(FALSE,rep(TRUE,q-1))  ## p.ident is an indicator of which coefficients must be positive (exponentiated) 
  object$rank <- ncol(object$X)  # penalty rank
  object$null.space.dim <- 2 ## m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  
  ## get model matrix for 1st and 2nd derivatives of the smooth...
  h <- (max(x)-min(x))/(q-m-1) ## distance between two adjacent knots
  object$Xdf1 <- splineDesign(xk,x,ord=m+1)[,1:(q-1)]/h ## ord is by one less for the 1st derivative
  object$Xdf2 <- splineDesign(xk,x,ord=m)[,1:(q-2)]/h^2 ## ord is by two less for the 2nd derivative

  object$knots <- xk;
  object$m <- m;
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  class(object)<-"mpdBy.smooth"  # Give object a class
  object
}


## Prediction matrix for the `mpdBy` smooth class 

Predict.matrix.mpdBy.smooth<-function(object,data)
## prediction method function for the `mpdBy' smooth class
{ m <- object$m+1; # spline order, m+1=3 default for cubic spline
  q <- object$df 
  # elements of matrix Sigma for decreasing smooth...
 # Sig <- matrix(as.numeric(rep(1:q,q)>=rep(1:q,each=q)),q,q) ## coef summation matrix
  Sig <- matrix(-1,q,q)  ## coef summation matrix
  Sig[upper.tri(Sig)] <-0
  Sig[,1] <- -Sig[,1] ## monotone decrease case
  ## find spline basis inner knot range...
  ll <- object$knots[m+1];ul <- object$knots[length(object$knots)-m]
  m <- m + 1
  x <- data[[object$term]]
  n <- length(x)
  ind <- x<=ul & x>=ll ## data in range
  if (sum(ind)==n) { ## all in range
     X <- spline.des(object$knots,x,m)$design
     X <- X%*%Sig 
  } else { ## some extrapolation needed 
     ## matrix mapping coefs to value and slope at end points...
     D <- spline.des(object$knots,c(ll,ll,ul,ul),m,c(0,1,0,1))$design
     X <- matrix(0,n,ncol(D)) ## full predict matrix
     if (sum(ind)> 0)  X[ind,] <- spline.des(object$knots,x[ind],m)$design ## interior rows
     ## Now add rows for linear extrapolation...
     ind <- x < ll 
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ll)%*%D[1:2,]
     ind <- x > ul
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ul)%*%D[3:4,]
     X <- X%*%Sig 
  }
  X
}


##############################################################
### Smooth constructor for the mixed constrainted smooths ......
##############################################################

###############################################################
### Adding Monotone decreasing & concave P-spline construction 
###############################################################

smooth.construct.mdcv.smooth.spec<- function(object, data, knots)
## construction of the monotone decreasing and concave smooth
{ 
  # require(splines)
  m <- object$p.order[1]
  if (is.na(m)) m <- 2 ## default for cubis spline
  if (m<1) stop("silly m supplied")
  if (object$bs.dim<0) object$bs.dim <- 10 ## default
  q <- object$bs.dim 
  nk <- q+m+2 ## number of knots
  if (nk<=0) stop("k too small for m")
  x <- data[[object$term]]  ## the data
  xk <- knots[[object$term]] ## will be NULL if none supplied
  if (is.null(xk)) # space knots through data
     { n<-length(x)
       xk<-rep(0,q+m+2)
       xk[(m+2):(q+1)]<-seq(min(x),max(x),length=q-m)
       for (i in 1:(m+1)) {xk[i]<-xk[m+2]-(m+2-i)*(xk[m+3]-xk[m+2])}
       for (i in (q+2):(q+m+2)) {xk[i]<-xk[q+1]+(i-q-1)*(xk[m+3]-xk[m+2])}
     }
  if (length(xk)!=nk) # right number of knots?
  stop(paste("there should be ",nk," supplied knots"))
  #  get model matrix-------------
  X1 <- splineDesign(xk,x,ord=m+2)
  # use matrix Sigma and remove the first column for the monotone smooth
  Sig <- matrix(0,(q-1),(q-1))   # Define Matrix Sigma
     # for monotone decreasing & concave smooth
  for (i in 1:(q-1)) Sig[i:(q-1),i]<--c(1:(q-i))
  X <- X1[,2:q]%*%Sig # model submatrix for the constrained term

  ## applying sum-to-zero (centering) constraint...
  cmx <- colMeans(X)
  X <- sweep(X,2,cmx) ## subtract cmx from columns 
  object$X <- X # the final model matrix
  object$cmX <- cmx

  object$Sigma <- Sig
  object$P <- list()
  object$S <- list()

  if (!object$fixed) # create the penalty matrix
  { P <- diff(diag(q-2),difference=1)
    object$P[[1]] <- P
    S <- matrix(0,q-1,q-1)
    S[2:(q-1),2:(q-1)] <- crossprod(P)
    object$S[[1]] <- S
  }
  object$p.ident <- rep(TRUE,q-1)  ## p.ident is an indicator of which coefficients must be positive (exponentiated)  
  object$rank <- ncol(object$X)-1  # penalty rank
  object$null.space.dim <- 2 ##m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  
  ## get model matrix for 1st and 2nd derivatives of the smooth...
  h <- (max(x)-min(x))/(q-m-1) ## distance between two adjacent knots
  object$Xdf1 <- splineDesign(xk,x,ord=m+1)[,2:(q-1)]/h ## ord is by one less for the 1st derivative
  object$Xdf2 <- splineDesign(xk,x,ord=m)[,2:(q-2)]/h^2 ## ord is by two less for the 2nd derivative

  object$knots <- xk;
  object$m <- m;
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  class(object)<-"mdcv.smooth"  # Give object a class
  object
}


## Prediction matrix for the `mdcv` smooth class.... 

Predict.matrix.mdcv.smooth<-function(object,data)
## prediction method function for the `mdcv' smooth class
{ m <- object$m+1; # spline order, m+1=3 default for cubic spline
  q <- object$df +1
  Sig <- matrix(0,q,q)   
  Sig[,1] <- rep(1,q)
  for (i in 2:q)  Sig[i:q,i] <- -c(1:(q-i+1))
  ## find spline basis inner knot range...
  ll <- object$knots[m+1];ul <- object$knots[length(object$knots)-m]
  m <- m + 1
  x <- data[[object$term]]
  n <- length(x)
  ind <- x<=ul & x>=ll ## data in range
  if (sum(ind)==n) { ## all in range
     X <- spline.des(object$knots,x,m)$design
     X <- X%*%Sig 
  } else { ## some extrapolation needed 
     ## matrix mapping coefs to value and slope at end points...
     D <- spline.des(object$knots,c(ll,ll,ul,ul),m,c(0,1,0,1))$design
     X <- matrix(0,n,ncol(D)) ## full predict matrix
     if (sum(ind)> 0)  X[ind,] <- spline.des(object$knots,x[ind],m)$design ## interior rows
     ## Now add rows for linear extrapolation...
     ind <- x < ll 
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ll)%*%D[1:2,]
     ind <- x > ul
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ul)%*%D[3:4,]
     X <- X%*%Sig 
    ## X <- sweep(X,2,c(0,object$cmX))
  }
  X <- sweep(X,2,c(0,object$cmX))
  X
}


###########################################################
### Adding decreasing & concave SCOP-spline construction without identifiability constraints
### to be used with numeric 'by' variable and linear functional terms...
##########################################################

smooth.construct.mdcvBy.smooth.spec<- function(object, data, knots)
## construction of the monotone decreasing and concave smooth
{ 
  m <- object$p.order[1]
  if (is.na(m)) m <- 2 ## default for cubis spline
  if (m<1) stop("silly m supplied")
  if (object$bs.dim<0) object$bs.dim <- 10 ## default
  q <- object$bs.dim 
  nk <- q+m+2 ## number of knots
  if (nk<=0) stop("k too small for m")
  x <- data[[object$term]]  ## the data
  xk <- knots[[object$term]] ## will be NULL if none supplied
  if (is.null(xk)) # space knots through data
     { n<-length(x)
       xk<-rep(0,q+m+2)
       xk[(m+2):(q+1)]<-seq(min(x),max(x),length=q-m)
       for (i in 1:(m+1)) {xk[i]<-xk[m+2]-(m+2-i)*(xk[m+3]-xk[m+2])}
       for (i in (q+2):(q+m+2)) {xk[i]<-xk[q+1]+(i-q-1)*(xk[m+3]-xk[m+2])}
     }
  if (length(xk)!=nk) # right number of knots?
  stop(paste("there should be ",nk," supplied knots"))
  #  get model matrix-------------
  X <- splineDesign(xk,x,ord=m+2)
  # Define matrix Sigma for monotone decreasing & concave smooth
  Sig <- matrix(0,q,q)   
  Sig[,1] <- rep(1,q)
  for (i in 2:q)  Sig[i:q,i] <- -c(1:(q-i+1))

  X <- X%*%Sig # model matrix for the constrained term
  object$X <- X 
  object$Sigma <- Sig
  object$P <- list()
  object$S <- list()

  if (!object$fixed) # create the penalty matrix
  { P <- diff(diag(q-2),difference=1)
    P <- rbind(matrix(0,2,q-2),P) ## adding first two rows of zeros
    P <- cbind(matrix(0,q-1,2),P) ## adding first two columns of zeros
    object$P[[1]] <- P
    object$S[[1]] <- crossprod(P)
  }
  object$p.ident <- c(FALSE,rep(TRUE,q-1))  ## p.ident is an indicator of which coefficients must be positive (exponentiated)  
  object$rank <- ncol(object$X)  # penalty rank
  object$null.space.dim <- 2 ## m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  
  ## get model matrix for 1st and 2nd derivatives of the smooth...
  h <- (max(x)-min(x))/(q-m-1) ## distance between two adjacent knots
  object$Xdf1 <- splineDesign(xk,x,ord=m+1)[,1:(q-1)]/h ## ord is by one less for the 1st derivative
  object$Xdf2 <- splineDesign(xk,x,ord=m)[,1:(q-2)]/h^2 ## ord is by two less for the 2nd derivative

  object$knots <- xk;
  object$m <- m;
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  class(object)<-"mdcvBy.smooth"  # Give object a class
  object
}


## Prediction matrix for the `mdcvBy` smooth class.... 

Predict.matrix.mdcvBy.smooth<-function(object,data)
## prediction method function for the `mdcvBy' smooth class
{ m <- object$m+1; # spline order, m+1=3 default for cubic spline
  q <- object$df 
  Sig <- matrix(0,q,q)   
  Sig[,1] <- rep(1,q)
  for (i in 2:q)  Sig[i:q,i] <- -c(1:(q-i+1))
  ## find spline basis inner knot range...
  ll <- object$knots[m+1];ul <- object$knots[length(object$knots)-m]
  m <- m + 1
  x <- data[[object$term]]
  n <- length(x)
  ind <- x<=ul & x>=ll ## data in range
  if (sum(ind)==n) { ## all in range
     X <- spline.des(object$knots,x,m)$design
     X <- X%*%Sig 
  } else { ## some extrapolation needed 
     ## matrix mapping coefs to value and slope at end points...
     D <- spline.des(object$knots,c(ll,ll,ul,ul),m,c(0,1,0,1))$design
     X <- matrix(0,n,ncol(D)) ## full predict matrix
     if (sum(ind)> 0)  X[ind,] <- spline.des(object$knots,x[ind],m)$design ## interior rows
     ## Now add rows for linear extrapolation...
     ind <- x < ll 
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ll)%*%D[1:2,]
     ind <- x > ul
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ul)%*%D[3:4,]
     X <- X%*%Sig 
  }
  X
}


###############################################################
### Adding decreasing & convex SCOP-spline construction 
################################################################

smooth.construct.mdcx.smooth.spec<- function(object, data, knots)
##  the constructor for the monotone decreasing and convex smooth
{ 
  #require(splines)
  m <- object$p.order[1]
  if (is.na(m)) m <- 2 ## default 
  if (m<1) stop("silly m supplied")
  if (object$bs.dim<0) object$bs.dim <- 10 ## default
  q <- object$bs.dim 
  nk <- q+m+2 ## number of knots
  if (nk<=0) stop("k too small for m")
  x <- data[[object$term]]  ## the data
  xk <- knots[[object$term]] ## will be NULL if none supplied
  if (is.null(xk)) # space knots through data
  { n<-length(x)
    xk<-rep(0,q+m+2)
    xk[(m+2):(q+1)]<-seq(min(x),max(x),length=q-m)
    for (i in 1:(m+1)) {xk[i]<-xk[m+2]-(m+2-i)*(xk[m+3]-xk[m+2])}
    for (i in (q+2):(q+m+2)) {xk[i]<-xk[q+1]+(i-q-1)*(xk[m+3]-xk[m+2])}
  }
  if (length(xk)!=nk) # right number of knots?
  stop(paste("there should be ",nk," supplied knots"))
  #  get model matrix-------------
  X1 <- splineDesign(xk,x,ord=m+2)
  # use matrix Sigma and remove the first column for the monotone smooth
  Sig <- matrix(0,(q-1),(q-1))   # Define Matrix Sigma
     # for monotone decreasing & convex smooth
  Sig[1,] <- -rep(1,q-1)
  for (i in 2:(q-1)) {
         Sig[i,1:(q-i)] <- -i;
         Sig[i,(q-i+1):(q-1)] <- -c((i-1):1)
  }
  X <- X1[,2:q]%*%Sig # model submatrix for the constrained term

  ## applying sum-to-zero (centering) constraint...
  cmx <- colMeans(X)
  X <- sweep(X,2,cmx) ## subtract cmx from columns 
  object$X <- X # the final model matrix
  object$cmX <- cmx

  object$Sigma <- Sig
  object$P <- list()
  object$S <- list()

  if (!object$fixed) # create the penalty matrix
  { P <- diff(diag(q-2),difference=1)
    object$P[[1]] <- P
    S <- matrix(0,q-1,q-1)
    S[2:(q-1),2:(q-1)] <- crossprod(P)
    object$S[[1]] <- S
  }
  object$p.ident <- rep(TRUE,q-1)  ## p.ident is an indicator of which coefficients must be positive (exponentiated)  
  object$rank <- ncol(object$X)-1  # penalty rank
  object$null.space.dim <- 2 ## m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  
  ## get model matrix for 1st and 2nd derivatives of the smooth...
  h <- (max(x)-min(x))/(q-m-1) ## distance between two adjacent knots
  object$Xdf1 <- splineDesign(xk,x,ord=m+1)[,2:(q-1)]/h ## ord is by one less for the 1st derivative
  object$Xdf2 <- splineDesign(xk,x,ord=m)[,2:(q-2)]/h^2 ## ord is by two less for the 2nd derivative

  object$knots <- xk;
  object$m <- m;
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  class(object)<-"mdcx.smooth"  # Give object a class
  object
}


## Prediction matrix for the `mdcx` smooth class.... 
Predict.matrix.mdcx.smooth<-function(object,data)
## the prediction method for the `mdcx' smooth class
{ x <- data[[object$term]]
  m <- object$m+1; # spline order, m+1=3 default for cubic spline
  q <- object$df +1
  Sig <- matrix(0,q,q)   
  Sig[,1] <- rep(1,q)
  Sig1 <- matrix(0,(q-1),(q-1))   
  Sig1[1,] <- -rep(1,q-1)
  for (i in 2:(q-1)) {
         Sig1[i,1:(q-i)]<--i;
         Sig1[i,(q-i+1):(q-1)]<--c((i-1):1)
  }
  Sig [2:q,2:q] <- Sig1

  ## find spline basis inner knot range...
  ll <- object$knots[m+1];ul <- object$knots[length(object$knots)-m]
  m <- m + 1
  n <- length(x)
  ind <- x<=ul & x>=ll ## data in range
  if (sum(ind)==n) { ## all in range
     X <- spline.des(object$knots,x,m)$design
     X <- X%*%Sig 
  } else { ## some extrapolation needed 
     ## matrix mapping coefs to value and slope at end points...
     D <- spline.des(object$knots,c(ll,ll,ul,ul),m,c(0,1,0,1))$design
     X <- matrix(0,n,ncol(D)) ## full predict matrix
     if (sum(ind)> 0)  X[ind,] <- spline.des(object$knots,x[ind],m)$design ## interior rows
     ## Now add rows for linear extrapolation...
     ind <- x < ll 
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ll)%*%D[1:2,]
     ind <- x > ul
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ul)%*%D[3:4,]
     X <- X%*%Sig 
   ##  X <- sweep(X,2,c(0,object$cmX))
  }
  X <- sweep(X,2,c(0,object$cmX))
  X
}

###########################################################
### Adding decreasing & convex SCOP-spline construction without identifiability constraints
### to be used with numeric 'by' variable and linear functional terms...
##########################################################

smooth.construct.mdcxBy.smooth.spec<- function(object, data, knots)
##  the constructor for the monotone decreasing and convex smooth
{ 
  m <- object$p.order[1]
  if (is.na(m)) m <- 2 ## default 
  if (m<1) stop("silly m supplied")
  if (object$bs.dim<0) object$bs.dim <- 10 ## default
  q <- object$bs.dim 
  nk <- q+m+2 ## number of knots
  if (nk<=0) stop("k too small for m")
  x <- data[[object$term]]  ## the data
  xk <- knots[[object$term]] ## will be NULL if none supplied
  if (is.null(xk)) # space knots through data
  { n<-length(x)
    xk<-rep(0,q+m+2)
    xk[(m+2):(q+1)]<-seq(min(x),max(x),length=q-m)
    for (i in 1:(m+1)) {xk[i]<-xk[m+2]-(m+2-i)*(xk[m+3]-xk[m+2])}
    for (i in (q+2):(q+m+2)) {xk[i]<-xk[q+1]+(i-q-1)*(xk[m+3]-xk[m+2])}
  }
  if (length(xk)!=nk) # right number of knots?
  stop(paste("there should be ",nk," supplied knots"))
  #  get model matrix-------------
  X <- splineDesign(xk,x,ord=m+2)
  # define matrix Sigma for monotone decreasing & convex smooth
  Sig <- matrix(0,q,q)   
  Sig[,1] <- rep(1,q)
  Sig1 <- matrix(0,(q-1),(q-1))   
  Sig1[1,] <- -rep(1,q-1)
  for (i in 2:(q-1)) {
         Sig1[i,1:(q-i)]<--i;
         Sig1[i,(q-i+1):(q-1)]<--c((i-1):1)
  }
  Sig [2:q,2:q] <- Sig1

  X <- X%*%Sig # model matrix for the constrained term
  object$X<-X 
  object$Sigma <- Sig
  object$P <- list()
  object$S <- list()

  if (!object$fixed) # create the penalty matrix
  { P <- diff(diag(q-2),difference=1)
    P <- rbind(matrix(0,2,q-2),P) ## adding first two rows of zeros
    P <- cbind(matrix(0,q-1,2),P) ## adding first two columns of zeros
    object$P[[1]] <- P
    object$S[[1]] <- crossprod(P)
  }
  object$p.ident <- c(FALSE,rep(TRUE,q-1))  ## p.ident is an indicator of which coefficients must be positive (exponentiated)  
  object$rank <- ncol(object$X)  # penalty rank
  object$null.space.dim <-2 ## m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  
  ## get model matrix for 1st and 2nd derivatives of the smooth...
  h <- (max(x)-min(x))/(q-m-1) ## distance between two adjacent knots
  object$Xdf1 <- splineDesign(xk,x,ord=m+1)[,1:(q-1)]/h ## ord is by one less for the 1st derivative
  object$Xdf2 <- splineDesign(xk,x,ord=m)[,1:(q-2)]/h^2 ## ord is by two less for the 2nd derivative

  object$knots <- xk;
  object$m <- m;
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  class(object)<-"mdcxBy.smooth"  # Give object a class
  object
}


## Prediction matrix for the `mdcxBy` smooth class.... 
Predict.matrix.mdcxBy.smooth<-function(object,data)
## the prediction method for the `mdcxBy' smooth class
{ x <- data[[object$term]]
  m <- object$m+1; # spline order, m+1=3 default for cubic spline
  q <- object$df
  Sig <- matrix(0,q,q)   
  Sig[,1] <- rep(1,q)
  Sig1 <- matrix(0,(q-1),(q-1))   
  Sig1[1,] <- -rep(1,q-1)
  for (i in 2:(q-1)) {
         Sig1[i,1:(q-i)]<--i;
         Sig1[i,(q-i+1):(q-1)]<--c((i-1):1)
  }
  Sig [2:q,2:q] <- Sig1

  ## find spline basis inner knot range...
  ll <- object$knots[m+1];ul <- object$knots[length(object$knots)-m]
  m <- m + 1
  n <- length(x)
  ind <- x<=ul & x>=ll ## data in range
  if (sum(ind)==n) { ## all in range
     X <- spline.des(object$knots,x,m)$design
     X <- X%*%Sig 
  } else { ## some extrapolation needed 
     ## matrix mapping coefs to value and slope at end points...
     D <- spline.des(object$knots,c(ll,ll,ul,ul),m,c(0,1,0,1))$design
     X <- matrix(0,n,ncol(D)) ## full predict matrix
     if (sum(ind)> 0)  X[ind,] <- spline.des(object$knots,x[ind],m)$design ## interior rows
     ## Now add rows for linear extrapolation...
     ind <- x < ll 
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ll)%*%D[1:2,]
     ind <- x > ul
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ul)%*%D[3:4,]
     X <- X%*%Sig 
  }
  X
}



################################################################
### Adding monotone increasing & concave SCOP-spline construction ################################################################

smooth.construct.micv.smooth.spec<- function(object, data, knots)
## construction of the monotone increasing and concave smooth
{ 
 # require(splines)
  m <- object$p.order[1]
  if (is.na(m)) m <- 2 ## default 
  if (m<1) stop("silly m supplied")
  if (object$bs.dim<0) object$bs.dim <- 10 ## default
  q <- object$bs.dim 
  nk <- q+m+2 ## number of knots
  if (nk<=0) stop("k too small for m")
  x <- data[[object$term]]  ## the data
  xk <- knots[[object$term]] ## will be NULL if none supplied
  if (is.null(xk)) # space knots through data
  { n<-length(x)
    xk<-rep(0,q+m+2)
    xk[(m+2):(q+1)]<-seq(min(x),max(x),length=q-m)
    for (i in 1:(m+1)) {xk[i]<-xk[m+2]-(m+2-i)*(xk[m+3]-xk[m+2])}
    for (i in (q+2):(q+m+2)) {xk[i]<-xk[q+1]+(i-q-1)*(xk[m+3]-xk[m+2])}
  }
  if (length(xk)!=nk) # right number of knots?
  stop(paste("there should be ",nk," supplied knots"))
  #  get model matrix-------------
  X1 <- splineDesign(xk,x,ord=m+2)
  # use matrix Sigma and remove the first column for the monotone smooth
  Sig <- matrix(0,(q-1),(q-1))   # Define Matrix Sigma
     # for monotone increasing & concave smooth
  Sig[1,]<-rep(1,q-1)
  for (i in 2:(q-1)) {
       Sig[i,1:(q-i)]<-i;
       Sig[i,(q-i+1):(q-1)]<-c((i-1):1)
  }
  X <- X1[,2:q]%*%Sig # model submatrix for the constrained term

  ## applying sum-to-zero (centering) constraint...
  cmx <- colMeans(X)
  X <- sweep(X,2,cmx) ## subtract cmx from columns 
  object$X <- X # the final model matrix
  object$cmX <- cmx

  object$Sigma <- Sig
  object$P <- list()
  object$S <- list()

  if (!object$fixed) # create the penalty matrix
  { P <- diff(diag(q-2),difference=1)
    object$P[[1]] <- P
    S <- matrix(0,q-1,q-1)
    S[2:(q-1),2:(q-1)] <- crossprod(P)
    object$S[[1]] <- S
  }
  object$p.ident <- rep(TRUE,q-1)  ## p.ident is an indicator of which coefficients must be positive (exponentiated) 
  object$rank <- ncol(object$X)-1  # penalty rank
  object$null.space.dim <- 2 ##m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  
  ## get model matrix for 1st and 2nd derivatives of the smooth...
  h <- (max(x)-min(x))/(q-m-1) ## distance between two adjacent knots
  object$Xdf1 <- splineDesign(xk,x,ord=m+1)[,2:(q-1)]/h ## ord is by one less for the 1st derivative
  object$Xdf2 <- splineDesign(xk,x,ord=m)[,2:(q-2)]/h^2 ## ord is by two less for the 2nd derivative

  object$knots <- xk;
  object$m <- m;
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  class(object)<-"micv.smooth"  # Give object a class
  object
}


## Prediction matrix for the `micv` smooth class....

Predict.matrix.micv.smooth<-function(object,data)
## prediction method function for the `micv' smooth class
{ m <- object$m+1; # spline order, m+1=3 default for cubic spline
  q <- object$df +1
  Sig <- matrix(0,q,q)   
  Sig[,1] <- rep(1,q)
  Sig1 <- matrix(0,(q-1),(q-1))   
  Sig1[1,] <- rep(1,q-1)
  for (i in 2:(q-1)) {
       Sig1[i,1:(q-i)] <- i;
       Sig1[i,(q-i+1):(q-1)] <- c((i-1):1)
  }
  Sig [2:q,2:q] <- Sig1

  ## find spline basis inner knot range...
  ll <- object$knots[m+1];ul <- object$knots[length(object$knots)-m]
  m <- m + 1
  x <- data[[object$term]]
  n <- length(x)
  ind <- x<=ul & x>=ll ## data in range
  if (sum(ind)==n) { ## all in range
     X <- spline.des(object$knots,x,m)$design
     X <- X%*%Sig 
  } else { ## some extrapolation needed 
     ## matrix mapping coefs to value and slope at end points...
     D <- spline.des(object$knots,c(ll,ll,ul,ul),m,c(0,1,0,1))$design
     X <- matrix(0,n,ncol(D)) ## full predict matrix
     if (sum(ind)> 0)  X[ind,] <- spline.des(object$knots,x[ind],m)$design ## interior rows
     ## Now add rows for linear extrapolation...
     ind <- x < ll 
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ll)%*%D[1:2,]
     ind <- x > ul
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ul)%*%D[3:4,]
     X <- X%*%Sig 
    ## X <- sweep(X,2,c(0,object$cmX))
  }
  X <- sweep(X,2,c(0,object$cmX))
  X
}


###########################################################
### Adding increasing & concave SCOP-spline construction without identifiability constraints
### to be used with numeric 'by' variable and linear functional terms...
##########################################################

smooth.construct.micvBy.smooth.spec<- function(object, data, knots)
## construction of the monotone increasing and concave smooth
{ 
  m <- object$p.order[1]
  if (is.na(m)) m <- 2 ## default 
  if (m<1) stop("silly m supplied")
  if (object$bs.dim<0) object$bs.dim <- 10 ## default
  q <- object$bs.dim 
  nk <- q+m+2 ## number of knots
  if (nk<=0) stop("k too small for m")
  x <- data[[object$term]]  ## the data
  xk <- knots[[object$term]] ## will be NULL if none supplied
  if (is.null(xk)) # space knots through data
  { n<-length(x)
    xk<-rep(0,q+m+2)
    xk[(m+2):(q+1)]<-seq(min(x),max(x),length=q-m)
    for (i in 1:(m+1)) {xk[i]<-xk[m+2]-(m+2-i)*(xk[m+3]-xk[m+2])}
    for (i in (q+2):(q+m+2)) {xk[i]<-xk[q+1]+(i-q-1)*(xk[m+3]-xk[m+2])}
  }
  if (length(xk)!=nk) # right number of knots?
  stop(paste("there should be ",nk," supplied knots"))
  ##  get model matrix-------------
  X <- splineDesign(xk,x,ord=m+2)
  ## define matrix Sigma for monotone increasing & concave smooth
  Sig <- matrix(0,q,q)   
  Sig[,1] <- rep(1,q)
  Sig1 <- matrix(0,(q-1),(q-1))   
  Sig1[1,] <- rep(1,q-1)
  for (i in 2:(q-1)) {
       Sig1[i,1:(q-i)] <- i;
       Sig1[i,(q-i+1):(q-1)] <- c((i-1):1)
  }
  Sig [2:q,2:q] <- Sig1

  X <- X%*%Sig ## model matrix for the constrained term
  object$X <- X 
  object$Sigma <- Sig
  object$P <- list()
  object$S <- list()

  if (!object$fixed) # create the penalty matrix
  { P <- diff(diag(q-2),difference=1)
    P <- rbind(matrix(0,2,q-2),P) ## adding first two rows of zeros
    P <- cbind(matrix(0,q-1,2),P) ## adding first two columns of zeros
    object$P[[1]] <- P
    object$S[[1]] <- crossprod(P)
  }
  object$p.ident <- c(FALSE,rep(TRUE,q-1))  ## p.ident is an indicator of which coefficients must be positive (exponentiated) 
  object$rank <- ncol(object$X)  # penalty rank
  object$null.space.dim <- 2 ##m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  
  ## get model matrix for 1st and 2nd derivatives of the smooth...
  h <- (max(x)-min(x))/(q-m-1) ## distance between two adjacent knots
  object$Xdf1 <- splineDesign(xk,x,ord=m+1)[,1:(q-1)]/h ## ord is by one less for the 1st derivative
  object$Xdf2 <- splineDesign(xk,x,ord=m)[,1:(q-2)]/h^2 ## ord is by two less for the 2nd derivative

  object$knots <- xk;
  object$m <- m;
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  class(object)<-"micvBy.smooth"  # Give object a class
  object
}


## Prediction matrix for the `micvBy` smooth class....
Predict.matrix.micvBy.smooth<-function(object,data)
 ## prediction method function for the `micvBy' smooth class
{ m <- object$m+1; # spline order, m+1=3 default for cubic spline
  q <- object$df
  Sig <- matrix(0,q,q)   
  Sig[,1] <- rep(1,q)
  Sig1 <- matrix(0,(q-1),(q-1))   
  Sig1[1,] <- rep(1,q-1)
  for (i in 2:(q-1)) {
       Sig1[i,1:(q-i)] <- i;
       Sig1[i,(q-i+1):(q-1)] <- c((i-1):1)
  }
  Sig [2:q,2:q] <- Sig1

  ## find spline basis inner knot range...
  ll <- object$knots[m+1];ul <- object$knots[length(object$knots)-m]
  m <- m + 1
  x <- data[[object$term]]
  n <- length(x)
  ind <- x<=ul & x>=ll ## data in range
  if (sum(ind)==n) { ## all in range
     X <- spline.des(object$knots,x,m)$design
     X <- X%*%Sig 
  } else { ## some extrapolation needed 
     ## matrix mapping coefs to value and slope at end points...
     D <- spline.des(object$knots,c(ll,ll,ul,ul),m,c(0,1,0,1))$design
     X <- matrix(0,n,ncol(D)) ## full predict matrix
     if (sum(ind)> 0)  X[ind,] <- spline.des(object$knots,x[ind],m)$design ## interior rows
     ## Now add rows for linear extrapolation...
     ind <- x < ll 
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ll)%*%D[1:2,]
     ind <- x > ul
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ul)%*%D[3:4,]
     X <- X%*%Sig 
  }
  X
}


###############################################################
### Adding monotone increasing & convex SCOP-spline construction 
################################################################

smooth.construct.micx.smooth.spec<- function(object, data, knots)
## construction of the monotone increasing and convex smooth
{ 
  # require(splines)
  m <- object$p.order[1]
  if (is.na(m)) m <- 2 ## default 
  if (m < 1) stop("silly m supplied")
  if (object$bs.dim<0) object$bs.dim <- 10 ## default
  q <- object$bs.dim 
  nk <- q+m+2 ## number of knots
  if (nk<=0) stop("k too small for m")
  x <- data[[object$term]]  ## the data
  xk <- knots[[object$term]] ## will be NULL if none supplied
  if (is.null(xk)) # space knots through data
  { n<-length(x)
    xk<-rep(0,q+m+2)
    xk[(m+2):(q+1)]<-seq(min(x),max(x),length=q-m)
    for (i in 1:(m+1)) {xk[i]<-xk[m+2]-(m+2-i)*(xk[m+3]-xk[m+2])}
    for (i in (q+2):(q+m+2)) {xk[i]<-xk[q+1]+(i-q-1)*(xk[m+3]-xk[m+2])}
  }
  if (length(xk)!=nk) # right number of knots?
  stop(paste("there should be ",nk," supplied knots"))
  #  get model matrix-------------
  X1 <- splineDesign(xk,x,ord=m+2)
  # use matrix Sigma and remove the first column for the monotone smooth
  Sig <- matrix(0,(q-1),(q-1))   # Define Matrix Sigma
     # for monotone increasing & convex smooth
  for (i in 1:(q-1)) Sig[i:(q-1),i]<-c(1:(q-i))
  X <- X1[,2:q]%*%Sig # model submatrix for the constrained term

  ## applying sum-to-zero (centering) constraint...
  cmx <- colMeans(X)
  X <- sweep(X,2,cmx) ## subtract cmx from columns 
  object$X <- X # the final model matrix
  object$cmX <- cmx

  object$Sigma <- Sig
  object$P <- list()
  object$S <- list()

  if (!object$fixed) # create the penalty matrix
  { P <- diff(diag(q-2),difference=1)
    object$P[[1]] <- P
    S <- matrix(0,q-1,q-1)
    S[2:(q-1),2:(q-1)] <- crossprod(P)
    object$S[[1]] <- S
  }
  object$p.ident <- rep(TRUE,q-1) ## p.ident is an indicator of which coefficients must be positive (exponentiated) 
  object$rank <- ncol(object$X)-1  # penalty rank
  object$null.space.dim <- 2 ##m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
 
  ## get model matrix for 1st and 2nd derivatives of the smooth...
  h <- (max(x)-min(x))/(q-m-1) ## distance between two adjacent knots
  object$Xdf1 <- splineDesign(xk,x,ord=m+1)[,2:(q-1)]/h ## ord is by one less for the 1st derivative
  object$Xdf2 <- splineDesign(xk,x,ord=m)[,2:(q-2)]/h^2 ## ord is by two less for the 2nd derivative

  object$knots <- xk;
  object$m <- m;
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  class(object)<-"micx.smooth"  # Give object a class
  object
}

## Prediction matrix for the `micx` smooth class... 

Predict.matrix.micx.smooth<-function(object,data)
## prediction method function for the `micx` smooth class
{ m <- object$m+1; # spline order, m+1=3 default for cubic spline
  q <- object$df +1
  Sig <- matrix(0,q,q)   
  Sig[,1] <- rep(1,q)
  for (i in 2:q) Sig[i:q,i] <- c(1:(q-i+1))
  
  ## find spline basis inner knot range...
  ll <- object$knots[m+1];ul <- object$knots[length(object$knots)-m]
  m <- m + 1
  x <- data[[object$term]]
  n <- length(x)
  ind <- x<=ul & x>=ll ## data in range
  if (sum(ind)==n) { ## all in range
     X <- spline.des(object$knots,x,m)$design
     X <- X%*%Sig 
  } else { ## some extrapolation needed 
     ## matrix mapping coefs to value and slope at end points...
     D <- spline.des(object$knots,c(ll,ll,ul,ul),m,c(0,1,0,1))$design
     X <- matrix(0,n,ncol(D)) ## full predict matrix
     if (sum(ind)> 0)  X[ind,] <- spline.des(object$knots,x[ind],m)$design ## interior rows
     ## Now add rows for linear extrapolation...
     ind <- x < ll 
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ll)%*%D[1:2,]
     ind <- x > ul
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ul)%*%D[3:4,]
     X <- X%*%Sig 
    ## X <- sweep(X,2,c(0,object$cmX))
  }
  X <- sweep(X,2,c(0,object$cmX))
  X
}

###########################################################
### Adding increasing & convex SCOP-spline construction without identifiability constraints
### to be used with numeric 'by' variable and linear functional terms...
##########################################################

smooth.construct.micxBy.smooth.spec<- function(object, data, knots)
## construction of the monotone increasing and convex smooth
{ 
  m <- object$p.order[1]
  if (is.na(m)) m <- 2 ## default 
  if (m < 1) stop("silly m supplied")
  if (object$bs.dim<0) object$bs.dim <- 10 ## default
  q <- object$bs.dim 
  nk <- q+m+2 ## number of knots
  if (nk<=0) stop("k too small for m")
  x <- data[[object$term]]  ## the data
  xk <- knots[[object$term]] ## will be NULL if none supplied
  if (is.null(xk)) # space knots through data
  { n<-length(x)
    xk<-rep(0,q+m+2)
    xk[(m+2):(q+1)]<-seq(min(x),max(x),length=q-m)
    for (i in 1:(m+1)) {xk[i]<-xk[m+2]-(m+2-i)*(xk[m+3]-xk[m+2])}
    for (i in (q+2):(q+m+2)) {xk[i]<-xk[q+1]+(i-q-1)*(xk[m+3]-xk[m+2])}
  }
  if (length(xk)!=nk) # right number of knots?
  stop(paste("there should be ",nk," supplied knots"))
  ##  get model matrix...
  X <- splineDesign(xk,x,ord=m+2)
  ## Define matrix Sigma for monotone increasing & convex smooth
  Sig <- matrix(0,q,q)   
  Sig[,1] <- rep(1,q)
  for (i in 2:q) Sig[i:q,i] <- c(1:(q-i+1))

  X <- X%*%Sig # model matrix for the constrained term
  object$X <- X 
  object$Sigma <- Sig
  object$P <- list()
  object$S <- list()

  if (!object$fixed) # create the penalty matrix
  { P <- diff(diag(q-2),difference=1)
    P <- rbind(matrix(0,2,q-2),P) ## adding first two rows of zeros
    P <- cbind(matrix(0,q-1,2),P) ## adding first two columns of zeros
    object$P[[1]] <- P
    object$S[[1]] <- crossprod(P)
  }
  object$p.ident <- c(FALSE,rep(TRUE,q-1)) ## p.ident is an indicator of which coefficients must be positive (exponentiated) 
  object$rank <- ncol(object$X)  # penalty rank
  object$null.space.dim <-2 ## m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
 
  ## get model matrix for 1st and 2nd derivatives of the smooth...
  h <- (max(x)-min(x))/(q-m-1) ## distance between two adjacent knots
  object$Xdf1 <- splineDesign(xk,x,ord=m+1)[,1:(q-1)]/h ## ord is by one less for the 1st derivative
  object$Xdf2 <- splineDesign(xk,x,ord=m)[,1:(q-2)]/h^2 ## ord is by two less for the 2nd derivative

  object$knots <- xk;
  object$m <- m;
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  class(object)<- "micxBy.smooth"  # Give object a class
  object
}


## Prediction matrix for the `micxBy` smooth class... 

Predict.matrix.micxBy.smooth<-function(object,data)
## prediction method function for the `micxBy` smooth class
{ m <- object$m+1; # spline order, m+1=3 default for cubic spline
  q <- object$df 
  Sig <- matrix(0,q,q)   
  Sig[,1] <- rep(1,q)
  for (i in 2:q) Sig[i:q,i] <- c(1:(q-i+1))
  
  ## find spline basis inner knot range...
  ll <- object$knots[m+1];ul <- object$knots[length(object$knots)-m]
  m <- m + 1
  x <- data[[object$term]]
  n <- length(x)
  ind <- x<=ul & x>=ll ## data in range
  if (sum(ind)==n) { ## all in range
     X <- spline.des(object$knots,x,m)$design
     X <- X%*%Sig 
  } else { ## some extrapolation needed 
     ## matrix mapping coefs to value and slope at end points...
     D <- spline.des(object$knots,c(ll,ll,ul,ul),m,c(0,1,0,1))$design
     X <- matrix(0,n,ncol(D)) ## full predict matrix
     if (sum(ind)> 0)  X[ind,] <- spline.des(object$knots,x[ind],m)$design ## interior rows
     ## Now add rows for linear extrapolation...
     ind <- x < ll 
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ll)%*%D[1:2,]
     ind <- x > ul
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ul)%*%D[3:4,]
     X <- X%*%Sig 
  }
  X
}


############################################################
### Smooth construct for the convex/concave smooths ......
###########################################################

###########################################################
### Adding concave SCOP-spline construction *************
###########################################################

smooth.construct.cv.smooth.spec<- function(object, data, knots)
## construction of the concave smooth
{ 
 # require(splines)
  m <- object$p.order[1]
  if (is.na(m)) m <- 2 ## default 
  if (m<1) stop("silly m supplied")
  if (object$bs.dim<0) object$bs.dim <- 10 ## default
  q <- object$bs.dim 
  nk <- q+m+2 ## number of knots
  if (nk<=0) stop("k too small for m")
  x <- data[[object$term]]  ## the data
  xk <- knots[[object$term]] ## will be NULL if none supplied
  if (is.null(xk)) # space knots through data
  { n<-length(x)
    xk<-rep(0,q+m+2)
    xk[(m+2):(q+1)]<-seq(min(x),max(x),length=q-m)
    for (i in 1:(m+1)) {xk[i]<-xk[m+2]-(m+2-i)*(xk[m+3]-xk[m+2])}
    for (i in (q+2):(q+m+2)) {xk[i]<-xk[q+1]+(i-q-1)*(xk[m+3]-xk[m+2])}
  }
  if (length(xk)!=nk) # right number of knots?
  stop(paste("there should be ",nk," supplied knots"))
  #  get model matrix-------------
  X1 <- splineDesign(xk,x,ord=m+2)
  # use matrix Sigma and remove the first column for the monotone smooth
  Sig <- matrix(0,(q-1),(q-1))  # Define Sigma for concave smooth
  Sig[1:(q-1),1]<- c(1:(q-1))
  for (i in 2:(q-1)) Sig[i:(q-1),i]<--c(1:(q-i))
  
  X <- X1[,2:q]%*%Sig # model submatrix for the constrained term

  ## applying sum-to-zero (centering) constraint...
  cmx <- colMeans(X)
  X <- sweep(X,2,cmx) ## subtract cmx from columns 
  object$X <- X # the final model matrix
  object$cmX <- cmx

  object$Sigma <- Sig
  object$P <- list()
  object$S <- list()

  if (!object$fixed) # create the penalty matrix
  { P <- diff(diag(q-2),difference=1)
    object$P[[1]] <- P
    S <- matrix(0,q-1,q-1)
    S[2:(q-1),2:(q-1)] <- crossprod(P)
    object$S[[1]] <- S
  }
  object$p.ident <- rep(TRUE,q-1) ## p.ident is an indicator of which coefficients must be positive (exponentiated) 
  object$rank <- ncol(object$X)-1  # penalty rank
  object$null.space.dim <- 2 ##m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  
  ## get model matrix for 1st and 2nd derivatives of the smooth...
  h <- (max(x)-min(x))/(q-m-1) ## distance between two adjacent knots
  object$Xdf1 <- splineDesign(xk,x,ord=m+1)[,2:(q-1)]/h ## ord is by one less for the 1st derivative
  object$Xdf2 <- splineDesign(xk,x,ord=m)[,2:(q-2)]/h^2 ## ord is by two less for the 2nd derivative

  object$knots <- xk;
  object$m <- m;
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  class(object)<-"cv.smooth"  # Give object a class
  object
}


## Prediction matrix for the `cv` smooth class******************

Predict.matrix.cv.smooth<-function(object,data)
## prediction method function for the `cv' smooth class
{ m <- object$m+1; # spline order, m+1=3 default for cubic spline
  q <- object$df +1
  Sig <- matrix(0,q,q)   
  Sig[,1] <- rep(1,q)
  Sig[2:q,2]<- c(1:(q-1))
  for (i in 3:q) Sig[i:q,i] <- -c(1:(q-i+1))

  ## find spline basis inner knot range...
  ll <- object$knots[m+1];ul <- object$knots[length(object$knots)-m]
  m <- m + 1
  x <- data[[object$term]]
  n <- length(x)
  ind <- x<=ul & x>=ll ## data in range
  if (sum(ind)==n) { ## all in range
     X <- spline.des(object$knots,x,m)$design
     X <- X%*%Sig 
  } else { ## some extrapolation needed 
     ## matrix mapping coefs to value and slope at end points...
     D <- spline.des(object$knots,c(ll,ll,ul,ul),m,c(0,1,0,1))$design
     X <- matrix(0,n,ncol(D)) ## full predict matrix
     if (sum(ind)> 0)  X[ind,] <- spline.des(object$knots,x[ind],m)$design ## interior rows
     ## Now add rows for linear extrapolation...
     ind <- x < ll 
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ll)%*%D[1:2,]
     ind <- x > ul
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ul)%*%D[3:4,]
     X <- X%*%Sig 
    ## X <- sweep(X,2,c(0,object$cmX))
  }
  X <- sweep(X,2,c(0,object$cmX))
  X
}

####################################################################################################
### Adding concave SCOP-spline construction without applying identifiability constraints
### to be used with numeric 'by' variable...
####################################################################################################
## when 'by' variable takes more than one value, the smooth terms are identifiable without a
## 'zero intercept' constraint, so they are left unconstrained...

smooth.construct.cvBy.smooth.spec<- function(object, data, knots)
## construction of the concave smooth
{ 
 # require(splines)
  m <- object$p.order[1]
  if (is.na(m)) m <- 2 ## default 
  if (m<1) stop("silly m supplied")
  if (object$bs.dim<0) object$bs.dim <- 10 ## default
  q <- object$bs.dim 
  nk <- q+m+2 ## number of knots
  if (nk<=0) stop("k too small for m")
  x <- data[[object$term]]  ## the data
  xk <- knots[[object$term]] ## will be NULL if none supplied
  if (is.null(xk)) # space knots through data
  { n<-length(x)
    xk<-rep(0,q+m+2)
    xk[(m+2):(q+1)]<-seq(min(x),max(x),length=q-m)
    for (i in 1:(m+1)) {xk[i]<-xk[m+2]-(m+2-i)*(xk[m+3]-xk[m+2])}
    for (i in (q+2):(q+m+2)) {xk[i]<-xk[q+1]+(i-q-1)*(xk[m+3]-xk[m+2])}
  }
  if (length(xk)!=nk) # right number of knots?
  stop(paste("there should be ",nk," supplied knots"))
  #  get model matrix-------------
  X <- splineDesign(xk,x,ord=m+2)
  # get matrix Sigma for concave smooth...
  Sig <- matrix(0,q,q)   
  Sig[,1] <- rep(1,q)
  Sig[2:q,2]<- c(1:(q-1))
  for (i in 3:q) Sig[i:q,i] <- -c(1:(q-i+1))
  
  X <- X%*%Sig # model matrix for the constrained term
  object$X<-X 
  object$Sigma <- Sig
  object$P <- list()
  object$S <- list()

  if (!object$fixed) # create the penalty matrix
  { P <- diff(diag(q-2),difference=1)
    P <- rbind(matrix(0,2,q-2),P) ## adding first two rows of zeros
    P <- cbind(matrix(0,q-1,2),P) ## adding first two columns of zeros
    object$P[[1]] <- P
    object$S[[1]] <- crossprod(P)
  }
  object$p.ident <- c(FALSE,rep(TRUE,q-1)) ## p.ident is an indicator of which coefficients must be positive (exponentiated) 
  object$rank <- ncol(object$X)  # penalty rank
  object$null.space.dim <- 2 ##m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  
  ## get model matrix for 1st and 2nd derivatives of the smooth...
  h <- (max(x)-min(x))/(q-m-1) ## distance between two adjacent knots
  object$Xdf1 <- splineDesign(xk,x,ord=m+1)[,1:(q-1)]/h ## ord is by one less for the 1st derivative
  object$Xdf2 <- splineDesign(xk,x,ord=m)[,1:(q-2)]/h^2 ## ord is by two less for the 2nd derivative

  object$knots <- xk;
  object$m <- m;
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  class(object)<-"cvBy.smooth"  # Give object a class
  object
}


## Prediction matrix for the `cvBy` smooth class******************

Predict.matrix.cvBy.smooth<-function(object,data)
## prediction method function for the `cvBy' smooth class
{ m <- object$m+1; # spline order, m+1=3 default for cubic spline
  q <- object$df 
  Sig <- matrix(0,q,q)   
  Sig[,1] <- rep(1,q)
  Sig[2:q,2]<- c(1:(q-1))
  for (i in 3:q) Sig[i:q,i] <- -c(1:(q-i+1))

  ## find spline basis inner knot range...
  ll <- object$knots[m+1];ul <- object$knots[length(object$knots)-m]
  m <- m + 1
  x <- data[[object$term]]
  n <- length(x)
  ind <- x<=ul & x>=ll ## data in range
  if (sum(ind)==n) { ## all in range
     X <- spline.des(object$knots,x,m)$design
     X <- X%*%Sig 
  } else { ## some extrapolation needed 
     ## matrix mapping coefs to value and slope at end points...
     D <- spline.des(object$knots,c(ll,ll,ul,ul),m,c(0,1,0,1))$design
     X <- matrix(0,n,ncol(D)) ## full predict matrix
     if (sum(ind)> 0)  X[ind,] <- spline.des(object$knots,x[ind],m)$design ## interior rows
     ## Now add rows for linear extrapolation...
     ind <- x < ll 
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ll)%*%D[1:2,]
     ind <- x > ul
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ul)%*%D[3:4,]
     X <- X%*%Sig 
  }
  X
}



###########################################################
### Adding convex SCOP-spline construction *************
##########################################################

smooth.construct.cx.smooth.spec<- function(object, data, knots)
## construction of the convex smooth
{ 
  # require(splines)
  m <- object$p.order[1]
  if (is.na(m)) m <- 2 ## default 
  if (m<1) stop("silly m supplied")
  if (object$bs.dim<0) object$bs.dim <- 10 ## default
  q <- object$bs.dim 
  nk <- q+m+2 ## number of knots
  if (nk<=0) stop("k too small for m")
  x <- data[[object$term]]  ## the data
  xk <- knots[[object$term]] ## will be NULL if none supplied
  if (is.null(xk)) # space knots through data
  { n<-length(x)
    xk<-rep(0,q+m+2)
    xk[(m+2):(q+1)]<-seq(min(x),max(x),length=q-m)
    for (i in 1:(m+1)) {xk[i]<-xk[m+2]-(m+2-i)*(xk[m+3]-xk[m+2])}
    for (i in (q+2):(q+m+2)) {xk[i]<-xk[q+1]+(i-q-1)*(xk[m+3]-xk[m+2])}
  }
  if (length(xk)!=nk) # right number of knots?
  stop(paste("there should be ",nk," supplied knots"))
  #  get model matrix-------------
  X1 <- splineDesign(xk,x,ord=m+2)
  # use matrix Sigma and remove the first column for the monotone smooth
  Sig <- matrix(0,(q-1),(q-1))   # Define Sigma for convex smooth
  Sig[1:(q-1),1]<- -c(1:(q-1))
  for (i in 2:(q-1)) Sig[i:(q-1),i]<- c(1:(q-i))
  
  X <- X1[,2:q]%*%Sig # model submatrix for the constrained term

  ## applying sum-to-zero (centering) constraint...
  cmx <- colMeans(X)
  X <- sweep(X,2,cmx) ## subtract cmx from columns 
  object$X <- X # the final model matrix
  object$cmX <- cmx

  object$Sigma <- Sig
  object$P <- list()
  object$S <- list()

  if (!object$fixed) # create the penalty matrix
  { P <- diff(diag(q-2),difference=1)
    object$P[[1]] <- P
    S <- matrix(0,q-1,q-1)
    S[2:(q-1),2:(q-1)] <- crossprod(P)
    object$S[[1]] <- S
  }
  object$p.ident <- rep(TRUE,q-1) ## p.ident is an indicator of which coefficients must be positive (exponentiated) 
  object$rank <- ncol(object$X)-1  # penalty rank
  object$null.space.dim <- 2 ##m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  
  ## get model matrix for 1st and 2nd derivatives of the smooth...
  h <- (max(x)-min(x))/(q-m-1) ## distance between two adjacent knots
  object$Xdf1 <- splineDesign(xk,x,ord=m+1)[,2:(q-1)]/h ## ord is by one less for the 1st derivative
  object$Xdf2 <- splineDesign(xk,x,ord=m)[,2:(q-2)]/h^2 ## ord is by two less for the 2nd derivative

  object$knots <- xk;
  object$m <- m;
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  class(object)<-"cx.smooth"  # Give object a class
  object
}


## Prediction matrix for the `cx` smooth class *************************

Predict.matrix.cx.smooth<-function(object,data)
## prediction method function for the `cx' smooth class
{ m <- object$m+1; # spline order, m+1=3 default for cubic spline
  q <- object$df +1
  Sig <- matrix(0,q,q)   
  Sig[,1] <- rep(1,q)
  Sig[2:q,2]<- -c(1:(q-1))
  for (i in 3:q) Sig[i:q,i] <- c(1:(q-i+1))
  
  ## find spline basis inner knot range...
  ll <- object$knots[m+1];ul <- object$knots[length(object$knots)-m]
  m <- m + 1
  x <- data[[object$term]]
  n <- length(x)
  ind <- x<=ul & x>=ll ## data in range
  if (sum(ind)==n) { ## all in range
     X <- spline.des(object$knots,x,m)$design
     X <- X%*%Sig 
  } else { ## some extrapolation needed 
     ## matrix mapping coefs to value and slope at end points...
     D <- spline.des(object$knots,c(ll,ll,ul,ul),m,c(0,1,0,1))$design
     X <- matrix(0,n,ncol(D)) ## full predict matrix
     if (sum(ind)> 0)  X[ind,] <- spline.des(object$knots,x[ind],m)$design ## interior rows
     ## Now add rows for linear extrapolation...
     ind <- x < ll 
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ll)%*%D[1:2,]
     ind <- x > ul
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ul)%*%D[3:4,]
     X <- X%*%Sig 
    ## X <- sweep(X,2,c(0,object$cmX))
  }
  X <- sweep(X,2,c(0,object$cmX))
  X
}

###########################################################
### Adding convex SCOP-spline construction without identifiability constraints
### to be used with numeric 'by' variable and linear functional terms...
##########################################################

smooth.construct.cxBy.smooth.spec<- function(object, data, knots)
## construction of the convex smooth
{ 
  m <- object$p.order[1]
  if (is.na(m)) m <- 2 ## default 
  if (m<1) stop("silly m supplied")
  if (object$bs.dim<0) object$bs.dim <- 10 ## default
  q <- object$bs.dim 
  nk <- q+m+2 ## number of knots
  if (nk<=0) stop("k too small for m")
  x <- data[[object$term]]  ## the data
  xk <- knots[[object$term]] ## will be NULL if none supplied
  if (is.null(xk)) # space knots through data
  { n<-length(x)
    xk<-rep(0,q+m+2)
    xk[(m+2):(q+1)]<-seq(min(x),max(x),length=q-m)
    for (i in 1:(m+1)) {xk[i]<-xk[m+2]-(m+2-i)*(xk[m+3]-xk[m+2])}
    for (i in (q+2):(q+m+2)) {xk[i]<-xk[q+1]+(i-q-1)*(xk[m+3]-xk[m+2])}
  }
  if (length(xk)!=nk) # right number of knots?
  stop(paste("there should be ",nk," supplied knots"))
  #  get model matrix-------------
  X <- splineDesign(xk,x,ord=m+2)
  # get matrix Sigma for convex  smooth...
  Sig <- matrix(0,q,q)   
  Sig[,1] <- 1
  Sig[2:q,2]<- -c(1:(q-1))
  for (i in 3:q) Sig[i:q,i] <- c(1:(q-i+1))
  
  X <- X%*%Sig # model submatrix for the constrained term
  object$X<-X # the final model matrix
  object$Sigma <- Sig
  object$P <- list()
  object$S <- list()

  if (!object$fixed) # create the penalty matrix
  { P <- diff(diag(q-2),difference=1)
    P <- rbind(matrix(0,2,q-2),P) ## adding first two rows of zeros
    P <- cbind(matrix(0,q-1,2),P) ## adding first two columns of zeros
    object$P[[1]] <- P
    object$S[[1]] <- crossprod(P)
  }
  object$p.ident <- c(FALSE,rep(TRUE,q-1)) ## p.ident is an indicator of which coefficients must be positive (exponentiated) 
  object$rank <- ncol(object$X)  # penalty rank
  object$null.space.dim <- 2 ##m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  
  ## get model matrix for 1st and 2nd derivatives of the smooth...
  h <- (max(x)-min(x))/(q-m-1) ## distance between two adjacent knots
  object$Xdf1 <- splineDesign(xk,x,ord=m+1)[,1:(q-1)]/h ## ord is by one less for the 1st derivative
  object$Xdf2 <- splineDesign(xk,x,ord=m)[,1:(q-2)]/h^2 ## ord is by two less for the 2nd derivative

  object$knots <- xk;
  object$m <- m;
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  class(object)<-"cxBy.smooth"  # Give object a class
  object
}


## Prediction matrix for the `cxBy` smooth class *************************

Predict.matrix.cxBy.smooth<-function(object,data)
## prediction method function for the `cxBy' smooth class
{ m <- object$m+1; # spline order, m+1=3 default for cubic spline
  q <- object$df 
  Sig <- matrix(0,q,q)   
  Sig[,1] <- rep(1,q)
  Sig[2:q,2]<- -c(1:(q-1))
  for (i in 3:q) Sig[i:q,i] <- c(1:(q-i+1))
  
  ## find spline basis inner knot range...
  ll <- object$knots[m+1];ul <- object$knots[length(object$knots)-m]
  m <- m + 1
  x <- data[[object$term]]
  n <- length(x)
  ind <- x<=ul & x>=ll ## data in range
  if (sum(ind)==n) { ## all in range
     X <- spline.des(object$knots,x,m)$design
     X <- X%*%Sig 
  } else { ## some extrapolation needed 
     ## matrix mapping coefs to value and slope at end points...
     D <- spline.des(object$knots,c(ll,ll,ul,ul),m,c(0,1,0,1))$design
     X <- matrix(0,n,ncol(D)) ## full predict matrix
     if (sum(ind)> 0)  X[ind,] <- spline.des(object$knots,x[ind],m)$design ## interior rows
     ## Now add rows for linear extrapolation...
     ind <- x < ll 
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ll)%*%D[1:2,]
     ind <- x > ul
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ul)%*%D[3:4,]
     X <- X%*%Sig 
  }
  X
}


#####################################################
### Adding positive SCOP-spline construction 
######################################################


smooth.construct.po.smooth.spec<- function(object, data, knots)
## construction of the positivelt costrained smooth
{ 
  m <- object$p.order[1]
  if (is.na(m)) m <- 2 ## default for cubic spline
  if (m<1) stop("silly m supplied")
  if (object$bs.dim<0) object$bs.dim <- 10 ## default
  q <- object$bs.dim 
  nk <- q+m+2 ## number of knots
  if (nk<=0) stop("k too small for m")
  x <- data[[object$term]]  ## the data
  xk <- knots[[object$term]] ## will be NULL if none supplied
  if (is.null(xk)) # space knots through data
     { n<-length(x)
       xk<-rep(0,q+m+2)
       xk[(m+2):(q+1)]<-seq(min(x),max(x),length=q-m)
       for (i in 1:(m+1)) {xk[i]<-xk[m+2]-(m+2-i)*(xk[m+3]-xk[m+2])}
       for (i in (q+2):(q+m+2)) {xk[i]<-xk[q+1]+(i-q-1)*(xk[m+3]-xk[m+2])}
     }
  if (length(xk)!=nk) # right number of knots?
  stop(paste("there should be ", nk," supplied knots"))
  #  get model matrix-------------
  X1 <- splineDesign(xk,x,ord=m+2)
  # use matrix Sigma and remove the first column for the positive smooth
  Sig <- diag(1,(q-1))  ## identity matrix 
 ## Sig <- diag(1,q)  ## identity matrix
     
  X <- X1[,2:q] 
 ## X <- X1
  object$X <- X # the final model matrix
  object$P <- list()
  object$S <- list()
  object$Sigma <- Sig

  if (!object$fixed) # create the penalty matrix
  { P <- diff(diag(q-1),difference=1)
   ## P <- diff(diag(q),difference=1)
    object$P[[1]] <- P
    object$S[[1]] <- crossprod(P)
  }
  object$p.ident <- rep(TRUE,q-1) ## p.ident is an indicator of which coefficients must be positive (exponentiated)
 ## object$p.ident <- rep(TRUE,q)
  object$rank <- ncol(object$X)-1  # penalty rank
 ## object$rank <- ncol(object$X)  # penalty rank
  object$null.space.dim <-2 ## m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  object$knots <- xk;
  object$m <- m;
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)

  ## get model matrix for 1st and 2nd derivatives of the smooth...
  h <- (max(x)-min(x))/(q-m-1) ## distance between two adjacent knots
  object$Xdf1 <- splineDesign(xk,x,ord=m+1)[,2:(q-1)]/h ## ord is by one less for the 1st derivative
  object$Xdf2 <- splineDesign(xk,x,ord=m)[,2:(q-2)]/h^2 ## ord is by two less for the 2nd derivative

  class(object)<-"po.smooth"  # Give object a class
  object
}

## NOte: maybe no need to set the first coefficient to zero, as no intercept is required in the model with 'po' smooth
## Prediction matrix for the `po` smooth class... 

Predict.matrix.po.smooth<-function(object,data)
## prediction method function for the `po' smooth class
{ m <- object$m+1; # spline order, m+1=3 default for cubic spline
  q <- object$df +1
 ## Sig <- diag(1,q)   # Define Matrix Sigma
  ## find spline basis inner knot range...
  ll <- object$knots[m+1];ul <- object$knots[length(object$knots)-m]
  m <- m + 1
  x <- data[[object$term]]
  n <- length(x)
  ind <- x<=ul & x>=ll ## data in range
  if (sum(ind)==n) { ## all in range
     X <- spline.des(object$knots,x,m)$design
  } else { ## some extrapolation needed 
     ## matrix mapping coefs to value and slope at end points...
     D <- spline.des(object$knots,c(ll,ll,ul,ul),m,c(0,1,0,1))$design
     X <- matrix(0,n,ncol(D)) ## full predict matrix
     if (sum(ind)> 0)  X[ind,] <- spline.des(object$knots,x[ind],m)$design ## interior rows
     ## Now add rows for linear extrapolation...
     ind <- x < ll 
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ll)%*%D[1:2,]
     ind <- x > ul
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ul)%*%D[3:4,]
  }
  X
}






