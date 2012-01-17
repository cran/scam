

### Adding Monotone increasing P-spline construction *************


smooth.construct.mpi.smooth.spec<- function(object, data, knots)
## construction of the monotone increasing smooth
{ 
  require(splines)
  m <- object$p.order[1]
  if (is.na(m)) m <- 2 ## default 
  if (m<1) stop("silly m supplied")
  if (object$bs.dim<0) object$bs.dim <- 10 ## default
  q <- object$bs.dim 
  nk <- q+m+2 ## number of knots
  if (nk<=0) stop("k too small for m")
  x <- data[[object$term]]  ## the data
 # x.shift <- mean(x) # shift used to enhance stability
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
#  x <- x - x.shift # basis stabilizing shift
#  xk <- xk - x.shift # knots treated the same
  #  get model matrix-------------
  X1 <- splineDesign(xk,x,ord=m+2)
  # use matrix Sigma and remove the first column for the monotone smooth
  Sig <- matrix(0,(q-1),(q-1))   # Define Matrix Sigma
     # elements of matrix Sigma for increasing smooth
  for (i in 1:(q-1))  Sig[i,1:i]<-1
  X <- X1[,2:q]%*%Sig # model submatrix for the monotone term
X <- sweep(X,2,colMeans(X))
  object$X<-X # the finished model matrix
  object$P <- list()
  object$S <- list()
  object$Sigma <- Sig

  if (!object$fixed) # create the penalty matrix
  { P <- diff(diag(q-1),difference=1)
    object$P[[1]] <- P
    object$S[[1]] <- t(P)%*%P
  }
  b<-rep(1,q-1) # define vector of 0's & 1's for model parameters identification
  object$p.ident <- b  
  object$rank <- ncol(object$X)-1  # penalty rank
  object$null.space.dim <- m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  
  ## store "mpi" specific stuff ...
  object$knots <- xk;
  object$m <- m;
#  object$x.shift <- x.shift
 
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  class(object)<-"mpi.smooth"  # Give object a class
  object
}


## Prediction matrix for the `mpi` smooth class *************************


Predict.matrix.mpi.smooth<-function(object,data)
## prediction method function for the `mpi' smooth class
{ x <- data[[object$term]]
 # x <- x - object$x.shift # stabilizing shift
  m <- object$m;     # spline order (3=cubic)
  q <- object$df + 1
  xk<-object$knots    # knot locations
  nk<-length(xk)      # number of knots
  X1<-splineDesign(xk,x,ord=m+2)
  Sig<-matrix(0,q,q)   # Define Matrix Sigma
   # elements of matrix Sigma for increasing smooth
  for (i in 1:q)  Sig[i,1:i]<-1
 # X<-X1[,2:ncol(X1)]%*%Sig 
  X<-X1%*%Sig 
# X <- sweep(X,2,colMeans(X)) # no identifiability constraints for predict method
  X # return the prediction matrix
}


### Adding Monotone decreasing P-spline construction *************


smooth.construct.mpd.smooth.spec<- function(object, data, knots)
## construction of the monotone decreasing smooth
{ 
  require(splines)
  m <- object$p.order[1]
  if (is.na(m)) m <- 2 ## default 
  if (m<1) stop("silly m supplied")
  if (object$bs.dim<0) object$bs.dim <- 10 ## default
  q <- object$bs.dim 
  nk <- q+m+2 ## number of knots
  if (nk<=0) stop("k too small for m")
  x <- data[[object$term]]  ## the data
#  x.shift <- mean(x) # shift used to enhance stability
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
#  x <- x - x.shift # basis stabilizing shift
#  xk <- xk - x.shift # knots treated the same
  #  get model matrix-------------
  X1 <- splineDesign(xk,x,ord=m+2)
  # use matrix Sigma and remove the first column for the monotone smooth
  Sig <- matrix(0,(q-1),(q-1))   # Define Matrix Sigma
      # elements of matrix Sigma for decreasing smooth
  for (i in 1:(q-1))  Sig[i,1:i]<- -1
  X <- X1[,2:q]%*%Sig # model submatrix for the monotone term
  X <- sweep(X,2,colMeans(X))
  object$X<-X # the finished model matrix
  object$Sigma <- Sig
  object$P <- list()
  object$S <- list()

  if (!object$fixed) # create the penalty matrix
  { P <- diff(diag(q-1),difference=1)
    object$P[[1]] <- P
    object$S[[1]] <- t(P)%*%P
  }
  b<-rep(1,q-1) # define vector of 0's & 1's for model parameters identification
  object$p.ident <- b  
  object$rank <- ncol(object$X)-1  # penalty rank
  object$null.space.dim <- m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  
  ## store "mpd" specific stuff ...
  object$knots <- xk;
  object$m <- m;
#  object$x.shift <- x.shift
 
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  class(object)<-"mpd.smooth"  # Give object a class
  object
}


## Prediction matrix for the `mpd` smooth class *************************

Predict.matrix.mpd.smooth<-function(object,data)
## prediction method function for the `mpd' smooth class
{ x <- data[[object$term]]
  m <- object$m;     # spline order (3=cubic)
  q <- object$df +1
  xk <- object$knots    # knot locations
  nk <- length(xk)      # number of knots
  X1 <- splineDesign(xk,x,ord=m+2)
   # use matrix Sigma and remove the first column for the monotone smooth
  Sig <- matrix(0,q,q)   # Define Matrix Sigma
      # elements of matrix Sigma for decreasing smooth
  Sig[,1] <- 1
  for (i in 2:q)  Sig[i,2:i]<- -1
  X <- X1%*%Sig 
#  X <- sweep(X,2,colMeans(X))
  X # return the prediction matrix
}






### Smooth construct for the mixed constrainted smooth ......

### Adding Monotone decreasing & concave P-spline construction *************


smooth.construct.mdcv.smooth.spec<- function(object, data, knots)
## construction of the monotone decreasing and concave smooth
{ 
  require(splines)
  m <- object$p.order[1]
  if (is.na(m)) m <- 2 ## default 
  if (m<1) stop("silly m supplied")
  if (object$bs.dim<0) object$bs.dim <- 10 ## default
  q <- object$bs.dim 
  nk <- q+m+2 ## number of knots
  if (nk<=0) stop("k too small for m")
  x <- data[[object$term]]  ## the data
#  x.shift <- mean(x) # shift used to enhance stability
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
#  x <- x - x.shift # basis stabilizing shift
#  xk <- xk - x.shift # knots treated the same
  #  get model matrix-------------
  X1 <- splineDesign(xk,x,ord=m+2)
  # use matrix Sigma and remove the first column for the monotone smooth
  Sig <- matrix(0,(q-1),(q-1))   # Define Matrix Sigma
     # for monotone decreasing & concave smooth
  for (i in 1:(q-1)) Sig[i:(q-1),i]<--c(1:(q-i))
  X <- X1[,2:q]%*%Sig # model submatrix for the constrained term
  X <- sweep(X,2,colMeans(X))
  object$X<-X # the finished model matrix
  object$Sigma <- Sig
  object$P <- list()
  object$S <- list()

  if (!object$fixed) # create the penalty matrix
  { P <- diff(diag(q-2),difference=1)
    object$P[[1]] <- P
    S <- matrix(0,q-1,q-1)
    S[2:(q-1),2:(q-1)] <- t(P)%*%P
    object$S[[1]] <- S
  }
  b<-rep(1,q-1) # define vector of 0's & 1's for model parameters identification
  object$p.ident <- b  
  object$rank <- ncol(object$X)-1  # penalty rank
  object$null.space.dim <- m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  
  ## store "mdcv" specific stuff ...
  object$knots <- xk;
  object$m <- m;
#  object$x.shift <- x.shift
 
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  class(object)<-"mdcv.smooth"  # Give object a class
  object
}


## Prediction matrix for the `mdcv` smooth class *************************

Predict.matrix.mdcv.smooth<-function(object,data)
## prediction method function for the `mdcv' smooth class
{ x <- data[[object$term]]
  m <- object$m;     # spline order (3=cubic)
  q <- object$df +1
  xk <- object$knots    # knot locations
  nk <- length(xk)      # number of knots
  X1 <- splineDesign(xk,x,ord=m+2)
  # use matrix Sigma and remove the first column for the monotone smooth
  Sig1 <- matrix(0,(q-1),(q-1))   # Define Matrix Sigma
     # for monotone decreasing & concave smooth
  for (i in 1:(q-1)) Sig1[i:(q-1),i] <- -c(1:(q-i))
  Sig <- matrix(0,q,q)
  Sig[,1] <- rep(1,q)
  Sig [2:q,2:q] <- Sig1 
  X <- X1%*%Sig 
#  X <- sweep(X,2,colMeans(X))
  X # return the prediction matrix
}



### Smooth construct for the mixed constrainted smooth ......
### Adding Monotone decreasing & convex P-spline construction *************


smooth.construct.mdcx.smooth.spec<- function(object, data, knots)
## construction of the monotone decreasing and convex smooth
{ 
  require(splines)
  m <- object$p.order[1]
  if (is.na(m)) m <- 2 ## default 
  if (m<1) stop("silly m supplied")
  if (object$bs.dim<0) object$bs.dim <- 10 ## default
  q <- object$bs.dim 
  nk <- q+m+2 ## number of knots
  if (nk<=0) stop("k too small for m")
  x <- data[[object$term]]  ## the data
#  x.shift <- mean(x) # shift used to enhance stability
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
#  x <- x - x.shift # basis stabilizing shift
#  xk <- xk - x.shift # knots treated the same
  #  get model matrix-------------
  X1 <- splineDesign(xk,x,ord=m+2)
  # use matrix Sigma and remove the first column for the monotone smooth
  Sig <- matrix(0,(q-1),(q-1))   # Define Matrix Sigma
     # for monotone decreasing & convex smooth
  Sig[1,]<--rep(1,q-1)
  for (i in 2:(q-1)) {
         Sig[i,1:(q-i)]<--i;
         Sig[i,(q-i+1):(q-1)]<--c((i-1):1)
  }
  X <- X1[,2:q]%*%Sig # model submatrix for the constrained term
  X <- sweep(X,2,colMeans(X))
  object$X<-X # the finished model matrix
  object$Sigma <- Sig
  object$P <- list()
  object$S <- list()

  if (!object$fixed) # create the penalty matrix
  { P <- diff(diag(q-2),difference=1)
    object$P[[1]] <- P
    S <- matrix(0,q-1,q-1)
    S[2:(q-1),2:(q-1)] <- t(P)%*%P
    object$S[[1]] <- S
  }
  b<-rep(1,q-1) # define vector of 0's & 1's for model parameters identification
  object$p.ident <- b  
  object$rank <- ncol(object$X)-1  # penalty rank
  object$null.space.dim <- m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  
  ## store "mdcx" specific stuff ...
  object$knots <- xk;
  object$m <- m;
#  object$x.shift <- x.shift
 
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  class(object)<-"mdcx.smooth"  # Give object a class
  object
}


## Prediction matrix for the `mdcx` smooth class *************************

Predict.matrix.mdcx.smooth<-function(object,data)
## prediction method function for the `mdcx' smooth class
{ x <- data[[object$term]]
  m <- object$m;     # spline order (3=cubic)
  q <- object$df +1
  xk<-object$knots    # knot locations
  nk<-length(xk)      # number of knots
  X1<-splineDesign(xk,x,ord=m+2)
  
  Sig1 <- matrix(0,(q-1),(q-1))   # Define Matrix Sigma
     # for monotone decreasing & convex smooth
  Sig1[1,] <- -rep(1,q-1)
  for (i in 2:(q-1)) {
         Sig1[i,1:(q-i)]<--i;
         Sig1[i,(q-i+1):(q-1)]<--c((i-1):1)
  }
  Sig <- matrix(0,q,q)
  Sig[,1] <- rep(1,q)
  Sig [2:q,2:q] <- Sig1
  X <- X1%*%Sig 
#  X <- sweep(X,2,colMeans(X))
  X # return the prediction matrix
}





### Smooth construct for the mixed constrainted smooth ......
### Adding Monotone increasing & concave P-spline construction *************


smooth.construct.micv.smooth.spec<- function(object, data, knots)
## construction of the monotone increasing and concave smooth
{ 
  require(splines)
  m <- object$p.order[1]
  if (is.na(m)) m <- 2 ## default 
  if (m<1) stop("silly m supplied")
  if (object$bs.dim<0) object$bs.dim <- 10 ## default
  q <- object$bs.dim 
  nk <- q+m+2 ## number of knots
  if (nk<=0) stop("k too small for m")
  x <- data[[object$term]]  ## the data
#  x.shift <- mean(x) # shift used to enhance stability
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
#  x <- x - x.shift # basis stabilizing shift
#  xk <- xk - x.shift # knots treated the same
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
  X <- sweep(X,2,colMeans(X))
  object$X<-X # the finished model matrix
  object$Sigma <- Sig
  object$P <- list()
  object$S <- list()

  if (!object$fixed) # create the penalty matrix
  { P <- diff(diag(q-2),difference=1)
    object$P[[1]] <- P
    S <- matrix(0,q-1,q-1)
    S[2:(q-1),2:(q-1)] <- t(P)%*%P
    object$S[[1]] <- S
  }
  b<-rep(1,q-1) # define vector of 0's & 1's for model parameters identification
  object$p.ident <- b  
  object$rank <- ncol(object$X)-1  # penalty rank
  object$null.space.dim <- m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  
  ## store "micv" specific stuff ...
  object$knots <- xk;
  object$m <- m;
#  object$x.shift <- x.shift
 
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  class(object)<-"micv.smooth"  # Give object a class
  object
}


## Prediction matrix for the `micv` smooth class *************************

Predict.matrix.micv.smooth<-function(object,data)
## prediction method function for the `micv' smooth class
{ x <- data[[object$term]]
  m <- object$m;     # spline order (3=cubic)
  q <- object$df +1
  xk <- object$knots    # knot locations
  nk <- length(xk)      # number of knots
  X1 <- splineDesign(xk,x,ord=m+2)
  
  Sig1 <- matrix(0,(q-1),(q-1))   # Define Matrix Sigma
     # for monotone increasing & concave smooth
  Sig1[1,] <- rep(1,q-1)
  for (i in 2:(q-1)) {
       Sig1[i,1:(q-i)] <- i;
       Sig1[i,(q-i+1):(q-1)] <- c((i-1):1)
  }
  Sig <- matrix(0,q,q)
  Sig[,1] <- rep(1,q)
  Sig [2:q,2:q] <- Sig1
  X <- X1%*%Sig 
 # X <- sweep(X,2,colMeans(X))
  X # return the prediction matrix
}



### Smooth construct for the mixed constrainted smooth ......
### Adding Monotone increasing & convex P-spline construction *************


smooth.construct.micx.smooth.spec<- function(object, data, knots)
## construction of the monotone increasing and convex smooth
{ 
  require(splines)
  m <- object$p.order[1]
  if (is.na(m)) m <- 2 ## default 
  if (m<1) stop("silly m supplied")
  if (object$bs.dim<0) object$bs.dim <- 10 ## default
  q <- object$bs.dim 
  nk <- q+m+2 ## number of knots
  if (nk<=0) stop("k too small for m")
  x <- data[[object$term]]  ## the data
#  x.shift <- mean(x) # shift used to enhance stability
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
#  x <- x - x.shift # basis stabilizing shift
#  xk <- xk - x.shift # knots treated the same
  #  get model matrix-------------
  X1 <- splineDesign(xk,x,ord=m+2)
  # use matrix Sigma and remove the first column for the monotone smooth
  Sig <- matrix(0,(q-1),(q-1))   # Define Matrix Sigma
     # for monotone increasing & convex smooth
  for (i in 1:(q-1)) Sig[i:(q-1),i]<-c(1:(q-i))
  X <- X1[,2:q]%*%Sig # model submatrix for the constrained term
  X <- sweep(X,2,colMeans(X))
  object$X<-X # the finished model matrix
  object$Sigma <- Sig
  object$P <- list()
  object$S <- list()

  if (!object$fixed) # create the penalty matrix
  { P <- diff(diag(q-2),difference=1)
    object$P[[1]] <- P
    S <- matrix(0,q-1,q-1)
    S[2:(q-1),2:(q-1)] <- t(P)%*%P
    object$S[[1]] <- S
  }
  b<-rep(1,q-1) # define vector of 0's & 1's for model parameters identification
  object$p.ident <- b  
  object$rank <- ncol(object$X)-1  # penalty rank
  object$null.space.dim <- m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
 
  ## store "micx" specific stuff ...
  object$knots <- xk;
  object$m <- m;
#  object$x.shift <- x.shift
 
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  class(object)<-"micx.smooth"  # Give object a class
  object
}


## Prediction matrix for the `micx` smooth class *************************

Predict.matrix.micx.smooth<-function(object,data)
## prediction method function for the `micx` smooth class
{ x <- data[[object$term]]
  m <- object$m;     # spline order (3=cubic)
  q <- object$df +1
  xk <- object$knots    # knot locations
  nk <- length(xk)      # number of knots
  X1 <- splineDesign(xk,x,ord=m+2)
  # use matrix Sigma and remove the first column for the monotone smooth
  Sig1 <- matrix(0,(q-1),(q-1))   # Define Matrix Sigma
     # for monotone increasing & convex smooth
  for (i in 1:(q-1)) Sig1[i:(q-1),i] <- c(1:(q-i))
  Sig <- matrix(0,q,q)
  Sig[,1] <- rep(1,q)
  Sig [2:q,2:q] <- Sig1
  X <- X1%*%Sig 
#  X <- sweep(X,2,colMeans(X))
  X # return the prediction matrix
}







