## (c) Natalya Pya (2023). Provided under GPL 2.
## routines for bivariate tensor product smooth 'interactions' with
## shape constrains



####################################################################################### 
## Tensor product 'interaction' P-spline construction with increasing constraint along the first
## covariate and unconstrained along the 2nd one...
#######################################################################################


smooth.construct.tismi.smooth.spec<- function(object, data, knots)
## construction of the tensor product 'iteraction' bivariate smooth with 
## single monotone increasing constraint wrt 1st covariate ...
{ 
  if (!is.null(object$xt)){
     if (!(object$xt %in% c("ps", "cc")) )
          stop("only 'ps' and 'cc' marginal basis are supported")
      else  bs2 <- object$xt ## basis for the marginal smooth along second direction
  } else bs2 <- "ps"
  if (object$dim !=2)
      stop("the number of covariates should be two")
  if (length(object$p.order)==1)
      {m <- rep(object$p.order, 2) # if a single number is supplied the same
             ## order of P-splines is provided for both marginal smooths
       object$p.order <- m
  }
  else m <- object$p.order
  m[is.na(m)] <- 2  # the default order is 2 (cubic P-spline)
  object$p.order[is.na(object$p.order)] <- 2
  if (object$bs.dim[1]==-1) { # set the default values for q1 and q2
      q1 <- object$bs.dim[1] <- 7
      q2 <- object$bs.dim[2] <- 7
  }
  else if (length(object$bs.dim)==1){
      q1 <- q2 <- object$bs.dim # if `k' is supplied as a single number, the same
             ## basis dimension is provided for both marginal smooths
      object$bs.dim <- rep(object$bs.dim, 2)
  }
  else {q1 <- object$bs.dim[1]; q2 <- object$bs.dim[2]}
  if (is.na(q1)) q1 <- object$bs.dim[1] <- 7  # the default basis dimension is 7
  if (is.na(q2)) q2 <- object$bs.dim[2] <- 7
    
  nk1 <- q1+m[1]+2 ## number of knots for the 1st smooth
  if (bs2=="cc") nk2 <- q2+1
    else nk2 <- q2+m[2]+2 ## number of knots for the 2nd smooth in case of p-splines
  if (nk1<=0 || nk2<=0) stop("either k[1] or k[2] too small for m")
  
  ## the values of the first covariate...   
  x <- data[[object$term[1]]]  
  xk <- knots[[object$term[1]]] ## will be NULL if none supplied
  z <- data[[object$term[2]]]  ## the values of the second covariate
  zk <- knots[[object$term[2]]] ## will be NULL if none supplied
  if (is.null(xk)) # space knots through the values of the 1st covariate
  { 
    xk<-rep(0,q1+m[1]+2)
    xk[(m[1]+2):(q1+1)]<-seq(min(x),max(x),length=q1-m[1])
    for (i in 1:(m[1]+1)) {xk[i]<-xk[m[1]+2]-(m[1]+2-i)*(xk[m[1]+3]-xk[m[1]+2])}
    for (i in (q1+2):(q1+m[1]+2)) {xk[i]<-xk[q1+1]+(i-q1-1)*(xk[m[1]+3]-xk[m[1]+2])}
    knots[[object$term[1]]] <- xk
  }
  n<-length(x)
  if (n != length(z))
     stop ("arguments of smooth not same dimension")
  
  if (is.null(zk)){ # space knots through the values of the 2nd covariate
    if (bs2=="cc") {
         zk <- place.knots(z,nk2) 
         if (length(zk)==2) {
             zk <- place.knots(c(zk,z),nk2)
         }  
    } else{
         zk<-rep(0,q2+m[2]+2)
         zk[(m[2]+2):(q2+1)]<-seq(min(z),max(z),length=q2-m[2])
         for (i in 1:(m[2]+1)) {zk[i]<-zk[m[2]+2]-(m[2]+2-i)*(zk[m[2]+3]-zk[m[2]+2])}
         for (i in (q2+2):(q2+m[2]+2)) {zk[i]<-zk[q2+1]+(i-q2-1)*(zk[m[2]+3]-zk[m[2]+2])}
         knots[[object$term[2]]] <- zk
      }
   }

  if (length(xk)!=nk1 ) # right number of knots?
      stop(paste("there should be ",nk1," supplied knotsfor the x"))
  if (length(zk)!=nk2) # right number of knots?
      stop(paste("there should be ",nk2," supplied knots for z"))
  
  #  get model matrix-------------
  # get marginal model matrices and penalties... 
  if (bs2=="cc") bm <- marginal.matrices.tesmi1.cc(x,z,xk,zk,m,q1,q2)
    else bm <- marginal.matrices.tesmi1.ps(x,z,xk,zk,m,q1,q2)
 
  X1 <- bm$X1
  X2 <- bm$X2
 ## S <- bm$S      
   
  # get a matrix Sigma -----------------------
  IS <- matrix(1,q1,q1)  ## coef summation matrix
  IS[upper.tri(IS)] <-0
  X1 <- X1%*%IS
  X1 <- X1[,-1] ## apply scop indentifiability constraint on the marginal
   
 # cmx <- colMeans(X2)
 # X2 <- sweep(X2,2,cmx) ## apply centering constraint for the unconstrained marginal
 # object$cmX <- cmx
  X2 <- X2[,-1]
  X <- matrix(0,n,(q1-1)*(q2-1)) ## tensor product model matrix
 # X <- matrix(0,n,q1*(q2-1)) 
  for (i in 1:n)
      X[i,] <- X1[i,]%x%X2[i,] # Kronecker product of two rows of marginal model matrices
 
 # IS <- matrix(1,q1,q1)  ## coef summation matrix
 # IS[upper.tri(IS)] <-0
 # I <- diag(q2)
 # Sig <- IS%x%I
 # X <- X%*%Sig

## RE-CONSIDERed D and Penalties, to have no extra constraints!!

  # apply scop identifiability constraint... 
 # D <- diag(q1*q2)
 # D <- D[,-q2]
 # D1 <- t(diff(diag(q2)))
 # D[1:q2,1:(q2-1)] <- D1
 # X <- X%*%D 
  # define identifiability constraint matrix to be used when predicting:
  # D removes the first q2 columns of X and also columns with ind as set below... 
  D <- diag(q1*q2)
  D <- D[,-c(1:q2)]
  ind <- rep(0,q1-1) # get index number for the columns of X to be removed
  for (i in 1:(q1-1))
        ind[i] <- (i-1)*q2+1
  D <- D[,-ind]   
     
  object$X <- X # the final model matrix with identifiability constraint
 
 ## create the penalty matrix... NEED here rather than from marginal.matrices.tesmi1.ps() used for X and S 
  ## since identifiability constraints were applied on the marginal smooths
  S <- list()
  # get the penalty matrix for the first monotone marginal...
  I2<- diag(q2-1)
  P <- diff(diag(q1-1),difference=1) 
  Pm1 <- matrix(0,q1-2,q1-1) # marginal sqrt penalty
  Pm1[1:(q1-2),1:(q1-1)] <- P
  S[[1]]<- Pm1%x%I2
  S[[1]] <- crossprod(S[[1]])  ## t(S[[1]])%*%S[[1]]

 # get penalty for the 2nd marginal...
  I2 <- diff(diag(q2-1),difference=1) 
  I21<- diff(diag(q2-1),difference=2)
  I1 <- diag(q1-1)
  S[[2]] <-matrix(0,q2-3+(q1-2)*(q2-2), (q1-1)*(q2-1))
  S[[2]][1:(q2-3),] <- t(I1[1,])%x%I21
  S[[2]][(q2-2):nrow(S[[2]]),] <- I1[2:(q1-1),]%x%I2
  S[[2]] <- crossprod(S[[2]])  ## t(S[[2]])%*%S[[2]]

  object$S <- list()
 # object$S[[1]] <- crossprod(D,S[[1]])%*%D ## t(D)%*%S[[1]]%*%D
 # object$S[[2]] <- crossprod(D,S[[2]])%*%D ## t(D)%*%S[[2]]%*%D
  object$S[[1]] <- S[[1]]
  object$S[[2]] <- S[[2]]
  object$p.ident <- rep(TRUE,(q1-1)*(q2-1)) ## p.ident is an indicator of which coefficients must be positive (exponentiated)  
  object$p.ident[1:(q2-1)] <- rep(FALSE, q2-1) 
  object$rank <- ncol(object$X)-1  # penalty rank
  object$null.space.dim <-  3 ##  dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  object$Zc <- D   # identifiability constraint matrix

  ## store "tismi" specific stuff ...
  object$knots <- list()
  object$knots[[1]] <- xk
  if (is.null(zk))
     object$knots[[2]] <- rep(0,0,0)
  else object$knots[[2]] <- zk
  object$m <- m
  object$margin.bs <- bs2
       
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  class(object)<-"tismi.smooth"  # Give object a class
  object
}

## Prediction matrix for the `tismi` smooth class *************************

Predict.matrix.tismi.smooth<-function(object,data)
## prediction method function for the `tesmi1' smooth class
{ if (length(object$bs.dim)==1)
      q1 <- q2 <- object$bs.dim # if `k' is supplied as a single number, the same
             ## basis dimension is provided for both marginal smooths
  else  {q1 <- object$bs.dim[1]; q2 <- object$bs.dim[2]}
  
  if (!is.null(object$xt)) bs2 <- object$xt ## basis for the marginal smooth along second direction
    else bs2 <- "ps"
  if (bs2=="cc") bm <- marginal.linear.extrapolation.tesmi1.cc(object, data)
    else bm <- marginal.linear.extrapolation(object, data)

  n <- length(data[[object$term[1]]])
  IS <- matrix(1,q1,q1)  ## coef summation matrix
  IS[upper.tri(IS)] <-0
  X1 <- bm$X1%*%IS
 
  X2 <- bm$X2
#  X2 <- sweep(X2,2,object$cmX) ## apply centering constraint for the unconstrained marginal
  X <- matrix(0,n,q1*q2) ## tensor product model matrix
  for (i in 1:n)
      X[i,] <- X1[i,]%x%X2[i,] # Kronecker product of two rows of marginal model matrices
 
  X  # return the prediction matrix
}

####################################################################################### 
## Tensor product 'interaction' P-spline construction with decreasing constraint along the first
## covariate and unconstrained along the 2nd one...
#######################################################################################

smooth.construct.tismd.smooth.spec<- function(object, data, knots)
## construction of the tensor product 'iteraction' bivariate smooth with 
## single monotone decreasing constraint wrt 1st covariate ...
{ 
  if (!is.null(object$xt)){
     if (!(object$xt %in% c("ps", "cc")) )
          stop("only 'ps' and 'cc' marginal basis are supported")
      else  bs2 <- object$xt ## basis for the marginal smooth along second direction
  } else bs2 <- "ps" ## (only "ps" and 'cc'  available currently)
 
  if (object$dim !=2)
      stop("the number of covariates should be two")
  if (length(object$p.order)==1)
      { m <- rep(object$p.order, 2) # if a single number is supplied the same
             ## order of P-splines is provided for both marginal smooths
        object$p.order <- m
      }
  else m <- object$p.order
  m[is.na(m)] <- 2  # the default order is 2 (cubic P-spline)
  object$p.order[is.na(object$p.order)] <- 2
  if (object$bs.dim[1]==-1)  # set the default values fro q1 and q2
     {  q1 <- object$bs.dim[1] <- 7
        q2 <- object$bs.dim[2] <- 7
     }
  else if (length(object$bs.dim)==1)
         {  q1 <- q2 <- object$bs.dim # if `k' is supplied as a single number, the same
                   ## basis dimension is provided for both marginal smooths
            object$bs.dim <- rep(object$bs.dim, 2)
         }
  else {q1 <- object$bs.dim[1]; q2 <- object$bs.dim[2]}
  if (is.na(q1)) q1 <- object$bs.dim[1] <- 7  # the default basis dimension is 7
  if (is.na(q2)) q2 <- object$bs.dim[2] <- 7
    
  nk1 <- q1+m[1]+2 ## number of knots for the 1st smooth
  if (bs2=="cc") nk2 <- q2+1
    else nk2 <- q2+m[2]+2 ## number of knots for the 2nd smooth in case of p-splines
  if (nk1<=0 || nk2<=0) stop("either k[1] or k[2] too small for m")
  
  ## the values of the first covariate...   
  x <- data[[object$term[1]]]  
  xk <- knots[[object$term[1]]] ## will be NULL if none supplied
  z <- data[[object$term[2]]]  ## the values of the second covariate
  zk <- knots[[object$term[2]]] ## will be NULL if none supplied
  if (is.null(xk)) # space knots through the values of the 1st covariate
     {  xk<-rep(0,q1+m[1]+2)
        xk[(m[1]+2):(q1+1)]<-seq(min(x),max(x),length=q1-m[1])
        for (i in 1:(m[1]+1)) {xk[i]<-xk[m[1]+2]-(m[1]+2-i)*(xk[m[1]+3]-xk[m[1]+2])}
        for (i in (q1+2):(q1+m[1]+2)) {xk[i]<-xk[q1+1]+(i-q1-1)*(xk[m[1]+3]-xk[m[1]+2])}
        knots[[object$term[1]]] <- xk
     }
  n<-length(x)
  if (n != length(z))
       stop ("arguments of smooth not same dimension")

  if (is.null(zk)){ # space knots through the values of the 2nd covariate
    if (bs2=="cc") {
         zk <- place.knots(z,nk2) 
         if (length(zk)==2) {
             zk <- place.knots(c(zk,z),nk2)
         }  
    } else{
         zk<-rep(0,q2+m[2]+2)
         zk[(m[2]+2):(q2+1)]<-seq(min(z),max(z),length=q2-m[2])
         for (i in 1:(m[2]+1)) {zk[i]<-zk[m[2]+2]-(m[2]+2-i)*(zk[m[2]+3]-zk[m[2]+2])}
         for (i in (q2+2):(q2+m[2]+2)) {zk[i]<-zk[q2+1]+(i-q2-1)*(zk[m[2]+3]-zk[m[2]+2])}
         knots[[object$term[2]]] <- zk
      }
   }

  if (length(xk)!=nk1 ) # right number of knots?
      stop(paste("there should be ",nk1," supplied knotsfor the x"))
  if (length(zk)!=nk2) # right number of knots?
      stop(paste("there should be ",nk2," supplied knots for z"))

  #  get model matrix-------------
  # get marginal model matrices and penalties... 
  if (bs2=="cc") bm <- marginal.matrices.tesmi1.cc(x,z,xk,zk,m,q1,q2)
    else bm <- marginal.matrices.tesmi1.ps(x,z,xk,zk,m,q1,q2)
  
  X1 <- bm$X1
  X2 <- bm$X2
     
  # get a matrix Sigma -----------------------
  IS <- matrix(-1,q1,q1)  ## coef summation matrix
  IS[upper.tri(IS)] <-0
  IS[,1] <- -IS[,1]
  X1 <- X1%*%IS
  ## apply scop indentifiability constraint on the marginal...
  X1 <- X1[,-1]

 # cmx <- colMeans(X2)
 # X2 <- sweep(X2,2,cmx) ## apply centering constraint for the unconstrained marginal
 # object$cmX <- cmx

  X2 <- X2[,-1] ## apply zero-intercept constraint 
  X <- matrix(0,n,(q1-1)*(q2-1)) ## tensor product model matrix
  for (i in 1:n)
      X[i,] <- X1[i,]%x%X2[i,] # Kronecker product of two rows of marginal model matrices
 
  # define matrix of identifiability constraints, D, to be used when predicting...
  # D removes the first q2 columns of X and also columns with ind as set below... 
  D <- diag(q1*q2)
  D <- D[,-c(1:q2)]
  ind <- rep(0,q1-1) # get index number for the columns of X to be removed
  for (i in 1:(q1-1))
        ind[i] <- (i-1)*q2+1
  D <- D[,-ind]   
  object$X <- X # the final model matrix with identifiability constraint
 
  ## create the penalty matrix... NEED here rather than from marginal.matrices.tesmi1.ps() used for X and S 
  ## since identifiability constraints were applied on the marginal smooths
  S <- list()
  # get the penalty matrix for the first monotone smooth...
  I2<- diag(q2-1)
  P <- diff(diag(q1-1),difference=1) 
  Pm1 <- matrix(0,q1-2,q1-1) # marginal sqrt penalty
  Pm1[1:(q1-2),1:(q1-1)] <- P
  S[[1]]<- Pm1%x%I2
  S[[1]] <- crossprod(S[[1]])  ## t(S[[1]])%*%S[[1]]

  ## get penalty for the 2nd marginal smooth...
  I2 <- diff(diag(q2-1),difference=1) 
  I21<- diff(diag(q2-1),difference=2)
  I1 <- diag(q1-1)
  S[[2]] <-matrix(0,q2-3+(q1-2)*(q2-2), (q1-1)*(q2-1))
  S[[2]][1:(q2-3),] <- t(I1[1,])%x%I21
  S[[2]][(q2-2):nrow(S[[2]]),] <- I1[2:(q1-1),]%x%I2
  S[[2]] <- crossprod(S[[2]])  ## t(S[[2]])%*%S[[2]]

  object$S <- list()
  object$S[[1]] <- S[[1]]
  object$S[[2]] <- S[[2]]
  object$p.ident <- rep(TRUE,(q1-1)*(q2-1)) ## p.ident is an indicator of which coefficients must be positive (exponentiated)  
  object$rank <- ncol(object$X)-1  # penalty rank
  object$null.space.dim <- 3 ## dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  object$Zc <- D   # identifiability constraint matrix

  ## store "tismd" specific stuff ...
  object$knots <- list()
  object$knots[[1]] <- xk
  if (is.null(zk))
     object$knots[[2]] <- rep(0,0,0)
  else object$knots[[2]] <- zk
  object$m <- m
  object$margin.bs <- bs2
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  
  class(object) <- "tismd.smooth"  # Give object a class
  object
}


###########################################################################

## Prediction matrix for the `tismd` smooth class *************************

Predict.matrix.tismd.smooth<-function(object,data)
{ ## prediction method function for the `tismd1' smooth class
  if (length(object$bs.dim)==1)
      q1 <- q2 <- object$bs.dim # if `k' is supplied as a single number, the same
             ## basis dimension is provided for both marginal smooths
   else  {q1 <- object$bs.dim[1]; q2 <- object$bs.dim[2]}
  
  if (!is.null(object$xt)) bs2 <- object$xt ## basis for the marginal smooth along second direction
    else bs2 <- "ps"
  if (bs2=="cc") bm <- marginal.linear.extrapolation.tesmi1.cc(object, data)
    else bm <- marginal.linear.extrapolation(object, data)

  n <- length(data[[object$term[1]]])
  # get a matrix Sigma -----------------------
  IS <- matrix(-1,q1,q1)  ## coef summation matrix
  IS[upper.tri(IS)] <-0
  IS[,1] <- -IS[,1]
  X1 <- bm$X1%*%IS
  
  X2 <- bm$X2
 # X2 <- sweep(X2,2,object$cmX) ## apply centering constraint for the unconstrained marginal
  X <- matrix(0,n,q1*q2)  # tensor product model matrix
  for (i in 1:n)
      X[i,] <- X1[i,]%x%X2[i,] # Kronecker product of two rows of marginal model matrices
  
  X  # return the prediction matrix
}



