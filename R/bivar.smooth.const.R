#########################################################################
### Shape constrained smooth construct for bivariate terms......       ##
#########################################################################


####################################################################################### 
### Tensor product P-spline construction with double monotone decreasing constraint  ##
#######################################################################################

smooth.construct.tedmd.smooth.spec<- function(object, data, knots)
{ ## construction of the double monotone decreasing smooth surface
 # require(splines)
  if (object$dim !=2)
      stop("the number of covariates should be two")
  if (length(object$p.order)==1)
      m <- rep(object$p.order, 2) # if a single number is supplied the same
             ## order of P-splines is provided for both marginal smooths
  else m <- object$p.order
  m[is.na(m)] <- 2  # the default order is 2 (cubic P-spline)
  if (object$bs.dim[1]==-1)  # set the default values for q1 and q2
     {  q1 <- object$bs.dim[1] <- 7
        q2 <- object$bs.dim[2] <- 7
     }
  else if (length(object$bs.dim)==1)
       q1 <- q2 <- object$bs.dim # if `k' is supplied as a single number, the same
             ## basis dimension is provided for both marginal smooths
  else  {q1 <- object$bs.dim[1]; q2 <- object$bs.dim[2]}
  if (is.na(q1))  q1 <- object$bs.dim[1] <- 7  # the default basis dimension is 7
  if (is.na(q2))  q2 <- object$bs.dim[2] <- 7
  nk1 <- q1+m[1]+2 ## number of knots for the 1st smooth
  nk2 <- q2+m[2]+2 ## number of knots for the 2nd smooth
  if (nk1<=0 || nk2<=0) stop("either k[1] or k[2] too small for m")
  ## the values of the first covariate...   
  x <- data[[object$term[1]]]  
  xk <- knots[[object$term[1]]] ## will be NULL if none supplied
  z <- data[[object$term[2]]]  ## the values of the second covariate
  zk <- knots[[object$term[2]]] ## will be NULL if none supplied
  if (is.null(xk)) # space knots through the values of the 1st covariate
     {  n<-length(x)
        xk<-rep(0,q1+m[1]+2)
        xk[(m[1]+2):(q1+1)]<-seq(min(x),max(x),length=q1-m[1])
        for (i in 1:(m[1]+1)) {xk[i]<-xk[m[1]+2]-(m[1]+2-i)*(xk[m[1]+3]-xk[m[1]+2])}
        for (i in (q1+2):(q1+m[1]+2)) {xk[i]<-xk[q1+1]+(i-q1-1)*(xk[m[1]+3]-xk[m[1]+2])}
     }
  if (n != length(z))
       stop ("arguments of smooth not same dimension")
  if (is.null(zk)) # space knots through the values of the 2nd covariate
     {  zk<-rep(0,q2+m[2]+2)
        zk[(m[2]+2):(q2+1)]<-seq(min(z),max(z),length=q2-m[2])
        for (i in 1:(m[2]+1)) {zk[i]<-zk[m[2]+2]-(m[2]+2-i)*(zk[m[2]+3]-zk[m[2]+2])}
        for (i in (q2+2):(q2+m[2]+2)) {zk[i]<-zk[q2+1]+(i-q2-1)*(zk[m[2]+3]-zk[m[2]+2])}
     }
  if (length(xk)!=nk1 || length(zk)!=nk2) # right number of knots?
       stop(paste("there should be ",nk1, " and ", nk2," supplied knots"))
  #  get model matrix-------------
  X1 <- splineDesign(xk,x,ord=m[1]+2)
  X2 <- splineDesign(zk,z,ord=m[2]+2)
  X <- matrix(0,n,q1*q2)  # model matrix
  for ( i in 1:n)
      {  X[i,] <- X1[i,]%x%X2[i,] # Kronecker product of two rows of marginal model matrices
      }
  # get a matrix Sigma -----------------------
 ## IS <- matrix(0,q2,q2)   # Define submatrix of Sigma
 ## for (j in 1:q2)  IS[j,1:j] <- -1
  IS <- matrix(-1,q2,q2)  ## Define submatrix of Sigma
  IS[upper.tri(IS)] <- 0
 ## IS1 <- matrix(0,q1,q1)   # Define submatrix of Sigma
 ## for (j in 1:q1)  IS1[j,1:j] <- 1
  IS1 <- matrix(1,q1,q1)  ## Define submatrix of Sigma
  IS1[upper.tri(IS1)] <- 0

  Sig <- IS1%x%IS # Knonecker product to get Sigma
  Sig[,1] <- rep(1,ncol(Sig))
 
  # apply identifiability constraint and get model matrix
  X <- X[,2:ncol(X)]%*%Sig[2:ncol(Sig),2:ncol(Sig)]  
  object$X <- X # the finished model matrix with identifiability constraint
 
  # create the penalty matrix
  S <- list()
  I2<- diag(q2)
  I1 <- diff(diag(q1-1),difference=1) 
  Pm1 <- matrix(0,q1-1,q1) # marginal sqrt penalty
  Pm1[2:(q1-1),2:q1] <- I1
  S[[1]]<- Pm1%x%I2

  I2 <- diff(diag(q2-1),difference=1) 
  Pm2 <- matrix(0,q2-1,q2)
  Pm2[2:(q2-1),2:q2] <- I2  # marginal sqrt penalty
  I1 <- diag(q1)
  S[[2]] <- I1%x%Pm2

  object$P <- list()
  object$P[[1]] <- S[[1]][2:nrow(S[[1]]),2:ncol(S[[1]])]
  object$P[[2]] <- S[[2]][2:nrow(S[[2]]),2:ncol(S[[2]])]
  object$S <- list()
  object$S[[1]] <- crossprod(object$P[[1]]) ## t(object$P[[1]])%*%object$P[[1]]
  object$S[[2]] <- crossprod(object$P[[2]]) ## t(object$P[[2]])%*%object$P[[2]]
  object$p.ident <- rep(TRUE,q1*q2-1)  ## p.ident is an indicator of which coefficients must be positive (exponentiated)  
  object$rank <- ncol(object$X)-1  # penalty rank
  object$null.space.dim <- m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  object$Zc <- diag(q1*q2-1) # identfiability constraint matrix
  object$Zc <- rbind(rep(0,ncol(object$Zc)),object$Zc)

  ## store "tedmd" specific stuff ...
  object$knots <- list()
  object$knots[[1]] <- xk
  object$knots[[2]] <- zk
  object$m <- m
   
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  class(object)<-"tedmd.smooth"  # Give object a class
  object
}

####################################################################


Predict.matrix.tedmd.smooth <- function(object, data)
{ ## prediction method function for the `tedmd' smooth class
  if (length(object$bs.dim)==1)
      q1 <- q2 <- object$bs.dim # if `k' is supplied as a single number, the same
             ## basis dimension is provided for both marginal smooths
  else  {q1 <- object$bs.dim[1]; q2 <- object$bs.dim[2]}
  
  bm <- marginal.linear.extrapolation(object, data)
  n <- length(data[[object$term[1]]])
  X <- matrix(0,n,q1*q2)  # model matrix
  for ( i in 1:n)
      {  X[i,] <- bm$X1[i,] %x% bm$X2[i,] # Kronecker product of two rows of marginal model matrices
      }
  # get a matrix Sigma -----------------------
  ## IS <- matrix(0,q2,q2)   # Define submatrix of Sigma
 ## for (j in 1:q2)  IS[j,1:j] <- -1
  IS <- matrix(-1,q2,q2)  ## Define submatrix of Sigma
  IS[upper.tri(IS)] <- 0
 ## IS1 <- matrix(0,q1,q1)   # Define submatrix of Sigma
 ## for (j in 1:q1)  IS1[j,1:j] <- 1
  IS1 <- matrix(1,q1,q1)  ## Define submatrix of Sigma
  IS1[upper.tri(IS1)] <- 0
  Sig <- IS1%x%IS # Knonecker product to get Sigma
  Sig[,1] <- rep(1,ncol(Sig))
  X <- X%*%Sig
  X # return the prediction matrix
}


########################################################################
## function used for predict method to get marginal model submatrices ## 
## with linear extrapolation if needed                                ##
########################################################################                                

marginal.linear.extrapolation <- function(object, data)
{ ## function to get marginal matrices used in predict method on bivariate SCOP-splines
  x <- data[[object$term[1]]]  
  z <- data[[object$term[2]]]  
  if (length(x) != length(z))
        stop ("arguments of smooth are not of the same dimension")
  m <- object$m + 1 ## vector of two components
  ## find spline basis inner knot range for 1st covariate, x...
  ll <- object$knots[[1]][m[1]+1];ul <- object$knots[[1]][length(object$knots[[1]])-m[1]]
  m[1] <- m[1] + 1
  n <- length(x)
  ind <- x<=ul & x>=ll ## data in range
  if (sum(ind)==n) { ## all in range
     X1 <- spline.des(object$knots[[1]],x,m[1])$design
  } else { ## some extrapolation needed 
     ## matrix mapping coefs to value and slope at end points...
     D <- spline.des(object$knots[[1]],c(ll,ll,ul,ul),m[1],c(0,1,0,1))$design
     X1 <- matrix(0,n,ncol(D)) ## full predict matrix
     if (sum(ind)> 0)  X1[ind,] <- spline.des(object$knots[[1]],x[ind],m[1])$design ## interior rows
     ## Now add rows for linear extrapolation...
     ind <- x < ll 
     if (sum(ind)>0) X1[ind,] <- cbind(1,x[ind]-ll)%*%D[1:2,]
     ind <- x > ul
     if (sum(ind)>0) X1[ind,] <- cbind(1,x[ind]-ul)%*%D[3:4,]
  }
  ## the same for 2nd maginal matrix...
  ## find spline basis inner knot range for 2nd covariate, z...
  ll <- object$knots[[2]][m[2]+1];ul <- object$knots[[2]][length(object$knots[[2]])-m[2]]
  m[2] <- m[2] + 1
  ind <- z<=ul & z>=ll ## data in range
  if (sum(ind)==n) { ## all in range
     X2 <- spline.des(object$knots[[2]],z,m[2])$design
  } else { ## some extrapolation needed 
     ## matrix mapping coefs to value and slope at end points...
     D <- spline.des(object$knots[[2]],c(ll,ll,ul,ul),m[2],c(0,1,0,1))$design
     X2 <- matrix(0,n,ncol(D)) ## full predict matrix
     if (sum(ind)> 0)  X2[ind,] <- spline.des(object$knots[[2]],z[ind],m[2])$design ## interior rows
     ## Now add rows for linear extrapolation...
     ind <- z < ll 
     if (sum(ind)>0) X2[ind,] <- cbind(1,z[ind]-ll)%*%D[1:2,]
     ind <- z > ul
     if (sum(ind)>0) X2[ind,] <- cbind(1,z[ind]-ul)%*%D[3:4,]
  }
  list(X1=X1, X2=X2)
}



####################################################################################### 
### Tensor product P-spline construction with double monotone increasing constraint  ##
#######################################################################################


smooth.construct.tedmi.smooth.spec <- function(object, data, knots)
{ ## construction of the double monotone increasing smooth surface
 # require(splines)
  if (object$dim !=2)
      stop("the number of covariates should be two")
  if (length(object$p.order)==1)
      m <- rep(object$p.order, 2) # if a single number is supplied the same
             ## order of P-splines is provided for both marginal smooths
  else m <- object$p.order
  m[is.na(m)] <- 2  # the default order is 2 (cubic P-spline)
  if (object$bs.dim[1]==-1)  # set the default values fro q1 and q2
     {  q1 <- object$bs.dim[1] <- 7
        q2 <- object$bs.dim[2] <- 7
     }
  else if (length(object$bs.dim)==1)
         q1 <- q2 <- object$bs.dim # if `k' is supplied as a single number, the same
             ## basis dimension is provided for both marginal smooths
  else  {q1 <- object$bs.dim[1]; q2 <- object$bs.dim[2]}
  if (is.na(q1)) q1 <- object$bs.dim[1] <- 7  # the default basis dimension is 7
  if (is.na(q2)) q2 <- object$bs.dim[2] <- 7
  nk1 <- q1+m[1]+2 ## number of knots for the 1st smooth
  nk2 <- q2+m[2]+2 ## number of knots for the 2nd smooth
  if (nk1<=0 || nk2<=0) stop("either k[1] or k[2] too small for m")
  ## the values of the first covariate...   
  x <- data[[object$term[1]]]  
  xk <- knots[[object$term[1]]] ## will be NULL if none supplied
  z <- data[[object$term[2]]]  ## the values of the second covariate
  zk <- knots[[object$term[2]]] ## will be NULL if none supplied
  if (is.null(xk)) # space knots through the values of the 1st covariate
     {  n<-length(x)
        xk<-rep(0,q1+m[1]+2)
        xk[(m[1]+2):(q1+1)]<-seq(min(x),max(x),length=q1-m[1])
        for (i in 1:(m[1]+1)) {xk[i]<-xk[m[1]+2]-(m[1]+2-i)*(xk[m[1]+3]-xk[m[1]+2])}
        for (i in (q1+2):(q1+m[1]+2)) {xk[i]<-xk[q1+1]+(i-q1-1)*(xk[m[1]+3]-xk[m[1]+2])}
     }
  if (n != length(z))
       stop ("arguments of smooth not same dimension")
  if (is.null(zk)) # space knots through the values of the 2nd covariate
      {  zk<-rep(0,q2+m[2]+2)
         zk[(m[2]+2):(q2+1)]<-seq(min(z),max(z),length=q2-m[2])
         for (i in 1:(m[2]+1)) {zk[i]<-zk[m[2]+2]-(m[2]+2-i)*(zk[m[2]+3]-zk[m[2]+2])}
         for (i in (q2+2):(q2+m[2]+2)) {zk[i]<-zk[q2+1]+(i-q2-1)*(zk[m[2]+3]-zk[m[2]+2])}
      }
  if (length(xk)!=nk1 || length(zk)!=nk2) # right number of knots?
         stop(paste("there should be ",nk1, " and ", nk2," supplied knots"))

  #  get model matrix-------------
  X1 <- splineDesign(xk,x,ord=m[1]+2)
  X2 <- splineDesign(zk,z,ord=m[2]+2)
  X <- matrix(0,n,q1*q2)  # model matrix
  for (i in 1:n)
      {  X[i,] <- X1[i,]%x%X2[i,] # Kronecker product of two rows of marginal model matrices
      }
  # get a matrix Sigma -----------------------
 # IS2 <- matrix(0,q2,q2)   # Define marginal matrix of Sigma
 # IS2[1:q2,1] <- rep(1,q2)
 #   for (j in 2:q2)  IS2[j,2:j] <- 1
 # IS1 <- matrix(0,q1,q1)   # Define marginal matrix of Sigma
 # IS1[1:q1,1] <- rep(1,q1)
 #   for (j in 2:q1)  IS1[j,2:j] <- 1
  IS1 <- matrix(1,q1,q1)  ## Define marginal matrix of Sigma
  IS1[upper.tri(IS1)] <- 0
  IS2 <- matrix(1,q2,q2)  ## Define marginal matrix of Sigma
  IS2[upper.tri(IS2)] <- 0
  Sig <- IS1 %x% IS2 
  
  # apply identifiability constraint and get model matrix
  X <- X[,2:ncol(X)]%*%Sig[2:ncol(Sig),2:ncol(Sig)]  
  object$X <- X # the finished model matrix with identifiability constraint
 
  # create the penalty matrix
  S <- list()
  I2<- diag(q2)
  I1 <- diff(diag(q1-1),difference=1) 
  Pm1 <- matrix(0,q1-1,q1) # marginal sqrt penalty
  Pm1[2:(q1-1),2:q1] <- I1
  S[[1]]<- Pm1%x%I2

  I2 <- diff(diag(q2-1),difference=1) 
  Pm2 <- matrix(0,q2-1,q2)
  Pm2[2:(q2-1),2:q2] <- I2  # marginal sqrt penalty
  I1 <- diag(q1)
  S[[2]] <- I1%x%Pm2
  
  object$P <- list()
  object$P[[1]] <- S[[1]][2:nrow(S[[1]]),2:ncol(S[[1]])]
  object$P[[2]] <- S[[2]][2:nrow(S[[2]]),2:ncol(S[[2]])]
  object$S <- list()
  object$S[[1]] <- crossprod(object$P[[1]]) ## t(object$P[[1]])%*%object$P[[1]]
  object$S[[2]] <- crossprod(object$P[[2]])  ## t(object$P[[2]])%*%object$P[[2]]
  object$p.ident <- rep(TRUE,q1*q2-1)   ## p.ident is an indicator of which coefficients must be positive (exponentiated)
  object$rank <- ncol(object$X)-1  # penalty rank
  object$null.space.dim <- m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  object$Zc <- diag(q1*q2-1) # identfiability constraint matrix
  object$Zc <- rbind(rep(0,ncol(object$Zc)),object$Zc)
  
  ## store "tedmi" specific stuff ...
  object$knots <- list()
  object$knots[[1]] <- xk
  object$knots[[2]] <- zk
  object$m <- m
   
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  class(object)<-"tedmi.smooth"  # Give object a class
  object
}


####################################################################


Predict.matrix.tedmi.smooth <- function(object, data)
{  ## prediction method function for the `tedmi' smooth class
  if (length(object$bs.dim)==1)
      q1 <- q2 <- object$bs.dim # if `k' is supplied as a single number, the same
             ## basis dimension is provided for both marginal smooths
  else  {q1 <- object$bs.dim[1]; q2 <- object$bs.dim[2]}
  
  bm <- marginal.linear.extrapolation(object, data)
  n <- length(data[[object$term[1]]])
  X <- matrix(0,n,q1*q2)  # model matrix
  for ( i in 1:n)
      {  X[i,] <- bm$X1[i,] %x% bm$X2[i,] # Kronecker product of two rows of marginal model matrices
      }
  # get a matrix Sigma -----------------------
  # IS2 <- matrix(0,q2,q2)   # Define marginal matrix of Sigma
 # IS2[1:q2,1] <- rep(1,q2)
 #   for (j in 2:q2)  IS2[j,2:j] <- 1
 # IS1 <- matrix(0,q1,q1)   # Define marginal matrix of Sigma
 # IS1[1:q1,1] <- rep(1,q1)
 #   for (j in 2:q1)  IS1[j,2:j] <- 1
  IS1 <- matrix(1,q1,q1)  ## Define marginal matrix of Sigma
  IS1[upper.tri(IS1)] <- 0
  IS2 <- matrix(1,q2,q2)  ## Define marginal matrix of Sigma
  IS2[upper.tri(IS2)] <- 0
  Sig <- IS1 %x% IS2 
  # get final model matrix
  X <- X %*%Sig
  X # return the prediction matrix
}



############################################################################ 
## Tensor product P-spline construction with single monotone decreasing   ##
## constraint wrt the first covariate                                     ##
############################################################################
                                                    
smooth.construct.tesmd1.smooth.spec<- function(object, data, knots)
{ ## construction of the single monotone decreasing smooth surface, deacreasing wrt the first covariate         
 # require(splines)
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
  S <- bm$S      
  
  # get the overall model matrix...
  X <- matrix(0,n,q1*q2)  # model matrix
  for (i in 1:n)
    {  X[i,] <- X1[i,]%x%X2[i,] # Kronecker product of two rows of marginal model matrices
    }
  # get a matrix Sigma -----------------------
 # IS <- matrix(0,q1,q1)   # Define marginal matrix of Sigma
 # IS[1:q1,1]<-1
 # for (j in 2:q1)  IS[j,2:j] <- -1
  IS <- matrix(-1,q1,q1)  ## coef summation matrix
  IS[upper.tri(IS)] <-0
  IS[,1] <- -IS[,1]
  I <- diag(q2)
  Sig <- IS%x%I

  # get model matrix
  X <- X%*%Sig
  # apply identifiability constraint 
  D <- diag(q1*q2)
  D <- D[,-q2]
  D1 <- t(diff(diag(q2)))
  D[1:q2,1:(q2-1)] <- D1
  X <- X%*%D 

  object$X <- X # the finished model matrix with identifiability constraint
  object$S <- list()
  object$S[[1]] <- crossprod(D,S[[1]])%*%D ##  t(D)%*%S[[1]]%*%D
  object$S[[2]] <- crossprod(D,S[[2]])%*%D ## t(D)%*%S[[2]]%*%D

  object$p.ident <- rep(TRUE,q1*q2-1)  
  object$p.ident[1:(q2-1)] <- rep(FALSE, q2-1) ## p.ident is an indicator of which coefficients must be positive (exponentiated)  
  object$rank <- ncol(object$X)-1  # penalty rank
  object$null.space.dim <- m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  object$Zc <- D   # identifiability constraint matrix

  ## store "tesmd1" specific stuff ...
  object$knots <- list()
  object$knots[[1]] <- xk
  if (is.null(zk))
     object$knots[[2]] <- rep(0,0,0)
  else object$knots[[2]] <- zk
  object$m <- m
  object$margin.bs <- bs2
      
  object$df <- ncol(object$X)     # maximum DoF (if unconstrained)
  class(object) <- "tesmd1.smooth"  # Give object a class
  object
}


###########################################################################

## Prediction matrix for the `tesmd1` smooth class *************************

Predict.matrix.tesmd1.smooth<-function(object,data)
{ ## prediction method function for the `tesmd1' smooth class
  if (length(object$bs.dim)==1)
      q1 <- q2 <- object$bs.dim # if `k' is supplied as a single number, the same
             ## basis dimension is provided for both marginal smooths
   else  {q1 <- object$bs.dim[1]; q2 <- object$bs.dim[2]}
  
  if (!is.null(object$xt)) bs2 <- object$xt ## basis for the marginal smooth along second direction
    else bs2 <- "ps"
  if (bs2=="cc") bm <- marginal.linear.extrapolation.tesmi1.cc(object, data)
    else bm <- marginal.linear.extrapolation(object, data)

  n <- length(data[[object$term[1]]])
  X <- matrix(0,n,q1*q2)  # model matrix
  for (i in 1:n)
    {  X[i,] <- bm$X1[i,] %x% bm$X2[i,] # Kronecker product of two rows of marginal model matrices
    }
  # get a matrix Sigma -----------------------
 # IS <- matrix(0,q1,q1)   # Define marginal matrix of Sigma
 # IS[1:q1,1]<-1
 # for (j in 2:q1)  IS[j,2:j] <- -1
  IS <- matrix(-1,q1,q1)  ## coef summation matrix
  IS[upper.tri(IS)] <-0
  IS[,1] <- -IS[,1]
  I <- diag(q2)
  Sig <- IS%x%I
  # get final model matrix
  X <- X%*%Sig
  X  # return the prediction matrix
}


############################################################################ 
## Tensor product P-spline construction with single monotone decreasing   ##
## constraint wrt the second covariate                                    ##
############################################################################

smooth.construct.tesmd2.smooth.spec<- function(object, data, knots)
## construction of the single monotone decreasing surface, decreasing wrt the second covariate        
{ 
  ## require(splines)
  #if (!is.null(object$xt)) bs2 <- object$xt ## basis for the marginal smooth along 1st direction
  #  else bs2 <- "ps" ## (only "ps" ia available currently)
  if (!is.null(object$xt)){
     if (!(object$xt %in% c("ps", "cc")) )
          stop("only 'ps' and 'cc' marginal basis are supported")
      else  bs2 <- object$xt ## basis for the marginal smooth along 1st direction
  } else bs2 <- "ps"

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
  if (object$bs.dim[1]==-1) { # set the default values fro q1 and q2
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
  
  nk2 <- q2+m[2]+2 ## number of knots for the 2nd smooth
  if (bs2=="cc") nk1 <- q1+1
    else nk1 <- q1+m[1]+2 ## number of knots for the 1st smooth in case of p-splines
  if (nk1<=0 || nk2<=0) stop("either k[1] or k[2] too small for m")
         
  ## the values of the first covariate...   
  x <- data[[object$term[1]]]  
  xk <- knots[[object$term[1]]] ## will be NULL if none supplied
  z <- data[[object$term[2]]]  ## the values of the second covariate
  zk <- knots[[object$term[2]]] ## will be NULL if none supplied
  if (is.null(xk)){ # space knots through the values of the 1st covariate
    if (bs2=="cc") {
         xk <- place.knots(x,nk1) 
         if (length(xk)==2) {
             xk <- place.knots(c(xk,x),nk1)
         }  
    } else{
         xk<-rep(0,q1+m[1]+2)
         xk[(m[1]+2):(q1+1)]<-seq(min(x),max(x),length=q1-m[1])
         for (i in 1:(m[1]+1)) {xk[i]<-xk[m[1]+2]-(m[1]+2-i)*(xk[m[1]+3]-xk[m[1]+2])}
         for (i in (q1+2):(q1+m[1]+2)) {xk[i]<-xk[q1+1]+(i-q1-1)*(xk[m[1]+3]-xk[m[1]+2])}
         knots[[object$term[1]]] <- xk
      }
   }
  n<-length(x)
  if (n != length(z))
     stop ("arguments of smooth not same dimension")
  if (is.null(zk)) # space knots through the values of the 2nd covariate
       { zk<-rep(0,q2+m[2]+2)
         zk[(m[2]+2):(q2+1)]<-seq(min(z),max(z),length=q2-m[2])
         for (i in 1:(m[2]+1)) {zk[i]<-zk[m[2]+2]-(m[2]+2-i)*(zk[m[2]+3]-zk[m[2]+2])}
         for (i in (q2+2):(q2+m[2]+2)) {zk[i]<-zk[q2+1]+(i-q2-1)*(zk[m[2]+3]-zk[m[2]+2])}
         knots[[object$term[2]]] <- zk
  }
  if (length(zk)!=nk2 ) # right number of knots?
      stop(paste("there should be ",nk2," supplied knots for z"))
  if (length(xk)!=nk1) # right number of knots?
      stop(paste("there should be ",nk1," supplied knots for x"))
  
  #  get model matrix-------------
 # get marginal model matrices and penalties... 
  if (bs2=="cc") bm <- marginal.matrices.tesmi2.cc(x,z,xk,zk,m,q1,q2)
    else bm <- marginal.matrices.tesmi2.ps(x,z,xk,zk,m,q1,q2)

  X1 <- bm$X1
  X2 <- bm$X2
  S <- bm$S  

  # get the overall model matrix...
  X <- matrix(0,n,q1*q2)  # model matrix
  for (i in 1:n)
    {  X[i,] <- X1[i,]%x%X2[i,] # Kronecker product of two rows of marginal model matrices
  }
  # get a matrix Sigma -----------------------
 # IS <- matrix(0,q2,q2)   # Define submatrix of Sigma
 # IS[1:q2,1]<-1
 # for (j in 2:q2)  IS[j,2:j] <- -1
  IS <- matrix(-1,q2,q2)  ## coef summation matrix
  IS[upper.tri(IS)] <-0
  IS[,1] <- -IS[,1]
  I <- diag(q1)
  Sig <- I%x%IS
    
  # get model matrix
  X <- X%*%Sig
  # apply identifiability constraint 
  D<- diag(q1*q2)
  D<-D[,-((q1-1)*q2+1)]
  ind <- rep(0,q1-1) # get index number for the cells to be changed to "-1"
  for (i in 1:(q1-1)){
        ind[i] <- (i-1)*q2+1
        D[ind[i],ind[i]] <- -1
  }
  for (i in 2:(q1-1))
        D[ind[i],ind[i]-q2] <- 1
  D[((q1-1)*q2+1),((q1-1)*q2+1-q2)] <- 1
  X <- X%*%D 

  object$X <- X # the finished model matrix with identifiability constraint
 
  # create the penalty matrix
  object$S <- list()
  object$S[[1]] <- t(D)%*%S[[1]]%*%D
  object$S[[2]] <- t(D)%*%S[[2]]%*%D

  object$p.ident <- rep(TRUE,q1*q2-1)  
  object$p.ident[ind] <- FALSE ## p.ident is an indicator of which coefficients must be positive (exponentiated)  
  object$rank <- ncol(object$X)-1  # penalty rank
  object$null.space.dim <- m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  object$Zc <- D   # identifiability constraint matrix

  ## store "tesmd2" specific stuff ...
  object$knots <- list()
  if (is.null(xk))
     object$knots[[1]] <- rep(0,0,0)
  else object$knots[[1]] <- xk
  object$knots[[2]] <- zk
  object$m <- m
  object$margin.bs <- bs2
   
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  class(object)<-"tesmd2.smooth"  # Give object a class
  object
}




####################################################################

marginal.matrices.tesmi2.ps <- function(x,z,xk,zk,m,q1,q2)
## function to get marginal model matrices and penalties in the overall 
## model coefficients in case of P-splines basis for the 1st unconstrained smooth
{
  # get marginal model matrix for the first unconstrained smooth... 
  X1 <- splineDesign(xk,x,ord=m[1]+2)
  # get marginal model matrix for the second monotonic smooth...
  X2 <- splineDesign(zk,z,ord=m[2]+2)

  # create the penalty matrix...
  S <- list()
  # get penalty matrix for the first unconstrained smooth...
  I2 <- diag(q2)
  P <- diff(diag(q1),difference=1)
  S[[1]]<- P %x% I2
  c1 <- c(1,-2,1)
  c2 <- c(1,rep(0,q2-1))
  c <- c1 %x% c2
  for (i in 1:(q1-2)){
       S[[1]][q2*(i-1)+1,(q2*(i-1)+1):(q2*(i-1)+length(c))] <- c
  }
  i <- q1-1
  S[[1]][q2*(i-1)+1,(q2*(i-1)+1):ncol(S[[1]])] <- rep(0,2*q2)

  S[[1]] <- crossprod(S[[1]])  ## t(S[[1]])%*%S[[1]]

  # get penalty for the 2nd monotonic smooth...
  I2 <- diff(diag(q2-1),difference=1) 
  P <- matrix(0,q2-1,q2)
  P[2:(q2-1),2:q2] <- I2  # marginal sqrt penalty
  I1 <- diag(q1)
  S[[2]] <- I1%x%P
  S[[2]] <- crossprod(S[[2]]) ## t(S[[2]])%*%S[[2]]

 list(X1=X1, X2=X2, S=S)
}


###########################################################################
## Prediction matrix for the `tesmd2` smooth class *************************


Predict.matrix.tesmd2.smooth<-function(object,data)
## prediction method function for the `tesmd2' smooth class
{ if (length(object$bs.dim)==1)
      q1 <- q2 <- object$bs.dim # if `k' is supplied as a single number, the same
             ## basis dimension is provided for both marginal smooths
  else  {q1 <- object$bs.dim[1]; q2 <- object$bs.dim[2]}
  
 # bm <- marginal.linear.extrapolation(object, data)
  if (!is.null(object$xt)) bs2 <- object$xt ## basis for the marginal smooth along second direction
    else bs2 <- "ps"
  if (bs2=="cc") bm <- marginal.linear.extrapolation.tesmi2.cc(object, data)
    else bm <- marginal.linear.extrapolation(object, data)

  n <- length(data[[object$term[1]]])
  X <- matrix(0,n,q1*q2)  # model matrix
  for (i in 1:n)
    {  X[i,] <- bm$X1[i,] %x% bm$X2[i,] # Kronecker product of two rows of marginal model matrices
  }
  # get a matrix Sigma -----------------------
 # IS <- matrix(0,q2,q2)   # Define submatrix of Sigma
 # IS[1:q2,1]<-1
 # for (j in 2:q2)  IS[j,2:j] <- -1
  IS <- matrix(-1,q2,q2)  ## coef summation matrix
  IS[upper.tri(IS)] <-0
  IS[,1] <- -IS[,1]
  I <- diag(q1)
  Sig <- I%x%IS
  # get final model matrix
  X <- X%*%Sig
  X  # return the prediction matrix
}



############################################################################ 
## Tensor product P-spline construction with single monotone increasing   ##
## constraint wrt the first covariate                                     ##
############################################################################

smooth.construct.tesmi1.smooth.spec<- function(object, data, knots)
## construction of the single monotone increasing surface, increasing wrt 1st covariate 
{ 
 # require(splines)
 # if (!is.null(object$xt)) bs2 <- object$xt ## basis for the marginal smooth along second direction
  #  else bs2 <- "ps"

 if (!is.null(object$xt)){
     if (!(object$xt %in% c("ps", "cc")) )
          stop("only 'ps' and 'cc' marginal basis are supported")
      else  bs2 <- object$xt ## basis for the marginal smooth along second direction
  } else bs2 <- "ps"

  ## if (bs2=="cc") print("tada") else print("nope")

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
  S <- bm$S      
  
  # get the overall model matrix...
  X <- matrix(0,n,q1*q2)  # model matrix
  for (i in 1:n)
    {  X[i,] <- X1[i,]%x%X2[i,] # Kronecker product of two rows of marginal model matrices
  }
  # get a matrix Sigma -----------------------
 # IS <- matrix(0,q1,q1)   # Define marginal matrix of Sigma
 # IS[1:q1,1]<-1
 # for (j in 2:q1)  IS[j,2:j] <- 1
  IS <- matrix(1,q1,q1)  ## coef summation matrix
  IS[upper.tri(IS)] <-0
  I <- diag(q2)
  Sig <- IS%x%I

  # get model matrix
  X <- X%*%Sig
  # apply identifiability constraint 
  D <- diag(q1*q2)
  D <- D[,-q2]
  D1 <- t(diff(diag(q2)))
  D[1:q2,1:(q2-1)] <- D1
  X <- X%*%D 

  object$X <- X # the finished model matrix with identifiability constraint
 
  object$S <- list()
  object$S[[1]] <- crossprod(D,S[[1]])%*%D ## t(D)%*%S[[1]]%*%D
  object$S[[2]] <- crossprod(D,S[[2]])%*%D ## t(D)%*%S[[2]]%*%D

  object$p.ident <- rep(TRUE,q1*q2-1)  
  object$p.ident[1:(q2-1)] <- rep(FALSE, q2-1) ## p.ident is an indicator of which coefficients must be positive (exponentiated)  
  object$rank <- ncol(object$X)-1  # penalty rank
  object$null.space.dim <- m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  object$Zc <- D   # identifiability constraint matrix

  ## store "tesmi1" specific stuff ...
  object$knots <- list()
  object$knots[[1]] <- xk
  if (is.null(zk))
     object$knots[[2]] <- rep(0,0,0)
  else object$knots[[2]] <- zk
  object$m <- m
  object$margin.bs <- bs2
       
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  class(object)<-"tesmi1.smooth"  # Give object a class
  object
}


###############################################################################################
## Cyclic cubic regression spline methods copied from mgcv:: smooth.construct.cc.smooth.spec())
################################################################################################
place.knots <- function(x,nk)
# knot placement code. x is a covariate array, nk is the number of knots,
# and this routine spaces nk knots evenly throughout the x values, with the 
# endpoints at the extremes of the data.
{ x<-sort(unique(x));n<-length(x)
  if (nk>n) stop("more knots than unique data values is not allowed")
  if (nk<2) stop("too few knots")
  if (nk==2) return(range(x))
  delta<-(n-1)/(nk-1) # how many data steps per knot
  lbi<-floor(delta*1:(nk-2))+1 # lower interval bound index
  frac<-delta*1:(nk-2)+1-lbi # left over proportion of interval  
  x.shift<-x[-1]
  knot<-array(0,nk)
  knot[nk]<-x[n];knot[1]<-x[1]
  knot[2:(nk-1)]<-x[lbi]*(1-frac)+x.shift[lbi]*frac
  knot
} ## place.knots

cwrap <- function(x0,x1,x) {
## map x onto [x0,x1] in manner suitable for cyclic smooth on
## [x0,x1].
  h <- x1-x0
  if (max(x)>x1) {
    ind <- x>x1
    x[ind] <- x0 + (x[ind]-x1)%%h
  }
  if (min(x)<x0) {
    ind <- x<x0
    x[ind] <- x1 - (x0-x[ind])%%h
  }
  x
} ## cwrap

pred.mat<-function(x,knots,BD)
  # BD is B^{-1}D. Basis as given in Lancaster and Salkauskas (1986)
  # Curve and Surface fitting, but wrapped to give periodic smooth.
  { j<-x
    n<-length(knots)
    h<-knots[2:n]-knots[1:(n-1)]
    if (max(x)>max(knots)||min(x)<min(knots)) x <- cwrap(min(knots),max(knots),x)
    ## stop("can't predict outside range of knots with periodic smoother")
    for (i in n:2) j[x<=knots[i]]<-i
    j1<-hj<-j-1
    j[j==n]<-1
    I<-diag(n-1)
    X<-BD[j1,,drop=FALSE]*as.numeric(knots[j1+1]-x)^3/as.numeric(6*h[hj])+
       BD[j,,drop=FALSE]*as.numeric(x-knots[j1])^3/as.numeric(6*h[hj])-
       BD[j1,,drop=FALSE]*as.numeric(h[hj]*(knots[j1+1]-x)/6)-
       BD[j,,drop=FALSE]*as.numeric(h[hj]*(x-knots[j1])/6) +
       I[j1,,drop=FALSE]*as.numeric((knots[j1+1]-x)/h[hj]) +
       I[j,,drop=FALSE]*as.numeric((x-knots[j1])/h[hj])
    X
} # end of pred.mat local function

getBD <- function(x){
  # matrices B and D in expression Bm=Dp where m are s"(x_i) and 
  # p are s(x_i) and the x_i are knots of periodic spline s(x)
  # B and D slightly modified (for periodicity) from Lancaster 
  # and Salkauskas (1986) Curve and Surface Fitting section 4.7.
    n<-length(x)
    h<-x[2:n]-x[1:(n-1)]
    n<-n-1
    D<-B<-matrix(0,n,n)
    B[1,1]<-(h[n]+h[1])/3;B[1,2]<-h[1]/6;B[1,n]<-h[n]/6
    D[1,1]<- -(1/h[1]+1/h[n]);D[1,2]<-1/h[1];D[1,n]<-1/h[n]
    for (i in 2:(n-1))
    { B[i,i-1]<-h[i-1]/6
      B[i,i]<-(h[i-1]+h[i])/3
      B[i,i+1]<-h[i]/6
      D[i,i-1]<-1/h[i-1]
      D[i,i]<- -(1/h[i-1]+1/h[i])
      D[i,i+1]<- 1/h[i]
    }
    B[n,n-1]<-h[n-1]/6;B[n,n]<-(h[n-1]+h[n])/3;B[n,1]<-h[n]/6
    D[n,n-1]<-1/h[n-1];D[n,n]<- -(1/h[n-1]+1/h[n]);D[n,1]<-1/h[n]
    list(B=B,D=D)
  } # end of getBD local function


marginal.matrices.tesmi1.cc <- function(x,z,xk,zk,m,q1,q2)
## function to get marginal model matrices and penalties in the overall 
## model coefficients in case of cyclic cubic splines for the 2nd marginal smooth
{
  # get marginal model matrix for the first monotonic smooth... 
  X1 <- splineDesign(xk,x,ord=m[1]+2)

  # get marginal model matrix for the second unconstrained smooth
  ## (copied from mgcv:: smooth.construct.cc.smooth.spec())...
   
  um <- getBD(zk)
  BD <- solve(um$B,um$D) # s"(k)=BD%*%s(k) where k are knots minus last knot
   
  
  X2 <- pred.mat(z,zk,BD)
   ## object$S<-list(t(um$D)%*%BD) # the penalty
   ##  object$S[[1]]<-(object$S[[1]]+t(object$S[[1]]))/2 # ensure exact symmetry
  ##============================

  # create the penalty matrices...
  S <- list()
  # get the penalty matrix for the first monotone smooth...
  I2<- diag(q2)
  P <- diff(diag(q1-1),difference=1) 
  Pm1 <- matrix(0,q1-1,q1) # marginal sqrt penalty
  Pm1[2:(q1-1),2:q1] <- P
  S[[1]]<- Pm1%x%I2
  S[[1]] <- crossprod(S[[1]])  ## t(S[[1]])%*%S[[1]]

  # get penalty for the 2nd 'cc' smooth...
  I2 <- diff(diag(q2),difference=1) 
  I21<- diff(diag(q2),difference=2)
  I1 <- diag(q1)
  S[[2]] <-matrix(0,q2-2+(q1-1)*(q2-1), q1*q2)
  S[[2]][1:(q2-2),] <- t(I1[1,])%x%I21
  S[[2]][(q2-1):nrow(S[[2]]),] <- I1[2:q1,]%x%I2
  S[[2]] <- crossprod(S[[2]])  ## t(S[[2]])%*%S[[2]]
  S[[2]] <- (S[[2]]+t(S[[2]]))/2 # ensure exact symmetry
 list(X1=X1, X2=X2, S=S)
}


####################################################################

marginal.matrices.tesmi1.ps <- function(x,z,xk,zk,m,q1,q2)
## function to get marginal model matrices and penalties in the overall 
## model coefficients in case of P-splines basis for the 2nd unconstrained marginal smooth
{
  # get marginal model matrix for the first monotonic smooth... 
  X1 <- splineDesign(xk,x,ord=m[1]+2)
  # get marginal model matrix for the second unconstrained smooth...
  X2 <- splineDesign(zk,z,ord=m[2]+2)
   # create the penalty matrix...
  S <- list()
  # get the penalty matrix for the first monotone smooth...
  I2<- diag(q2)
  P <- diff(diag(q1-1),difference=1) 
  Pm1 <- matrix(0,q1-1,q1) # marginal sqrt penalty
  Pm1[2:(q1-1),2:q1] <- P
  S[[1]]<- Pm1%x%I2
  S[[1]] <- crossprod(S[[1]])  ## t(S[[1]])%*%S[[1]]

  # get penalty for the 2nd smooth
  I2 <- diff(diag(q2),difference=1) 
  I21<- diff(diag(q2),difference=2)
  I1 <- diag(q1)
  S[[2]] <-matrix(0,q2-2+(q1-1)*(q2-1), q1*q2)
  S[[2]][1:(q2-2),] <- t(I1[1,])%x%I21
  S[[2]][(q2-1):nrow(S[[2]]),] <- I1[2:q1,]%x%I2
  S[[2]] <- crossprod(S[[2]])  ## t(S[[2]])%*%S[[2]]
 list(X1=X1, X2=X2, S=S)
}


## Prediction matrix for the `tesmi1` smooth class *************************

Predict.matrix.tesmi1.smooth<-function(object,data)
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
  X <- matrix(0,n,q1*q2)  # model matrix
  for (i in 1:n)
    {  X[i,] <- bm$X1[i,] %x% bm$X2[i,] # Kronecker product of two rows of marginal model matrices
  }
  # get a matrix Sigma -----------------------
 # IS <- matrix(0,q1,q1)   # Define marginal matrix of Sigma
 # IS[1:q1,1]<-1
 # for (j in 2:q1)  IS[j,2:j] <- 1
  IS <- matrix(1,q1,q1)  ## coef summation matrix
  IS[upper.tri(IS)] <-0
  I <- diag(q2)
  Sig <- IS%x%I

  # get final model matrix
  X <- X%*%Sig
  X  # return the prediction matrix
}


## function used for predict method to get marginal model submatrices for tesmi1.cc
## with linear extrapolation along x1 if needed                                          

marginal.linear.extrapolation.tesmi1.cc <- function(object, data)
{ ## function to get marginal matrices used in predict method on bivariate SCOP-splines
  x <- data[[object$term[1]]]  
  z <- data[[object$term[2]]]  
  if (length(x) != length(z))
        stop ("arguments of smooth are not of the same dimension")
  m <- object$m + 1 ## vector of two components

  ## find spline basis inner knot range for 1st covariate, x...
  ll <- object$knots[[1]][m[1]+1];ul <- object$knots[[1]][length(object$knots[[1]])-m[1]]
  m[1] <- m[1] + 1
  n <- length(x)
  ind <- x<=ul & x>=ll ## data in range
  if (sum(ind)==n) { ## all in range
     X1 <- spline.des(object$knots[[1]],x,m[1])$design
  } else { ## some extrapolation needed 
     ## matrix mapping coefs to value and slope at end points...
     D <- spline.des(object$knots[[1]],c(ll,ll,ul,ul),m[1],c(0,1,0,1))$design
     X1 <- matrix(0,n,ncol(D)) ## full predict matrix
     if (sum(ind)> 0)  X1[ind,] <- spline.des(object$knots[[1]],x[ind],m[1])$design ## interior rows
     ## Now add rows for linear extrapolation...
     ind <- x < ll 
     if (sum(ind)>0) X1[ind,] <- cbind(1,x[ind]-ll)%*%D[1:2,]
     ind <- x > ul
     if (sum(ind)>0) X1[ind,] <- cbind(1,x[ind]-ul)%*%D[3:4,]
  }
 
  ## "cc" basis for the 2nd marginal matrix...
  um <- getBD(object$knots[[2]])
  BD <- solve(um$B,um$D) # s"(k)=BD%*%s(k) where k are knots minus last knot

  X2 <- pred.mat(z,object$knots[[2]],BD)
  list(X1=X1, X2=X2)
}


############################################################################ 
## Tensor product P-spline construction with single monotone increasing   ##
## constraint wrt the second covariate                                    ##
############################################################################

smooth.construct.tesmi2.smooth.spec<- function(object, data, knots)
## construction of the single monotone increasing smooth surface, increasing wrt the second covariate 
{ 
 ## require(splines)
 # if (!is.null(object$xt)) bs2 <- object$xt ## basis for the marginal smooth along 1st direction
  #  else bs2 <- "ps" ## (only "ps" is available currently)

 if (!is.null(object$xt)){
     if (!(object$xt %in% c("ps", "cc")) )
          stop("only 'ps' and 'cc' marginal basis are supported")
      else  bs2 <- object$xt ## basis for the marginal smooth along 1st direction
  } else bs2 <- "ps"
  
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
  
  nk2 <- q2+m[2]+2 ## number of knots for the 2nd SCOP-smooth
  ## nk1 <- q1+m[1]+2 ## number of knots for the 1st smooth
  if (bs2=="cc") nk1 <- q1+1
    else nk1 <- q1+m[1]+2 ## number of knots for the 1st smooth in case of p-splines

  if (nk1<=0 || nk2<=0) stop("either k[1] or k[2] too small for m")
          
  ## the values of the first covariate...   
  x <- data[[object$term[1]]]  
  xk <- knots[[object$term[1]]] ## will be NULL if none supplied
  z <- data[[object$term[2]]]  ## the values of the second covariate
  zk <- knots[[object$term[2]]] ## will be NULL if none supplied
  n<-length(x)
  if (is.null(xk)){ # space knots through the values of the 1st covariate
    if (bs2=="cc") {
         xk <- place.knots(x,nk1) 
         if (length(xk)==2) {
             xk <- place.knots(c(xk,x),nk1)
         }  
    } else{
         xk<-rep(0,q1+m[1]+2)
         xk[(m[1]+2):(q1+1)]<-seq(min(x),max(x),length=q1-m[1])
         for (i in 1:(m[1]+1)) {xk[i]<-xk[m[1]+2]-(m[1]+2-i)*(xk[m[1]+3]-xk[m[1]+2])}
         for (i in (q1+2):(q1+m[1]+2)) {xk[i]<-xk[q1+1]+(i-q1-1)*(xk[m[1]+3]-xk[m[1]+2])}
         knots[[object$term[1]]] <- xk
      }
   }

  if (n != length(z))
     stop ("arguments of smooth not same dimension")
  if (is.null(zk)) # space knots through the values of the 2nd covariate
       { zk<-rep(0,q2+m[2]+2)
         zk[(m[2]+2):(q2+1)]<-seq(min(z),max(z),length=q2-m[2])
         for (i in 1:(m[2]+1)) {zk[i]<-zk[m[2]+2]-(m[2]+2-i)*(zk[m[2]+3]-zk[m[2]+2])}
         for (i in (q2+2):(q2+m[2]+2)) {zk[i]<-zk[q2+1]+(i-q2-1)*(zk[m[2]+3]-zk[m[2]+2])}
         knots[[object$term[2]]] <- zk
  }
  if (length(zk)!=nk2 ) # right number of knots?
      stop(paste("there should be ",nk2," supplied knots for z"))
  if (length(xk)!=nk1) # right number of knots?
      stop(paste("there should be ",nk1," supplied knots for x"))
  
  #  get model matrix-------------
  # get marginal model matrices and penalties... 
  if (bs2=="cc") bm <- marginal.matrices.tesmi2.cc(x,z,xk,zk,m,q1,q2)
    else bm <- marginal.matrices.tesmi2.ps(x,z,xk,zk,m,q1,q2)

  X1 <- bm$X1
  X2 <- bm$X2
  S <- bm$S  
  # get the overall model matrix...
  X <- matrix(0,n,q1*q2)  # model matrix
  for (i in 1:n)
    {  X[i,] <- X1[i,]%x%X2[i,] # Kronecker product of two rows of marginal model matrices
    }
  # get a matrix Sigma -----------------------  
 # IS <- matrix(0,q2,q2)   # Define submatrix of Sigma
 # IS[1:q2,1]<-1
 # for (j in 2:q2)  IS[j,2:j] <- 1
  IS <- matrix(1,q2,q2)  ## coef summation matrix
  IS[upper.tri(IS)] <-0
  I <- diag(q1)
  Sig <- I%x%IS
    
  # get model matrix
  X <- X%*%Sig
  # apply identifiability constraint 
  D<- diag(q1*q2)
  D<-D[,-((q1-1)*q2+1)]
  ind <- rep(0,q1-1) # get index number for the cells to be changed to "-1"
  for (i in 1:(q1-1)){
        ind[i] <- (i-1)*q2+1
        D[ind[i],ind[i]] <- -1
      }
  for (i in 2:(q1-1))
        D[ind[i],ind[i]-q2] <- 1
  D[((q1-1)*q2+1),((q1-1)*q2+1-q2)] <- 1
  X <- X%*%D 

  object$X <- X # the finished model matrix with identifiability constraint
 
  # create the penalty matrix
  object$S <- list()
  object$S[[1]] <- t(D)%*%S[[1]]%*%D
  object$S[[2]] <- t(D)%*%S[[2]]%*%D

  object$p.ident <- rep(TRUE,q1*q2-1)  
  object$p.ident[ind]<-FALSE ## p.ident is an indicator of which coefficients must be positive (exponentiated)
  object$rank <- ncol(object$X)-1  # penalty rank
  object$null.space.dim <- m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  object$Zc <- D   # identifiability constraint matrix

  ## store "tesmi2" specific stuff ...
  object$knots <- list()
  if (is.null(xk))
     object$knots[[1]] <- rep(0,0,0)
  else object$knots[[1]] <- xk
  object$knots[[2]] <- zk
  object$m <- m
  object$margin.bs <- bs2
   
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  class(object)<-"tesmi2.smooth"  # Give object a class
  object
}


####################################################################

marginal.matrices.tesmi2.ps <- function(x,z,xk,zk,m,q1,q2)
## function to get marginal model matrices and penalties in the overall 
## model coefficients in case of P-splines basis for the 1st unconstrained smooth
{
  # get marginal model matrix for the first unconstrained smooth... 
  X1 <- splineDesign(xk,x,ord=m[1]+2)
  # get marginal model matrix for the second monotonic smooth...
  X2 <- splineDesign(zk,z,ord=m[2]+2)

  # create the penalty matrix...
  S <- list()
  # get penalty matrix for the first unconstrained smooth...
  I2 <- diag(q2)
  P <- diff(diag(q1),difference=1)
  S[[1]]<- P %x% I2
  c1 <- c(1,-2,1)
  c2 <- c(1,rep(0,q2-1))
  c <- c1 %x% c2
  for (i in 1:(q1-2)){
       S[[1]][q2*(i-1)+1,(q2*(i-1)+1):(q2*(i-1)+length(c))] <- c
  }
  i <- q1-1
  S[[1]][q2*(i-1)+1,(q2*(i-1)+1):ncol(S[[1]])] <- rep(0,2*q2)
  S[[1]] <- crossprod(S[[1]])  ## t(S[[1]])%*%S[[1]]

  # get penalty for the 2nd monotonic smooth...
  I2 <- diff(diag(q2-1),difference=1) 
  P <- matrix(0,q2-1,q2)
  P[2:(q2-1),2:q2] <- I2  # marginal sqrt penalty
  I1 <- diag(q1)
  S[[2]] <- I1%x%P
  S[[2]] <- crossprod(S[[2]])  ## t(S[[2]])%*%S[[2]]

 list(X1=X1, X2=X2, S=S)
}


marginal.matrices.tesmi2.cc <- function(x,z,xk,zk,m,q1,q2)
## function to get marginal model matrices and penalties in the overall 
## model coefficients in case of cyclic cubic splines for the 1st marginal smooth
{
  # get marginal model matrix for the 2nd monotonic smooth... 
  X2 <- splineDesign(zk,z,ord=m[2]+2)

  # get marginal model matrix for the 1st 'cc' smooth
  ## (copied from mgcv:: smooth.construct.cc.smooth.spec())...
  um <- getBD(xk)
  BD <- solve(um$B,um$D) # s"(k)=BD%*%s(k) where k are knots minus last knot
  
  X1 <- pred.mat(x,xk,BD)
   ## object$S<-list(t(um$D)%*%BD) # the penalty
   ##  object$S[[1]]<-(object$S[[1]]+t(object$S[[1]]))/2 # ensure exact symmetry
  ##============================

  # create the penalty matrices...
  S <- list()
  # get penalty matrix for the 1st 'cc' smooth...
  I2 <- diag(q2)
  P <- diff(diag(q1),difference=1)
  S[[1]]<- P %x% I2
  c1 <- c(1,-2,1)
  c2 <- c(1,rep(0,q2-1))
  c <- c1 %x% c2
  for (i in 1:(q1-2)){
       S[[1]][q2*(i-1)+1,(q2*(i-1)+1):(q2*(i-1)+length(c))] <- c
  }
  i <- q1-1
  S[[1]][q2*(i-1)+1,(q2*(i-1)+1):ncol(S[[1]])] <- rep(0,2*q2)
  S[[1]] <- crossprod(S[[1]])  ## t(S[[1]])%*%S[[1]]
  S[[1]] <- (S[[1]]+t(S[[1]]))/2 # ensure exact symmetry

  # get penalty for the 2nd monotonic smooth...
  I2 <- diff(diag(q2-1),difference=1) 
  P <- matrix(0,q2-1,q2)
  P[2:(q2-1),2:q2] <- I2  # marginal sqrt penalty
  I1 <- diag(q1)
  S[[2]] <- I1%x%P
  S[[2]] <- crossprod(S[[2]])  ## t(S[[2]])%*%S[[2]]

 list(X1=X1, X2=X2, S=S)
}


############################################################
## Prediction matrix for the `tesmi2` smooth class ....   ##
############################################################

Predict.matrix.tesmi2.smooth<-function(object,data)
## prediction method function for the `tesmi2' smooth class
{ if (length(object$bs.dim)==1)
      q1 <- q2 <- object$bs.dim # if `k' is supplied as a single number, the same
             ## basis dimension is provided for both marginal smooths
  else  {q1 <- object$bs.dim[1]; q2 <- object$bs.dim[2]}
  
  if (!is.null(object$xt)) bs2 <- object$xt ## basis for the marginal smooth along second direction
    else bs2 <- "ps"
  if (bs2=="cc") bm <- marginal.linear.extrapolation.tesmi2.cc(object, data)
    else bm <- marginal.linear.extrapolation(object, data)

  n <- length(data[[object$term[1]]])
  X <- matrix(0,n,q1*q2)  # model matrix
  for (i in 1:n)
    {  X[i,] <- bm$X1[i,] %x% bm$X2[i,] # Kronecker product of two rows of marginal model matrices
    }
  # get a matrix Sigma -----------------------
 # IS <- matrix(0,q2,q2)   # Define submatrix of Sigma
 # IS[1:q2,1]<-1
 # for (j in 2:q2)  IS[j,2:j] <- 1
  IS <- matrix(1,q2,q2)  ## coef summation matrix
  IS[upper.tri(IS)] <-0
  I <- diag(q1)
  Sig <- I%x%IS
    
  # get final model matrix
  X <- X%*%Sig
  X  # return the prediction matrix
}


## function used for predict method to get marginal model submatrices for tesmi2.cc and tesmd2.cc (along x1)
## with linear extrapolation along x2 (SCOP-spline marginal) if needed...                                           

marginal.linear.extrapolation.tesmi2.cc <- function(object, data)
{ ## function to get marginal matrices used in predict method on bivariate SCOP-splines
  x <- data[[object$term[1]]]  
  z <- data[[object$term[2]]]  
  if (length(x) != length(z))
        stop ("arguments of smooth are not of the same dimension")
  m <- object$m + 1 ## vector of two components

  ## find spline basis inner knot range for 2nd covariate, z...
  ll <- object$knots[[2]][m[2]+1];ul <- object$knots[[2]][length(object$knots[[2]])-m[2]]
  m[2] <- m[2] + 1
  n <- length(z)
  ind <- z<=ul & z>=ll ## data in range
  if (sum(ind)==n) { ## all in range
     X2 <- spline.des(object$knots[[2]],z,m[2])$design
  } else { ## some extrapolation needed 
     ## matrix mapping coefs to value and slope at end points...
     D <- spline.des(object$knots[[2]],c(ll,ll,ul,ul),m[2],c(0,1,0,1))$design
     X2 <- matrix(0,n,ncol(D)) ## full predict matrix
     if (sum(ind)> 0)  X2[ind,] <- spline.des(object$knots[[2]],z[ind],m[2])$design ## interior rows
     ## Now add rows for linear extrapolation...
     ind <- z < ll 
     if (sum(ind)>0) X2[ind,] <- cbind(1,z[ind]-ll)%*%D[1:2,]
     ind <- z > ul
     if (sum(ind)>0) X2[ind,] <- cbind(1,z[ind]-ul)%*%D[3:4,]
  }
 
  ## "cc" basis for the 1st marginal matrix...
  um <- getBD(object$knots[[1]])
  BD <- solve(um$B,um$D) # s"(k)=BD%*%s(k) where k are knots minus last knot

  X1 <- pred.mat(x,object$knots[[1]],BD)
  list(X1=X1, X2=X2)
}

###############################################################
## Tensor product P-spline construction with mixed constraints:
## increasing wrt the 1st covariate and convex wrt the 2nd covariate ...
###################################


smooth.construct.temicx.smooth.spec<- function(object, data, knots)
## construction of the bivariate smooth surface with mixed constraints: increasing
## wrt the 1st covariate and convex wrt the 2nd one...
{ 
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
  nk2 <- q2+m[2]+2 ## number of knots for the 2nd smooth
  if (nk1<=0 || nk2<=0) stop("either k[1] or k[2] too small for m")
  
  ## the values of the first covariate...   
  x <- data[[object$term[1]]]  
  xk <- knots[[object$term[1]]] ## will be NULL if none supplied
  z <- data[[object$term[2]]]  ## the values of the second covariate
  zk <- knots[[object$term[2]]] ## will be NULL if none supplied
  if (is.null(xk)) # space knots through the values of the 1st covariate
  { n<-length(x)
    xk<-rep(0,q1+m[1]+2)
    xk[(m[1]+2):(q1+1)]<-seq(min(x),max(x),length=q1-m[1])
    for (i in 1:(m[1]+1)) {xk[i]<-xk[m[1]+2]-(m[1]+2-i)*(xk[m[1]+3]-xk[m[1]+2])}
    for (i in (q1+2):(q1+m[1]+2)) {xk[i]<-xk[q1+1]+(i-q1-1)*(xk[m[1]+3]-xk[m[1]+2])}
    knots[[object$term[1]]] <- xk
  }
  if (n != length(z))
     stop ("arguments of smooth not same dimension")
  if (is.null(zk)) # space knots through the values of the 2nd covariate
  {      zk<-rep(0,q2+m[2]+2)
         zk[(m[2]+2):(q2+1)]<-seq(min(z),max(z),length=q2-m[2])
         for (i in 1:(m[2]+1)) {zk[i]<-zk[m[2]+2]-(m[2]+2-i)*(zk[m[2]+3]-zk[m[2]+2])}
         for (i in (q2+2):(q2+m[2]+2)) {zk[i]<-zk[q2+1]+(i-q2-1)*(zk[m[2]+3]-zk[m[2]+2])}
         knots[[object$term[2]]] <- zk
   }
  if (length(xk)!=nk1 ) # right number of knots?
      stop(paste("there should be ",nk1," supplied knotsfor the x"))
  if (length(zk)!=nk2) # right number of knots?
      stop(paste("there should be ",nk2," supplied knots for z"))
  
 #  get model matrix-------------
  X1 <- splineDesign(xk,x,ord=m[1]+2)
  X2 <- splineDesign(zk,z,ord=m[2]+2)
  X <- matrix(0,n,q1*q2)  # model matrix
  for (i in 1:n)
      {  X[i,] <- X1[i,]%x%X2[i,] # Kronecker product of two rows of marginal model matrices
      }
  # get a matrix Sigma -----------------------
  IS2 <- matrix(0,q2,q2)   # Define marginal matrix of Sigma for convexity
  IS2[1:q2,1] <- rep(1,q2)
  IS2[2:q2,2]<- -c(1:(q2-1))
  for (i in 3:q2) IS2[i:q2,i] <- c(1:(q2-i+1))  

#  IS1 <- matrix(0,q1,q1)   # Define marginal matrix of Sigma for increasing constraint
#  IS1[1:q1,1] <- rep(1,q1)
#    for (j in 2:q1)  IS1[j,2:j] <- 1
  IS1 <- matrix(1,q1,q1)  ## Define marginal matrix of Sigma for increasing constraint
  IS1[upper.tri(IS1)] <-0

  Sig <- IS1 %x% IS2 
  
  # apply identifiability constraint and get model matrix
  X <- X[,2:ncol(X)]%*%Sig[2:ncol(Sig),2:ncol(Sig)]  
  object$X <- X # the finished model matrix with identifiability constraint
 
  # create the penalty matrix
  S <- list()
  I2<- diag(q2)
  I1 <- diff(diag(q1-1),difference=1) 
  Pm1 <- matrix(0,q1-1,q1) # marginal sqrt penalty
  Pm1[2:(q1-1),2:q1] <- I1
  S[[1]]<- Pm1%x%I2

  I2 <- diff(diag(q2-2),difference=1) 
  Pm2 <- matrix(0,q2-1,q2)
  Pm2[3:(q2-1),3:q2] <- I2  # marginal sqrt penalty
  I1 <- diag(q1)
  S[[2]] <- I1%x%Pm2
  
  object$P <- list()
  object$P[[1]] <- S[[1]][2:nrow(S[[1]]),2:ncol(S[[1]])]
  object$P[[2]] <- S[[2]][2:nrow(S[[2]]),2:ncol(S[[2]])]
  object$S <- list()
  object$S[[1]] <- crossprod(object$P[[1]]) ## t(object$P[[1]])%*%object$P[[1]]
  object$S[[2]] <- crossprod(object$P[[2]])  ## t(object$P[[2]])%*%object$P[[2]]

  object$p.ident <- rep(TRUE,q1*q2-1)  ## p.ident is an indicator of which coefficients must be positive (exponentiated)
  object$rank <- ncol(object$X)-1  # penalty rank
  object$null.space.dim <- m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  object$Zc <- diag(q1*q2-1) # identfiability constraint matrix
  object$Zc <- rbind(rep(0,ncol(object$Zc)),object$Zc)
  
  ## store "temicx" specific stuff ...
  object$knots <- list()
  object$knots[[1]] <- xk
  object$knots[[2]] <- zk
  object$m <- m
   
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  class(object)<-"temicx.smooth"  # Give object a class
  object
}

####################################################################


Predict.matrix.temicx.smooth <- function(object, data)
{  ## prediction method function for the `temicx' smooth class
  if (length(object$bs.dim)==1)
      q1 <- q2 <- object$bs.dim # if `k' is supplied as a single number, the same
             ## basis dimension is provided for both marginal smooths
  else  {q1 <- object$bs.dim[1]; q2 <- object$bs.dim[2]}
  
  bm <- marginal.linear.extrapolation(object, data)
  n <- length(data[[object$term[1]]])
  X <- matrix(0,n,q1*q2)  # model matrix
  for ( i in 1:n)
      {  X[i,] <- bm$X1[i,] %x% bm$X2[i,] # Kronecker product of two rows of marginal model matrices
      }
  # get a matrix Sigma -----------------------
  IS2 <- matrix(0,q2,q2)   # Define marginal matrix of Sigma for convexity
  IS2[1:q2,1] <- rep(1,q2)
  IS2[2:q2,2]<- -c(1:(q2-1))
  for (i in 3:q2) IS2[i:q2,i] <- c(1:(q2-i+1))  

 #  IS1 <- matrix(0,q1,q1)   # Define marginal matrix of Sigma for increasing constraint
 #  IS1[1:q1,1] <- rep(1,q1)
 #    for (j in 2:q1)  IS1[j,2:j] <- 1
  IS1 <- matrix(1,q1,q1)  ## Define marginal matrix of Sigma for increasing constraint
  IS1[upper.tri(IS1)] <-0

  Sig <- IS1 %x% IS2 
  # get final model matrix
  X <- X %*%Sig
  X # return the prediction matrix
}



################################################################
## Tensor product P-spline construction with mixed constraints:
## increasing wrt the 1st covariate and concave wrt the 2nd covariate ...
#################################################################


smooth.construct.temicv.smooth.spec<- function(object, data, knots)
## construction of the bivariate smooth surface with mixed constraints: increasing
## wrt the 1st covariate and concave wrt the 2nd one...
{ 
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
  nk2 <- q2+m[2]+2 ## number of knots for the 2nd smooth
  if (nk1<=0 || nk2<=0) stop("either k[1] or k[2] too small for m")
  
  ## the values of the first covariate...   
  x <- data[[object$term[1]]]  
  xk <- knots[[object$term[1]]] ## will be NULL if none supplied
  z <- data[[object$term[2]]]  ## the values of the second covariate
  zk <- knots[[object$term[2]]] ## will be NULL if none supplied
  if (is.null(xk)) # space knots through the values of the 1st covariate
  { n<-length(x)
    xk<-rep(0,q1+m[1]+2)
    xk[(m[1]+2):(q1+1)]<-seq(min(x),max(x),length=q1-m[1])
    for (i in 1:(m[1]+1)) {xk[i]<-xk[m[1]+2]-(m[1]+2-i)*(xk[m[1]+3]-xk[m[1]+2])}
    for (i in (q1+2):(q1+m[1]+2)) {xk[i]<-xk[q1+1]+(i-q1-1)*(xk[m[1]+3]-xk[m[1]+2])}
    knots[[object$term[1]]] <- xk
  }
  if (n != length(z))
     stop ("arguments of smooth not same dimension")
  if (is.null(zk)) # space knots through the values of the 2nd covariate
  {      zk<-rep(0,q2+m[2]+2)
         zk[(m[2]+2):(q2+1)]<-seq(min(z),max(z),length=q2-m[2])
         for (i in 1:(m[2]+1)) {zk[i]<-zk[m[2]+2]-(m[2]+2-i)*(zk[m[2]+3]-zk[m[2]+2])}
         for (i in (q2+2):(q2+m[2]+2)) {zk[i]<-zk[q2+1]+(i-q2-1)*(zk[m[2]+3]-zk[m[2]+2])}
         knots[[object$term[2]]] <- zk
   }
  if (length(xk)!=nk1 ) # right number of knots?
      stop(paste("there should be ",nk1," supplied knotsfor the x"))
  if (length(zk)!=nk2) # right number of knots?
      stop(paste("there should be ",nk2," supplied knots for z"))
  
 #  get model matrix-------------
  X1 <- splineDesign(xk,x,ord=m[1]+2)
  X2 <- splineDesign(zk,z,ord=m[2]+2)
  X <- matrix(0,n,q1*q2)  # model matrix
  for (i in 1:n)
      {  X[i,] <- X1[i,]%x%X2[i,] # Kronecker product of two rows of marginal model matrices
      }
  # get a matrix Sigma -----------------------
  IS2 <- matrix(0,q2,q2)   # Define marginal matrix of Sigma for concavity
  IS2[1:q2,1] <- rep(1,q2)
  IS2[2:q2,2]<- c(1:(q2-1))
  for (i in 3:q2) IS2[i:q2,i] <- -c(1:(q2-i+1))  

 #  IS1 <- matrix(0,q1,q1)   # Define marginal matrix of Sigma for increasing constraint
 #  IS1[1:q1,1] <- rep(1,q1)
 #    for (j in 2:q1)  IS1[j,2:j] <- 1
  IS1 <- matrix(1,q1,q1)  ## Define marginal matrix of Sigma for increasing constraint
  IS1[upper.tri(IS1)] <-0

  Sig <- IS1 %x% IS2 
  
  # apply identifiability constraint and get model matrix
  X <- X[,2:ncol(X)]%*%Sig[2:ncol(Sig),2:ncol(Sig)]  
  object$X <- X # the finished model matrix with identifiability constraint
 
  # create the penalty matrix
  S <- list()
  I2<- diag(q2)
  I1 <- diff(diag(q1-1),difference=1) 
  Pm1 <- matrix(0,q1-1,q1) # marginal sqrt penalty
  Pm1[2:(q1-1),2:q1] <- I1
  S[[1]]<- Pm1%x%I2

  I2 <- diff(diag(q2-2),difference=1) 
  Pm2 <- matrix(0,q2-1,q2)
  Pm2[3:(q2-1),3:q2] <- I2  # marginal sqrt penalty
  I1 <- diag(q1)
  S[[2]] <- I1%x%Pm2
  
  object$P <- list()
  object$P[[1]] <- S[[1]][2:nrow(S[[1]]),2:ncol(S[[1]])]
  object$P[[2]] <- S[[2]][2:nrow(S[[2]]),2:ncol(S[[2]])]
  object$S <- list()
  object$S[[1]] <- crossprod(object$P[[1]]) ## t(object$P[[1]])%*%object$P[[1]]
  object$S[[2]] <- crossprod(object$P[[2]])  ## t(object$P[[2]])%*%object$P[[2]]

  object$p.ident <- rep(TRUE,q1*q2-1)  ## p.ident is an indicator of which coefficients must be positive (exponentiated)
  object$rank <- ncol(object$X)-1  # penalty rank
  object$null.space.dim <- m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  object$Zc <- diag(q1*q2-1) # identfiability constraint matrix
  object$Zc <- rbind(rep(0,ncol(object$Zc)),object$Zc)
  
  ## store "temicv" specific stuff ...
  object$knots <- list()
  object$knots[[1]] <- xk
  object$knots[[2]] <- zk
  object$m <- m
   
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  class(object)<-"temicv.smooth"  # Give object a class
  object
}

####################################################################


Predict.matrix.temicv.smooth <- function(object, data)
{  ## prediction method function for the `temicx' smooth class
  if (length(object$bs.dim)==1)
      q1 <- q2 <- object$bs.dim # if `k' is supplied as a single number, the same
             ## basis dimension is provided for both marginal smooths
  else  {q1 <- object$bs.dim[1]; q2 <- object$bs.dim[2]}
  
  bm <- marginal.linear.extrapolation(object, data)
  n <- length(data[[object$term[1]]])
  X <- matrix(0,n,q1*q2)  # model matrix
  for ( i in 1:n)
      {  X[i,] <- bm$X1[i,] %x% bm$X2[i,] # Kronecker product of two rows of marginal model matrices
      }
  # get a matrix Sigma -----------------------
  IS2 <- matrix(0,q2,q2)   # Define marginal matrix of Sigma for concavity
  IS2[1:q2,1] <- rep(1,q2)
  IS2[2:q2,2]<- c(1:(q2-1))
  for (i in 3:q2) IS2[i:q2,i] <- -c(1:(q2-i+1))  

 #  IS1 <- matrix(0,q1,q1)   # Define marginal matrix of Sigma for increasing constraint
 #  IS1[1:q1,1] <- rep(1,q1)
 #    for (j in 2:q1)  IS1[j,2:j] <- 1
  IS1 <- matrix(1,q1,q1)  ## Define marginal matrix of Sigma for increasing constraint
  IS1[upper.tri(IS1)] <-0

  Sig <- IS1 %x% IS2 
  # get final model matrix
  X <- X %*%Sig
  X # return the prediction matrix
}





## BELOW TO CORRECT ...



#################################################################
## Tensor product P-spline construction with mixed constraints: 
## decreasing wrt the 1st covariate and concave wrt the second covariate ...
##############################################################


smooth.construct.tedecv.smooth.spec<- function(object, data, knots)
## construction of the bivariate smooth surface with mixed constraints: decreasing
## wrt the 1st covariate and convex wrt the 2nd one...
{ 
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
  nk2 <- q2+m[2]+2 ## number of knots for the 2nd smooth
  if (nk1<=0 || nk2<=0) stop("either k[1] or k[2] too small for m")
  
  ## the values of the first covariate...   
  x <- data[[object$term[1]]]  
  xk <- knots[[object$term[1]]] ## will be NULL if none supplied
  z <- data[[object$term[2]]]  ## the values of the second covariate
  zk <- knots[[object$term[2]]] ## will be NULL if none supplied
  if (is.null(xk)) # space knots through the values of the 1st covariate
  { n<-length(x)
    xk<-rep(0,q1+m[1]+2)
    xk[(m[1]+2):(q1+1)]<-seq(min(x),max(x),length=q1-m[1])
    for (i in 1:(m[1]+1)) {xk[i]<-xk[m[1]+2]-(m[1]+2-i)*(xk[m[1]+3]-xk[m[1]+2])}
    for (i in (q1+2):(q1+m[1]+2)) {xk[i]<-xk[q1+1]+(i-q1-1)*(xk[m[1]+3]-xk[m[1]+2])}
    knots[[object$term[1]]] <- xk
  }
  if (n != length(z))
     stop ("arguments of smooth not same dimension")
  if (is.null(zk)) # space knots through the values of the 2nd covariate
  {      zk<-rep(0,q2+m[2]+2)
         zk[(m[2]+2):(q2+1)]<-seq(min(z),max(z),length=q2-m[2])
         for (i in 1:(m[2]+1)) {zk[i]<-zk[m[2]+2]-(m[2]+2-i)*(zk[m[2]+3]-zk[m[2]+2])}
         for (i in (q2+2):(q2+m[2]+2)) {zk[i]<-zk[q2+1]+(i-q2-1)*(zk[m[2]+3]-zk[m[2]+2])}
         knots[[object$term[2]]] <- zk
   }
  if (length(xk)!=nk1 ) # right number of knots?
      stop(paste("there should be ",nk1," supplied knotsfor the x"))
  if (length(zk)!=nk2) # right number of knots?
      stop(paste("there should be ",nk2," supplied knots for z"))
  
 #  get model matrix-------------
  X1 <- splineDesign(xk,x,ord=m[1]+2)
  X2 <- splineDesign(zk,z,ord=m[2]+2)
  X <- matrix(0,n,q1*q2)  # model matrix
  for (i in 1:n)
      {  X[i,] <- X1[i,]%x%X2[i,] # Kronecker product of two rows of marginal model matrices
      }
  # get a matrix Sigma -----------------------
  IS2 <- matrix(0,q2,q2)   # Define marginal matrix of Sigma for convexity
  IS2[1:q2,1] <- rep(1,q2)
  IS2[2:q2,2]<- -c(1:(q2-1))
  for (i in 3:q2) IS2[i:q2,i] <- c(1:(q2-i+1))  

  IS1 <- matrix(0,q1,q1)   # Define marginal matrix of Sigma for decreasing constraint
  IS1[1:q1,1] <- -rep(1,q1)
    for (j in 2:q1)  IS1[j,2:j] <- -1
  Sig <- IS1 %x% IS2 
  
  # apply identifiability constraint and get model matrix
  X <- X[,2:ncol(X)]%*%Sig[2:ncol(Sig),2:ncol(Sig)]  
  object$X <- X # the finished model matrix with identifiability constraint
 
  # create the penalty matrix
  S <- list()
  I2<- diag(q2)
  I1 <- diff(diag(q1-1),difference=1) 
  Pm1 <- matrix(0,q1-1,q1) # marginal sqrt penalty
  Pm1[2:(q1-1),2:q1] <- I1
  S[[1]]<- Pm1%x%I2

  I2 <- diff(diag(q2-2),difference=1) 
  Pm2 <- matrix(0,q2-1,q2)
  Pm2[3:(q2-1),3:q2] <- I2  # marginal sqrt penalty
  I1 <- diag(q1)
  S[[2]] <- I1%x%Pm2
  
  object$P <- list()
  object$P[[1]] <- S[[1]][2:nrow(S[[1]]),2:ncol(S[[1]])]
  object$P[[2]] <- S[[2]][2:nrow(S[[2]]),2:ncol(S[[2]])]
  object$S <- list()
  object$S[[1]] <- crossprod(object$P[[1]]) ## t(object$P[[1]])%*%object$P[[1]]
  object$S[[2]] <- crossprod(object$P[[2]])  ## t(object$P[[2]])%*%object$P[[2]]

  object$p.ident <- rep(TRUE,q1*q2-1) ## p.ident is an indicator of which coefficients must be positive (exponentiated)
  object$rank <- ncol(object$X)-1  # penalty rank
  object$null.space.dim <- m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  object$Zc <- diag(q1*q2-1) # identfiability constraint matrix
  object$Zc <- rbind(rep(0,ncol(object$Zc)),object$Zc)
  
  ## store "tedecv" specific stuff ...
  object$knots <- list()
  object$knots[[1]] <- xk
  object$knots[[2]] <- zk
  object$m <- m
   
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  class(object)<-"tedecv.smooth"  # Give object a class
  object
}

####################################################################

Predict.matrix.tedecv.smooth <- function(object, data)
{  ## prediction method function for the `tedecv' smooth class
  if (length(object$bs.dim)==1)
      q1 <- q2 <- object$bs.dim # if `k' is supplied as a single number, the same
             ## basis dimension is provided for both marginal smooths
  else  {q1 <- object$bs.dim[1]; q2 <- object$bs.dim[2]}
  
  bm <- marginal.linear.extrapolation(object, data)
  n <- length(data[[object$term[1]]])
  X <- matrix(0,n,q1*q2)  # model matrix
  for ( i in 1:n)
      {  X[i,] <- bm$X1[i,] %x% bm$X2[i,] # Kronecker product of two rows of marginal model matrices
      }
  # get a matrix Sigma -----------------------
  IS2 <- matrix(0,q2,q2)   # Define marginal matrix of Sigma for convexity
  IS2[1:q2,1] <- rep(1,q2)
  IS2[2:q2,2]<- -c(1:(q2-1))
  for (i in 3:q2) IS2[i:q2,i] <- c(1:(q2-i+1))  

  IS1 <- matrix(0,q1,q1)   # Define marginal matrix of Sigma for decreasing constraint
  IS1[1:q1,1] <- -rep(1,q1)
    for (j in 2:q1)  IS1[j,2:j] <- -1
  Sig <- IS1 %x% IS2 
  # get final model matrix
  X <- X %*%Sig
  X # return the prediction matrix
}


################################################################
## Tensor product P-spline construction with mixed constraints:
## decreasing wrt the 1st covariate and convex wrt the 2nd covariate ...
#################################################################


smooth.construct.tedecx.smooth.spec<- function(object, data, knots)
## construction of the bivariate smooth surface with mixed constraints: decreasing
## wrt the first covariate and concave wrt the 2nd one...
{ 
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
  nk2 <- q2+m[2]+2 ## number of knots for the 2nd smooth
  if (nk1<=0 || nk2<=0) stop("either k[1] or k[2] too small for m")
  
  ## the values of the first covariate...   
  x <- data[[object$term[1]]]  
  xk <- knots[[object$term[1]]] ## will be NULL if none supplied
  z <- data[[object$term[2]]]  ## the values of the second covariate
  zk <- knots[[object$term[2]]] ## will be NULL if none supplied
  if (is.null(xk)) # space knots through the values of the 1st covariate
  { n<-length(x)
    xk<-rep(0,q1+m[1]+2)
    xk[(m[1]+2):(q1+1)]<-seq(min(x),max(x),length=q1-m[1])
    for (i in 1:(m[1]+1)) {xk[i]<-xk[m[1]+2]-(m[1]+2-i)*(xk[m[1]+3]-xk[m[1]+2])}
    for (i in (q1+2):(q1+m[1]+2)) {xk[i]<-xk[q1+1]+(i-q1-1)*(xk[m[1]+3]-xk[m[1]+2])}
    knots[[object$term[1]]] <- xk
  }
  if (n != length(z))
     stop ("arguments of smooth not same dimension")
  if (is.null(zk)) # space knots through the values of the 2nd covariate
  {      zk<-rep(0,q2+m[2]+2)
         zk[(m[2]+2):(q2+1)]<-seq(min(z),max(z),length=q2-m[2])
         for (i in 1:(m[2]+1)) {zk[i]<-zk[m[2]+2]-(m[2]+2-i)*(zk[m[2]+3]-zk[m[2]+2])}
         for (i in (q2+2):(q2+m[2]+2)) {zk[i]<-zk[q2+1]+(i-q2-1)*(zk[m[2]+3]-zk[m[2]+2])}
         knots[[object$term[2]]] <- zk
   }
  if (length(xk)!=nk1 ) # right number of knots?
      stop(paste("there should be ",nk1," supplied knotsfor the x"))
  if (length(zk)!=nk2) # right number of knots?
      stop(paste("there should be ",nk2," supplied knots for z"))
  
 #  get model matrix-------------
  X1 <- splineDesign(xk,x,ord=m[1]+2)
  X2 <- splineDesign(zk,z,ord=m[2]+2)
  X <- matrix(0,n,q1*q2)  # model matrix
  for (i in 1:n)
      {  X[i,] <- X1[i,]%x%X2[i,] # Kronecker product of two rows of marginal model matrices
      }
  # get a matrix Sigma -----------------------
  IS2 <- matrix(0,q2,q2)   # Define marginal matrix of Sigma for concavity
  IS2[1:q2,1] <- rep(1,q2)
  IS2[2:q2,2]<- c(1:(q2-1))
  for (i in 3:q2) IS2[i:q2,i] <- -c(1:(q2-i+1))  

  IS1 <- matrix(0,q1,q1)   # Define marginal matrix of Sigma for decreasing constraint
  IS1[1:q1,1] <- -rep(1,q1)
    for (j in 2:q1)  IS1[j,2:j] <- -1
  Sig <- IS1 %x% IS2 
  
  # apply identifiability constraint and get model matrix
  X <- X[,2:ncol(X)]%*%Sig[2:ncol(Sig),2:ncol(Sig)]  
  object$X <- X # the finished model matrix with identifiability constraint
 
  # create the penalty matrix
  S <- list()
  I2<- diag(q2)
  I1 <- diff(diag(q1-1),difference=1) 
  Pm1 <- matrix(0,q1-1,q1) # marginal sqrt penalty
  Pm1[2:(q1-1),2:q1] <- I1
  S[[1]]<- Pm1%x%I2

  I2 <- diff(diag(q2-2),difference=1) 
  Pm2 <- matrix(0,q2-1,q2)
  Pm2[3:(q2-1),3:q2] <- I2  # marginal sqrt penalty
  I1 <- diag(q1)
  S[[2]] <- I1%x%Pm2
  
  object$P <- list()
  object$P[[1]] <- S[[1]][2:nrow(S[[1]]),2:ncol(S[[1]])]
  object$P[[2]] <- S[[2]][2:nrow(S[[2]]),2:ncol(S[[2]])]
  object$S <- list()
  object$S[[1]] <- crossprod(object$P[[1]]) ## t(object$P[[1]])%*%object$P[[1]]
  object$S[[2]] <- crossprod(object$P[[2]])  ## t(object$P[[2]])%*%object$P[[2]]

  object$p.ident <- rep(TRUE,q1*q2-1)  ## p.ident is an indicator of which coefficients must be positive (exponentiated)
  object$rank <- ncol(object$X)-1  # penalty rank
  object$null.space.dim <- m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  object$Zc <- diag(q1*q2-1) # identfiability constraint matrix
  object$Zc <- rbind(rep(0,ncol(object$Zc)),object$Zc)
  
  ## store "tedecx" specific stuff ...
  object$knots <- list()
  object$knots[[1]] <- xk
  object$knots[[2]] <- zk
  object$m <- m
   
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  class(object)<-"tedecx.smooth"  # Give object a class
  object
}

####################################################################

Predict.matrix.tedecx.smooth <- function(object, data)
{  ## prediction method function for the `tedecx' smooth class
  if (length(object$bs.dim)==1)
      q1 <- q2 <- object$bs.dim # if `k' is supplied as a single number, the same
             ## basis dimension is provided for both marginal smooths
  else  {q1 <- object$bs.dim[1]; q2 <- object$bs.dim[2]}
  
  bm <- marginal.linear.extrapolation(object, data)
  n <- length(data[[object$term[1]]])
  X <- matrix(0,n,q1*q2)  # model matrix
  for ( i in 1:n)
      {  X[i,] <- bm$X1[i,] %x% bm$X2[i,] # Kronecker product of two rows of marginal model matrices
      }
  # get a matrix Sigma -----------------------
  IS2 <- matrix(0,q2,q2)   # Define marginal matrix of Sigma for concavity
  IS2[1:q2,1] <- rep(1,q2)
  IS2[2:q2,2]<- c(1:(q2-1))
  for (i in 3:q2) IS2[i:q2,i] <- -c(1:(q2-i+1))  

  IS1 <- matrix(0,q1,q1)   # Define marginal matrix of Sigma for decreasing constraint
  IS1[1:q1,1] <- -rep(1,q1)
    for (j in 2:q1)  IS1[j,2:j] <- -1
  Sig <- IS1 %x% IS2 
  # get final model matrix
  X <- X %*%Sig
  X # return the prediction matrix
}





################################################################
## Tensor product P-spline construction with single concavity constraint
## wrt the 2nd covariate ...
#################################################################

####################################################################

marginal.matrices.tescv.ps <- function(x,z,xk,zk,m,q1,q2)
## function to get marginal model matrices and penalties in the overall 
## model coefficients in case of P-splines basis for the 1st unconstrained smooth
{
  # get marginal model matrix for the first unconstrained smooth... 
  X1 <- splineDesign(xk,x,ord=m[1]+2)
  # get marginal model matrix for the second concave smooth...
  X2 <- splineDesign(zk,z,ord=m[2]+2)

  # create the penalty matrix...
  S <- list()
  # get penalty matrix for the first unconstrained smooth...
  I2 <- diag(q2)
  P <- diff(diag(q1),difference=1)
  S[[1]]<- P %x% I2
  c1 <- c(1,-2,1)
  c2 <- c(1,rep(0,q2-1))
  c <- c1 %x% c2
  for (i in 1:(q1-2)){
       S[[1]][q2*(i-1)+1,(q2*(i-1)+1):(q2*(i-1)+length(c))] <- c
  }
  i <- q1-1
  S[[1]][q2*(i-1)+1,(q2*(i-1)+1):ncol(S[[1]])] <- rep(0,2*q2)

  S[[1]] <- crossprod(S[[1]])  ## t(S[[1]])%*%S[[1]]

  # get penalty for the 2nd concave smooth...
  I2 <- diff(diag(q2-2),difference=1) 
  P <- matrix(0,q2-1,q2)
  P[3:(q2-1),3:q2] <- I2  # marginal sqrt penalty
  I1 <- diag(q1)
  S[[2]] <- I1%x%P
  S[[2]] <- crossprod(S[[2]])  ## t(S[[2]])%*%S[[2]]
   

 list(X1=X1, X2=X2, S=S)
}


smooth.construct.tescv.smooth.spec<- function(object, data, knots)
## construction of the bivariate smooth surface with single concavity constraint
## wrt the 2nd covariate ...
{ 
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
  nk2 <- q2+m[2]+2 ## number of knots for the 2nd smooth
  if (nk1<=0 || nk2<=0) stop("either k[1] or k[2] too small for m")
  
  ## the values of the first covariate...   
  x <- data[[object$term[1]]]  
  xk <- knots[[object$term[1]]] ## will be NULL if none supplied
  z <- data[[object$term[2]]]  ## the values of the second covariate
  zk <- knots[[object$term[2]]] ## will be NULL if none supplied
  if (is.null(xk)) # space knots through the values of the 1st covariate
  { n<-length(x)
    xk<-rep(0,q1+m[1]+2)
    xk[(m[1]+2):(q1+1)]<-seq(min(x),max(x),length=q1-m[1])
    for (i in 1:(m[1]+1)) {xk[i]<-xk[m[1]+2]-(m[1]+2-i)*(xk[m[1]+3]-xk[m[1]+2])}
    for (i in (q1+2):(q1+m[1]+2)) {xk[i]<-xk[q1+1]+(i-q1-1)*(xk[m[1]+3]-xk[m[1]+2])}
    knots[[object$term[1]]] <- xk
  }
  if (n != length(z))
     stop ("arguments of smooth not same dimension")
  if (is.null(zk)) # space knots through the values of the 2nd covariate
  {      zk<-rep(0,q2+m[2]+2)
         zk[(m[2]+2):(q2+1)]<-seq(min(z),max(z),length=q2-m[2])
         for (i in 1:(m[2]+1)) {zk[i]<-zk[m[2]+2]-(m[2]+2-i)*(zk[m[2]+3]-zk[m[2]+2])}
         for (i in (q2+2):(q2+m[2]+2)) {zk[i]<-zk[q2+1]+(i-q2-1)*(zk[m[2]+3]-zk[m[2]+2])}
         knots[[object$term[2]]] <- zk
   }
  if (length(xk)!=nk1 ) # right number of knots?
      stop(paste("there should be ",nk1," supplied knotsfor the x"))
  if (length(zk)!=nk2) # right number of knots?
      stop(paste("there should be ",nk2," supplied knots for z"))
  
 #  get model matrix-------------
  bm <- marginal.matrices.tescv.ps(x,z,xk,zk,m,q1,q2)
  X1 <- bm$X1
  X2 <- bm$X2
  S <- bm$S  
  # get the overall model matrix...
  X <- matrix(0,n,q1*q2)  # model matrix
  for (i in 1:n)
    {  X[i,] <- X1[i,]%x%X2[i,] # Kronecker product of two rows of marginal model matrices
    }
  # get a matrix Sigma -----------------------
  IS2 <- matrix(0,q2,q2)   # Define marginal matrix of Sigma for concavity
  IS2[1:q2,1] <- rep(1,q2)
  IS2[2:q2,2]<- c(1:(q2-1))
  for (i in 3:q2) IS2[i:q2,i] <- -c(1:(q2-i+1))  

  I <- diag(q1) ## identity matrix for the unconstrained marginal
  Sig <- I%x%IS2
  # get model matrix
  X <- X%*%Sig
  # apply identifiability constraint and get model matrix
  D<- diag(q1*q2)
  D<- D[,-1]  ## D[,-((q1-1)*q2+1)]
  ind <- rep(0,q1-1) # get index number for the cells to be changed to "-1"
  for (i in 1:(q1-1)){
        ind[i] <- (i-1)*q2+1
        D[ind[i],ind[i]] <- -1
      }
  for (i in 2:(q1-1))
        D[ind[i],ind[i]-q2] <- 1
  D[((q1-1)*q2+1),((q1-1)*q2+1-q2)] <- 1
  X <- X%*%D 

  object$X <- X # the final model matrix with identifiability constraint
 
  # create the penalty matrix
  object$S <- list()
  object$S[[1]] <- crossprod(D,S[[1]])%*%D  ## t(D)%*%S[[1]]%*%D
  object$S[[2]] <- crossprod(D,S[[2]])%*%D  ## t(D)%*%S[[2]]%*%D

  object$p.ident <- rep(TRUE,q1*q2-1)  
  object$p.ident[ind] <- FALSE ## p.ident is an indicator of which coefficients must be positive (exponentiated)
  object$rank <- ncol(object$X)-1  # penalty rank
  object$null.space.dim <- m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  object$Zc <- D   # identifiability constraint matrix

  
  ## store "tescv" specific stuff ...
  object$knots <- list()
  if (is.null(xk))
     object$knots[[1]] <- rep(0,0,0)
  else object$knots[[1]] <- xk
  object$knots[[2]] <- zk
  object$m <- m
   
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  class(object)<-"tescv.smooth"  # Give object a class
  object
}

####################################################################


Predict.matrix.tescv.smooth <- function(object, data)
{  ## prediction method function for the `tescv' smooth class
  if (length(object$bs.dim)==1)
      q1 <- q2 <- object$bs.dim # if `k' is supplied as a single number, the same
             ## basis dimension is provided for both marginal smooths
  else  {q1 <- object$bs.dim[1]; q2 <- object$bs.dim[2]}
  
  bm <- marginal.linear.extrapolation(object, data)
  n <- length(data[[object$term[1]]])
  X <- matrix(0,n,q1*q2)  # model matrix
  for ( i in 1:n)
      {  X[i,] <- bm$X1[i,] %x% bm$X2[i,] # Kronecker product of two rows of marginal model matrices
      }
  # get a matrix Sigma -----------------------
  IS2 <- matrix(0,q2,q2)   # Define marginal matrix of Sigma for concavity
  IS2[1:q2,1] <- rep(1,q2)
  IS2[2:q2,2]<- c(1:(q2-1))
  for (i in 3:q2) IS2[i:q2,i] <- -c(1:(q2-i+1))  
  
  I <- diag(q1) ## identity matrix for the unconstrained marginal
  Sig <- I%x%IS2
  # get final model matrix
  X <- X %*%Sig
  X # return the prediction matrix
}


################################################################
## Tensor product P-spline construction with single convexity 
## constraint wrt the second covariate ...
################################################################


smooth.construct.tescx.smooth.spec<- function(object, data, knots)
## construction of the bivariate smooth surface with single convexity 
## constraint wrt the second covariate ...
{ 
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
  nk2 <- q2+m[2]+2 ## number of knots for the 2nd smooth
  if (nk1<=0 || nk2<=0) stop("either k[1] or k[2] too small for m")
  
  ## the values of the first covariate...   
  x <- data[[object$term[1]]]  
  xk <- knots[[object$term[1]]] ## will be NULL if none supplied
  z <- data[[object$term[2]]]  ## the values of the second covariate
  zk <- knots[[object$term[2]]] ## will be NULL if none supplied
  if (is.null(xk)) # space knots through the values of the 1st covariate
  { n<-length(x)
    xk<-rep(0,q1+m[1]+2)
    xk[(m[1]+2):(q1+1)]<-seq(min(x),max(x),length=q1-m[1])
    for (i in 1:(m[1]+1)) {xk[i]<-xk[m[1]+2]-(m[1]+2-i)*(xk[m[1]+3]-xk[m[1]+2])}
    for (i in (q1+2):(q1+m[1]+2)) {xk[i]<-xk[q1+1]+(i-q1-1)*(xk[m[1]+3]-xk[m[1]+2])}
    knots[[object$term[1]]] <- xk
  }
  if (n != length(z))
     stop ("arguments of smooth not same dimension")
  if (is.null(zk)) # space knots through the values of the 2nd covariate
  {      zk<-rep(0,q2+m[2]+2)
         zk[(m[2]+2):(q2+1)]<-seq(min(z),max(z),length=q2-m[2])
         for (i in 1:(m[2]+1)) {zk[i]<-zk[m[2]+2]-(m[2]+2-i)*(zk[m[2]+3]-zk[m[2]+2])}
         for (i in (q2+2):(q2+m[2]+2)) {zk[i]<-zk[q2+1]+(i-q2-1)*(zk[m[2]+3]-zk[m[2]+2])}
         knots[[object$term[2]]] <- zk
   }
  if (length(xk)!=nk1 ) # right number of knots?
      stop(paste("there should be ",nk1," supplied knotsfor the x"))
  if (length(zk)!=nk2) # right number of knots?
      stop(paste("there should be ",nk2," supplied knots for z"))
  
 #  get model matrix-------------
  bm <- marginal.matrices.tescv.ps(x,z,xk,zk,m,q1,q2)
  X1 <- bm$X1
  X2 <- bm$X2
  S <- bm$S  
  # get the overall model matrix...
  X <- matrix(0,n,q1*q2)  # model matrix
  for (i in 1:n)
    {  X[i,] <- X1[i,]%x%X2[i,] # Kronecker product of two rows of marginal model matrices
    }
  # get a matrix Sigma -----------------------
  IS2 <- matrix(0,q2,q2)   # Define marginal matrix of Sigma for convexity
  IS2[1:q2,1] <- rep(1,q2)
  IS2[2:q2,2]<- -c(1:(q2-1))
  for (i in 3:q2) IS2[i:q2,i] <- c(1:(q2-i+1))  

  I <- diag(q1) ## identity matrix for the unconstrained marginal
  Sig <- I%x%IS2
  
  # get model matrix
  X <- X%*%Sig
  # apply identifiability constraint 
  D<- diag(q1*q2)
  D<- D[,-1]  ##  D[,-((q1-1)*q2+1)]
  ind <- rep(0,q1-1) # get index number for the cells to be changed to "-1"
  for (i in 1:(q1-1)){
        ind[i] <- (i-1)*q2+1
        D[ind[i],ind[i]] <- -1
      }
  for (i in 2:(q1-1))
        D[ind[i],ind[i]-q2] <- 1
  D[((q1-1)*q2+1),((q1-1)*q2+1-q2)] <- 1
  X <- X%*%D 

  object$X <- X # the finished model matrix with identifiability constraint
 
  # create the penalty matrix
  object$S <- list()
  object$S[[1]] <- crossprod(D,S[[1]])%*%D  ## t(D)%*%S[[1]]%*%D
  object$S[[2]] <- crossprod(D,S[[2]])%*%D  ## t(D)%*%S[[2]]%*%D

  object$p.ident <- rep(TRUE,q1*q2-1)  
  object$p.ident[ind] <- FALSE ## p.ident is an indicator of which coefficients must be positive (exponentiated)
  object$rank <- ncol(object$X)-1  # penalty rank
  object$null.space.dim <- m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  object$Zc <- D   # identifiability constraint matrix
  
  ## store "tescx" specific stuff ...
  object$knots <- list()
  if (is.null(xk))
     object$knots[[1]] <- rep(0,0,0)
  else object$knots[[1]] <- xk
  object$knots[[2]] <- zk
  object$m <- m
      
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  class(object)<-"tescx.smooth"  # Give object a class
  object
}

####################################################################


Predict.matrix.tescx.smooth <- function(object, data)
{  ## prediction method function for the `tescx' smooth class
  if (length(object$bs.dim)==1)
      q1 <- q2 <- object$bs.dim # if `k' is supplied as a single number, the same
             ## basis dimension is provided for both marginal smooths
  else  {q1 <- object$bs.dim[1]; q2 <- object$bs.dim[2]}
  
  bm <- marginal.linear.extrapolation(object, data)
  n <- length(data[[object$term[1]]])
  X <- matrix(0,n,q1*q2)  # model matrix
  for ( i in 1:n)
      {  X[i,] <- bm$X1[i,] %x% bm$X2[i,] # Kronecker product of two rows of marginal model matrices
      }
  # get a matrix Sigma -----------------------
  IS2 <- matrix(0,q2,q2)   # Define marginal matrix of Sigma for convexity
  IS2[1:q2,1] <- rep(1,q2)
  IS2[2:q2,2]<- -c(1:(q2-1))
  for (i in 3:q2) IS2[i:q2,i] <- c(1:(q2-i+1))  

  I <- diag(q1) ## identity matrix for the unconstrained marginal
  Sig <- I%x%IS2
  # get final model matrix
  X <- X %*%Sig
  X # return the prediction matrix
}


#############################################################
## Tensor product P-spline construction with double concavity 
## constraint ...
##############################################################

smooth.construct.tecvcv.smooth.spec<- function(object, data, knots)
## construction of the bivariate smooth surface with double concavity constraint...
{ 
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
  nk2 <- q2+m[2]+2 ## number of knots for the 2nd smooth
  if (nk1<=0 || nk2<=0) stop("either k[1] or k[2] too small for m")
  
  ## the values of the first covariate...   
  x <- data[[object$term[1]]]  
  xk <- knots[[object$term[1]]] ## will be NULL if none supplied
  z <- data[[object$term[2]]]  ## the values of the second covariate
  zk <- knots[[object$term[2]]] ## will be NULL if none supplied
  if (is.null(xk)) # space knots through the values of the 1st covariate
  { n<-length(x)
    xk<-rep(0,q1+m[1]+2)
    xk[(m[1]+2):(q1+1)]<-seq(min(x),max(x),length=q1-m[1])
    for (i in 1:(m[1]+1)) {xk[i]<-xk[m[1]+2]-(m[1]+2-i)*(xk[m[1]+3]-xk[m[1]+2])}
    for (i in (q1+2):(q1+m[1]+2)) {xk[i]<-xk[q1+1]+(i-q1-1)*(xk[m[1]+3]-xk[m[1]+2])}
    knots[[object$term[1]]] <- xk
  }
  if (n != length(z))
     stop ("arguments of smooth not same dimension")
  if (is.null(zk)) # space knots through the values of the 2nd covariate
  {      zk<-rep(0,q2+m[2]+2)
         zk[(m[2]+2):(q2+1)]<-seq(min(z),max(z),length=q2-m[2])
         for (i in 1:(m[2]+1)) {zk[i]<-zk[m[2]+2]-(m[2]+2-i)*(zk[m[2]+3]-zk[m[2]+2])}
         for (i in (q2+2):(q2+m[2]+2)) {zk[i]<-zk[q2+1]+(i-q2-1)*(zk[m[2]+3]-zk[m[2]+2])}
         knots[[object$term[2]]] <- zk
   }
  if (length(xk)!=nk1 ) # right number of knots?
      stop(paste("there should be ",nk1," supplied knotsfor the x"))
  if (length(zk)!=nk2) # right number of knots?
      stop(paste("there should be ",nk2," supplied knots for z"))
  
 #  get model matrix-------------
  X1 <- splineDesign(xk,x,ord=m[1]+2)
  X2 <- splineDesign(zk,z,ord=m[2]+2)
  X <- matrix(0,n,q1*q2)  # model matrix
  for (i in 1:n)
      {  X[i,] <- X1[i,]%x%X2[i,] # Kronecker product of two rows of marginal model matrices
      }
  # get a matrix Sigma -----------------------
  IS2 <- matrix(0,q2,q2)   # Define marginal matrix of Sigma for concavity
  IS2[1:q2,1] <- rep(1,q2)
  IS2[2:q2,2]<- c(1:(q2-1))
  for (i in 3:q2) IS2[i:q2,i] <- -c(1:(q2-i+1))  

  IS1 <- matrix(0,q1,q1)   # Define marginal matrix of Sigma for concavity constraint wrt the first covariate
  IS1[1:q1,1] <- rep(1,q1)
  IS1[2:q1,2]<- c(1:(q1-1))
  for (i in 3:q1) IS1[i:q1,i] <- -c(1:(q1-i+1))
  Sig <- IS1 %x% IS2 
  
  # apply identifiability constraint and get model matrix
  X <- X[,2:ncol(X)]%*%Sig[2:ncol(Sig),2:ncol(Sig)]  
  object$X <- X # the finished model matrix with identifiability constraint
 
  # create the penalty matrix
  S <- list()
  I2<- diag(q2)
  I1 <- diff(diag(q1-1),difference=1) 
  Pm1 <- matrix(0,q1-1,q1) # marginal sqrt penalty
  Pm1[2:(q1-1),2:q1] <- I1
  S[[1]]<- Pm1%x%I2

  I2 <- diff(diag(q2-2),difference=1) 
  Pm2 <- matrix(0,q2-1,q2)
  Pm2[3:(q2-1),3:q2] <- I2  # marginal sqrt penalty
  I1 <- diag(q1)
  S[[2]] <- I1%x%Pm2
  
  object$P <- list()
  object$P[[1]] <- S[[1]][2:nrow(S[[1]]),2:ncol(S[[1]])]
  object$P[[2]] <- S[[2]][2:nrow(S[[2]]),2:ncol(S[[2]])]
  object$S <- list()
  object$S[[1]] <- crossprod(object$P[[1]]) ## t(object$P[[1]])%*%object$P[[1]]
  object$S[[2]] <- crossprod(object$P[[2]])  ## t(object$P[[2]])%*%object$P[[2]]

  object$p.ident <- rep(TRUE,q1*q2-1) ## p.ident is an indicator of which coefficients must be positive (exponentiated) 
  object$rank <- ncol(object$X)-1  # penalty rank
  object$null.space.dim <- m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  object$Zc <- diag(q1*q2-1) # identfiability constraint matrix
  object$Zc <- rbind(rep(0,ncol(object$Zc)),object$Zc)
  
  ## store "tecvcv" specific stuff ...
  object$knots <- list()
  object$knots[[1]] <- xk
  object$knots[[2]] <- zk
  object$m <- m
   
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  class(object)<-"tecvcv.smooth"  # Give object a class
  object
}

####################################################################

Predict.matrix.tecvcv.smooth <- function(object, data)
{  ## prediction method function for the `tecvcv' smooth class
  if (length(object$bs.dim)==1)
      q1 <- q2 <- object$bs.dim # if `k' is supplied as a single number, the same
             ## basis dimension is provided for both marginal smooths
  else  {q1 <- object$bs.dim[1]; q2 <- object$bs.dim[2]}
  
  bm <- marginal.linear.extrapolation(object, data)
  n <- length(data[[object$term[1]]])
  X <- matrix(0,n,q1*q2)  # model matrix
  for ( i in 1:n)
      {  X[i,] <- bm$X1[i,] %x% bm$X2[i,] # Kronecker product of two rows of marginal model matrices
      }
  # get a matrix Sigma -----------------------
  IS2 <- matrix(0,q2,q2)   # Define marginal matrix of Sigma for concavity
  IS2[1:q2,1] <- rep(1,q2)
  IS2[2:q2,2]<- c(1:(q2-1))
  for (i in 3:q2) IS2[i:q2,i] <- -c(1:(q2-i+1))  

  IS1 <- matrix(0,q1,q1)   # Define marginal matrix of Sigma for concavity constraint wrt the 1st covariate......
  IS1[1:q1,1] <- rep(1,q1)
  IS1[2:q1,2]<- c(1:(q1-1))
  for (i in 3:q1) IS1[i:q1,i] <- -c(1:(q1-i+1))
  Sig <- IS1 %x% IS2 
  # get final model matrix
  X <- X %*%Sig
  X # return the prediction matrix
}

################################################################
## Tensor product P-spline construction with double convexity 
## constraint...
################################################################


smooth.construct.tecxcx.smooth.spec<- function(object, data, knots)
## construction of the bivariate smooth surface with double
## convexity constraint...
{ 
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
  nk2 <- q2+m[2]+2 ## number of knots for the 2nd smooth
  if (nk1<=0 || nk2<=0) stop("either k[1] or k[2] too small for m")
  
  ## the values of the first covariate...   
  x <- data[[object$term[1]]]  
  xk <- knots[[object$term[1]]] ## will be NULL if none supplied
  z <- data[[object$term[2]]]  ## the values of the second covariate
  zk <- knots[[object$term[2]]] ## will be NULL if none supplied
  if (is.null(xk)) # space knots through the values of the 1st covariate
  { n<-length(x)
    xk<-rep(0,q1+m[1]+2)
    xk[(m[1]+2):(q1+1)]<-seq(min(x),max(x),length=q1-m[1])
    for (i in 1:(m[1]+1)) {xk[i]<-xk[m[1]+2]-(m[1]+2-i)*(xk[m[1]+3]-xk[m[1]+2])}
    for (i in (q1+2):(q1+m[1]+2)) {xk[i]<-xk[q1+1]+(i-q1-1)*(xk[m[1]+3]-xk[m[1]+2])}
    knots[[object$term[1]]] <- xk
  }
  if (n != length(z))
     stop ("arguments of smooth not same dimension")
  if (is.null(zk)) # space knots through the values of the 2nd covariate
  {      zk<-rep(0,q2+m[2]+2)
         zk[(m[2]+2):(q2+1)]<-seq(min(z),max(z),length=q2-m[2])
         for (i in 1:(m[2]+1)) {zk[i]<-zk[m[2]+2]-(m[2]+2-i)*(zk[m[2]+3]-zk[m[2]+2])}
         for (i in (q2+2):(q2+m[2]+2)) {zk[i]<-zk[q2+1]+(i-q2-1)*(zk[m[2]+3]-zk[m[2]+2])}
         knots[[object$term[2]]] <- zk
   }
  if (length(xk)!=nk1 ) # right number of knots?
      stop(paste("there should be ",nk1," supplied knotsfor the x"))
  if (length(zk)!=nk2) # right number of knots?
      stop(paste("there should be ",nk2," supplied knots for z"))
  
 #  get model matrix-------------
  X1 <- splineDesign(xk,x,ord=m[1]+2)
  X2 <- splineDesign(zk,z,ord=m[2]+2)
  X <- matrix(0,n,q1*q2)  # model matrix
  for (i in 1:n)
      {  X[i,] <- X1[i,]%x%X2[i,] # Kronecker product of two rows of marginal model matrices
      }
  # get a matrix Sigma -----------------------
  IS2 <- matrix(0,q2,q2)   # Define marginal matrix of Sigma for convexity
  IS2[1:q2,1] <- rep(1,q2)
  IS2[2:q2,2]<- -c(1:(q2-1))
  for (i in 3:q2) IS2[i:q2,i] <- c(1:(q2-i+1))  

  IS1 <- matrix(0,q1,q1)   # Define marginal matrix of Sigma for convexity
  IS1[1:q1,1] <- rep(1,q1)
  IS1[2:q1,2]<- -c(1:(q1-1))
  for (i in 3:q1) IS1[i:q1,i] <- c(1:(q1-i+1))
  Sig <- IS1 %x% IS2 
  
  # apply identifiability constraint and get model matrix
  X <- X[,2:ncol(X)]%*%Sig[2:ncol(Sig),2:ncol(Sig)]  
  object$X <- X # the finished model matrix with identifiability constraint
 
  # create the penalty matrix
  S <- list()
  I2<- diag(q2)
  I1 <- diff(diag(q1-1),difference=1) 
  Pm1 <- matrix(0,q1-1,q1) # marginal sqrt penalty
  Pm1[2:(q1-1),2:q1] <- I1
  S[[1]]<- Pm1%x%I2

  I2 <- diff(diag(q2-2),difference=1) 
  Pm2 <- matrix(0,q2-1,q2)
  Pm2[3:(q2-1),3:q2] <- I2  # marginal sqrt penalty
  I1 <- diag(q1)
  S[[2]] <- I1%x%Pm2
  
  object$P <- list()
  object$P[[1]] <- S[[1]][2:nrow(S[[1]]),2:ncol(S[[1]])]
  object$P[[2]] <- S[[2]][2:nrow(S[[2]]),2:ncol(S[[2]])]
  object$S <- list()
  object$S[[1]] <- crossprod(object$P[[1]]) ## t(object$P[[1]])%*%object$P[[1]]
  object$S[[2]] <- crossprod(object$P[[2]])  ## t(object$P[[2]])%*%object$P[[2]]
  object$p.ident <- rep(TRUE,q1*q2-1)    ## p.ident is an indicator of which coefficients must be positive (exponentiated) 
  object$rank <- ncol(object$X)-1  # penalty rank
  object$null.space.dim <- m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  object$Zc <- diag(q1*q2-1) # identfiability constraint matrix
  object$Zc <- rbind(rep(0,ncol(object$Zc)),object$Zc)
  
  ## store "tecxcx" specific stuff ...
  object$knots <- list()
  object$knots[[1]] <- xk
  object$knots[[2]] <- zk
  object$m <- m
   
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  class(object)<-"tecxcx.smooth"  # Give object a class
  object
}

####################################################################


Predict.matrix.tecxcx.smooth <- function(object, data)
{  ## prediction method function for the `tecxcx' smooth class
  if (length(object$bs.dim)==1)
      q1 <- q2 <- object$bs.dim # if `k' is supplied as a single number, the same
             ## basis dimension is provided for both marginal smooths
  else  {q1 <- object$bs.dim[1]; q2 <- object$bs.dim[2]}
  
  bm <- marginal.linear.extrapolation(object, data)
  n <- length(data[[object$term[1]]])
  X <- matrix(0,n,q1*q2)  # model matrix
  for ( i in 1:n)
      {  X[i,] <- bm$X1[i,] %x% bm$X2[i,] # Kronecker product of two rows of marginal model matrices
      }
  # get a matrix Sigma -----------------------
  IS2 <- matrix(0,q2,q2)   # Define marginal matrix of Sigma for convexity
  IS2[1:q2,1] <- rep(1,q2)
  IS2[2:q2,2]<- -c(1:(q2-1))
  for (i in 3:q2) IS2[i:q2,i] <- c(1:(q2-i+1))  

  IS1 <- matrix(0,q1,q1)   # Define marginal matrix of Sigma for convexity
  IS1[1:q1,1] <- rep(1,q1)
  IS1[2:q1,2]<- -c(1:(q1-1))
  for (i in 3:q1) IS1[i:q1,i] <- c(1:(q1-i+1))
  Sig <- IS1 %x% IS2 
  # get final model matrix
  X <- X %*%Sig
  X # return the prediction matrix
}


################################################################
## Tensor product P-spline construction with convexity 
## constraint along the 1st covarite and concavity along the 2nd...
################################################################


smooth.construct.tecxcv.smooth.spec<- function(object, data, knots)
## construction of the bivariate smooth surface with convexity 
## constraint along the 1st covarite and concavity along the 2nd...
{ 
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
  nk2 <- q2+m[2]+2 ## number of knots for the 2nd smooth
  if (nk1<=0 || nk2<=0) stop("either k[1] or k[2] too small for m")
  
  ## the values of the first covariate...   
  x <- data[[object$term[1]]]  
  xk <- knots[[object$term[1]]] ## will be NULL if none supplied
  z <- data[[object$term[2]]]  ## the values of the second covariate
  zk <- knots[[object$term[2]]] ## will be NULL if none supplied
  if (is.null(xk)) # space knots through the values of the 1st covariate
  { n<-length(x)
    xk<-rep(0,q1+m[1]+2)
    xk[(m[1]+2):(q1+1)]<-seq(min(x),max(x),length=q1-m[1])
    for (i in 1:(m[1]+1)) {xk[i]<-xk[m[1]+2]-(m[1]+2-i)*(xk[m[1]+3]-xk[m[1]+2])}
    for (i in (q1+2):(q1+m[1]+2)) {xk[i]<-xk[q1+1]+(i-q1-1)*(xk[m[1]+3]-xk[m[1]+2])}
    knots[[object$term[1]]] <- xk
  }
  if (n != length(z))
     stop ("arguments of smooth not same dimension")
  if (is.null(zk)) # space knots through the values of the 2nd covariate
  {      zk<-rep(0,q2+m[2]+2)
         zk[(m[2]+2):(q2+1)]<-seq(min(z),max(z),length=q2-m[2])
         for (i in 1:(m[2]+1)) {zk[i]<-zk[m[2]+2]-(m[2]+2-i)*(zk[m[2]+3]-zk[m[2]+2])}
         for (i in (q2+2):(q2+m[2]+2)) {zk[i]<-zk[q2+1]+(i-q2-1)*(zk[m[2]+3]-zk[m[2]+2])}
         knots[[object$term[2]]] <- zk
   }
  if (length(xk)!=nk1 ) # right number of knots?
      stop(paste("there should be ",nk1," supplied knotsfor the x"))
  if (length(zk)!=nk2) # right number of knots?
      stop(paste("there should be ",nk2," supplied knots for z"))
  
 #  get model matrix-------------
  X1 <- splineDesign(xk,x,ord=m[1]+2)
  X2 <- splineDesign(zk,z,ord=m[2]+2)
  X <- matrix(0,n,q1*q2)  # model matrix
  for (i in 1:n)
      {  X[i,] <- X1[i,]%x%X2[i,] # Kronecker product of two rows of marginal model matrices
      }
  # get a matrix Sigma -----------------------
  IS1 <- matrix(0,q1,q1)   # Define marginal matrix of Sigma for convexity
  IS1[1:q1,1] <- rep(1,q1)
  IS1[2:q1,2]<- -c(1:(q1-1))
  for (i in 3:q1) IS1[i:q1,i] <- c(1:(q1-i+1))  

  IS2 <- matrix(0,q2,q2)   # Define marginal matrix of Sigma for concavity
  IS2[1:q2,1] <- rep(1,q2)
  IS2[2:q2,2]<- c(1:(q2-1))
  for (i in 3:q2) IS2[i:q2,i] <- -c(1:(q2-i+1))
  Sig <- IS1 %x% IS2 
  
  # apply identifiability constraint and get model matrix
  X <- X[,2:ncol(X)]%*%Sig[2:ncol(Sig),2:ncol(Sig)]  
  object$X <- X # the finished model matrix with identifiability constraint
 
  # create the penalty matrix
  S <- list()
  I2<- diag(q2)
  I1 <- diff(diag(q1-1),difference=1) 
  Pm1 <- matrix(0,q1-1,q1) # marginal sqrt penalty
  Pm1[2:(q1-1),2:q1] <- I1
  S[[1]]<- Pm1%x%I2

  I2 <- diff(diag(q2-2),difference=1) 
  Pm2 <- matrix(0,q2-1,q2)
  Pm2[3:(q2-1),3:q2] <- I2  # marginal sqrt penalty
  I1 <- diag(q1)
  S[[2]] <- I1%x%Pm2
  
  object$P <- list()
  object$P[[1]] <- S[[1]][2:nrow(S[[1]]),2:ncol(S[[1]])]
  object$P[[2]] <- S[[2]][2:nrow(S[[2]]),2:ncol(S[[2]])]
  object$S <- list()
  object$S[[1]] <- crossprod(object$P[[1]]) ## t(object$P[[1]])%*%object$P[[1]]
  object$S[[2]] <- crossprod(object$P[[2]])  ## t(object$P[[2]])%*%object$P[[2]]
  object$p.ident <- rep(TRUE,q1*q2-1)   ## p.ident is an indicator of which coefficients must be positive (exponentiated) 
  object$rank <- ncol(object$X)-1  # penalty rank
  object$null.space.dim <- m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  object$Zc <- diag(q1*q2-1) # identfiability constraint matrix
  object$Zc <- rbind(rep(0,ncol(object$Zc)),object$Zc)
  
  ## store "tecxcv" specific stuff ...
  object$knots <- list()
  object$knots[[1]] <- xk
  object$knots[[2]] <- zk
  object$m <- m
   
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  class(object)<-"tecxcv.smooth"  # Give object a class
  object
}

####################################################################


Predict.matrix.tecxcv.smooth <- function(object, data)
{  ## prediction method function for the `tecxcv' smooth class
  if (length(object$bs.dim)==1)
      q1 <- q2 <- object$bs.dim # if `k' is supplied as a single number, the same
             ## basis dimension is provided for both marginal smooths
  else  {q1 <- object$bs.dim[1]; q2 <- object$bs.dim[2]}
  
  bm <- marginal.linear.extrapolation(object, data)
  n <- length(data[[object$term[1]]])
  X <- matrix(0,n,q1*q2)  # model matrix
  for ( i in 1:n)
      {  X[i,] <- bm$X1[i,] %x% bm$X2[i,] # Kronecker product of two rows of marginal model matrices
      }
  # get a matrix Sigma -----------------------
  IS1 <- matrix(0,q1,q1)   # Define marginal matrix of Sigma for convexity
  IS1[1:q1,1] <- rep(1,q1)
  IS1[2:q1,2]<- -c(1:(q1-1))
  for (i in 3:q1) IS1[i:q1,i] <- c(1:(q1-i+1))  

  IS2 <- matrix(0,q2,q2)   # Define marginal matrix of Sigma for concavity
  IS2[1:q2,1] <- rep(1,q2)
  IS2[2:q2,2]<- c(1:(q2-1))
  for (i in 3:q2) IS2[i:q2,i] <- -c(1:(q2-i+1))
  Sig <- IS1 %x% IS2 
  # get final model matrix
  X <- X %*%Sig
  X # return the prediction matrix
}

