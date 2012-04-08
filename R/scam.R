

## the wrapper overall Function to fit scam...

scam <- function(formula,family=gaussian(),data=list(),gamma=1,sp=NULL,
            weights=NULL,offset=NULL,optimizer="bfgs", 
            optim.method=c("Nelder-Mead","fd"),
           scale=0,epsilon=1e-8, check.analytical=FALSE, del=1e-4,
              start=NULL, etastart, mustart)
{  
  ## scale - scale parameter of the exponential deistribution as in gam(mgcv)
  ## min.edf - minimum allowed value of the edf for each parameter 
  ## epsilon - convergence control in the full Newton method
  ## optimizer - numerical optimization method to use to optimize the smoothing 
  ##             parameter estimation criterion: "bfgs", "optim", "nlm", "nlm.fd"
  ## optim.method - if optimizer=="optim" then the first argument of optim.method specifies the method,
  ##             and the second can be either "fd" for finite-difference approximation of the gradient or
  ##             "grad" - to use analytical gradient of gcv/ubre
  ## check.analytical - logical whether the analytical gradient of GCV/UBRE should be checked
  ## del - increment for finite differences when checking analytical gradients
 ### ------------------------------------------------------------------

  ## Setting from mgcv(gam).......
  require(mgcv)
  G<-gam(formula,family,data,fit=FALSE)
  n.terms <- length(G$smooth)  # number of smooths in the model
  n <- nrow(G$X)
  intercept <- G$intercept ## TRUE or FALSE
  if (is.null(offset)) 
        offset <- rep.int(0, n)
  if (is.null(weights)) 
        weights <- rep.int(1,n)
  fam.name <- G$family[1]
  if (scale == 0) {
        if (fam.name == "binomial" || fam.name == "poisson") 
            sig2 <- 1
        else sig2 <- -1
  }
  else {
        sig2 <- scale
  }
 if (sig2 > 0) scale.known <- TRUE else scale.known <- FALSE


  ## get Penalty matrices and 
  ## vector of identifications for exponentiated model coefficients...
  Q <- penalty_pident(G)

  ## checking sp...
  if (!is.null(sp)) {
       neg <- FALSE
       if (length(sp)!= length(G$off)){
            warning("Supplied smoothing parameter vector is too short - ignored.")
            sp <- NULL       
       }
       else if (sum(is.na(sp))) {
              warning("NA's in supplied smoothing parameter vector - ignoring.")
              sp <- NULL
           } 
       else {
             for ( i in 1:length(sp)) { # cheking negative values... 
                 if (sp[i] < 0){
                     warning("Supplied smoothing parameter vector has negative values - ignored.")
                     neg <- TRUE  
                  } 
            }
       }
       if (neg) sp <- NULL
  }
  
  ## Create new environments with `start' initially empty
  ee <- new.env()
  assign("start",rep(0,0),envir=ee)
  eb <- new.env()
  assign("dbeta.start",rep(0,0),envir=eb)
  esp <- new.env()
  assign("sp.last",rep(0,0),envir=esp)
  
  q.f <- rep(0,n.terms)
  for (i in 1:n.terms) {
             q.f[i] <- ncol(G$smooth[[i]]$S[[1]]) +1
  }
  G$S <- Q$S
  G$q.f <- q.f
  G$q0 <- G$off[1]-1  ## number of the parameters of the strictly parametric model
  G$p.ident <- Q$p.ident  # vector of 0's & 1's for the model parameters identification: 
  G$n.terms <- n.terms   ## number of the smooth terms in the mono-GAM
  G$intercept <- intercept
  G$weights <- weights
  G$sig2 <- sig2
  G$scale.known <- scale.known
  object <- list() 
  if (is.null(sp)) {
      ## get initial estimates of the smoothing parameter...
      ## ADDED code......  initialize y...
      start <- etastart <- mustart <- NULL
      y <- G$y; family <- G$family
      nobs <- NROW(y)
      eval(family$initialize)
      G$y <- y  ## needed to set factor response values as numeric
      ## end of ADDED CODE...
      def.sp <- initial.sp.scam (G,Q,q.f=q.f,n.terms=n.terms,family=family,SVD=TRUE,
              ee,eb,esp,intercept=intercept,offset=G$offset,
              weights=weights,epsilon=epsilon)
      rho <- log(def.sp) ## get initial log(sp) ...

     ## minimize GCV/UBRE by optimizer....
      ptm <- proc.time()
      re <- estimate.scam(G=G,optimizer=optimizer,optim.method=optim.method,
               rho=rho, gamma=gamma,ee=ee,eb=eb,esp=esp, 
              check.analytical=check.analytical, del=del,epsilon=epsilon)
      CPU.time <- proc.time()-ptm
      best <- re
      object$gcv.ubre <- re$gcv.ubre
      object$dgcv.ubre <- re$dgcv.ubre
      object$aic <- re$aic
      best$p.ident <- Q$p.ident
      best$S <- Q$S
      object$termcode <- re$termcode
      if (optimizer == "bfgs"){
            object$check.grad <- re$check.grad
            object$dgcv.ubre.check <- re$dgcv.ubre.check
      }
  }
  else { ## no GCV minimization if sp is given...
      best <- scam.fit(G=G, sp=sp,SVD=TRUE,ee,eb,esp,epsilon=epsilon)
      object$aic <- best$aic
  }
 
  ## post-fitting values...
  
  best$n.smooth <- object$n.smooth <- n.terms
  best$formula <- object$formula <- formula
  best$family <- object$family <- G$family
  best$smooth <- object$smooth <- G$smooth
  best$model <- object$model <- G$mf
  object$trA <- best$trA
  object$edf <- best$edf
  if (is.null(sp))
       object$sp <- best$sp
  object$deviance <- best$deviance
  object$residuals <- best$residuals
  object$X <- best$X
  object$conv <- best$conv # whether or not the full Newton method converged
  
  post <- scam.fit.post(y=G$y,X=G$X,object=best,sig2=sig2,offset = G$offset,
                    intercept=intercept,weights=weights,scale.known=scale.known)
  

  object$null.deviance <- post$nulldev 
  object$Ve <- post$Ve
  object$Vp <- post$Vb
  object$Ve.t <- post$Ve.t
  object$Vp.t <- post$Vb.t
  object$sig2 <- post$sig2
  object$coefficients <- best$beta
  object$coefficients.t <- best$beta.t
  object$beta <- best$beta
  object$beta.t <- best$beta.t
  object$pterms <- G$pterms
  object$terms <- G$terms
  object$assign <- G$assign
  object$nsdf <- G$nsdf
  object$y <- G$y
  object$data <- G$mf
  object$offset <- G$offset
  object$scale.known <- scale.known # to be passed in the summary function
  object$prior.weights <-weights # prior weights on observations
  object$weights <- best$w  # final weights used in full Newton iteration
  object$fitted.values <- best$mu
  object$linear.predictors <- best$eta
  object$cmX <- G$cmX
  object$p.ident <- Q$p.ident
  object$intercept <- intercept
  object$gamma <- gamma
  object$iter <- best$iter  # number of iterations of the Full Newton
  object$optimizer <- optimizer
  if (is.null(sp)) 
        object$CPU.time <- CPU.time
  else 
        object$CPU.time <- NULL

  ## get the optimizer info (smoothing parameter selection).....
  if (is.null(sp)){ 
      if (optimizer == "bfgs"){  
     ## get the bfgs info in case of sp selection... 
           object$bfgs.info <- list()
           object$bfgs.info$conv <- re$conv.bfgs  
           object$bfgs.info$iter <- re$iterations 
           object$bfgs.info$grad <- re$dgcv.ubre
      }
      else if (optimizer == "nlm.fd" || optimizer == "nlm"){
           object$nlm.info <- list()
           object$nlm.info$conv <- re$conv 
           object$nlm.info$iter <- re$iterations 
           object$nlm.info$grad <- re$dgcv.ubre
      }
      else if (optimizer=="optim"){
           object$optim.info <- list()
           object$optim.info$conv <- re$conv 
           object$optim.info$iter <- re$iterations 
           object$optim.method <- re$optim.method  
      }
  }
   

  if (scale.known)
      object$method <- "UBRE"
  else
      object$method <- "GCV" 

  if (G$nsdf > 0) 
        term.names <- colnames(G$X)[1:G$nsdf]
  else term.names <- array("", 0)
  if (n.terms) 
        for (i in 1:n.terms) {
            k <- 1
            for (j in G$smooth[[i]]$first.para:G$smooth[[i]]$last.para) {
                term.names[j] <- paste(G$smooth[[i]]$label, ".", 
                  as.character(k), sep = "")
                k <- k + 1
            }
        }
  names(object$coefficients) <- term.names
  names(object$coefficients.t) <- term.names
  ynames <- if (is.matrix(G$y)) 
                rownames(G$y)
            else names(G$y)
  names(object$residuals) <- ynames
  class(object) <- "scam"   
  object
}


##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## function to get initial estimates of smoothing parameters...

initial.sp.scam <- function(G,Q,q.f,n.terms,family,SVD,
              ee,eb,esp,intercept,offset,weights,epsilon)  {
   ## function to get initial estimates of smoothing parameters
   ## step 1 ## set sp=rep(0.5,p) and estimate hessian...
    b <- scam.fit(G=G,sp=rep(0.5,length(G$off)),SVD=TRUE,
              ee,eb,esp,epsilon=epsilon)

   svd.Xw <- try(svd(b$wX1),silent=TRUE)
   if (inherits(svd.Xw, "try-error"))
          svd.Xw <- try(svd(b$wX1,LINPACK = TRUE),silent=TRUE)
   H <- svd.Xw$v%*%(svd.Xw$d^2*t(svd.Xw$v)) - b$E   ## unpenalized Hessian
   ## step 2:...
   n.p <- length(Q$S) ## number of penalty matrices
  def.sp <- array(0,n.p) ## initialize the initial sp values
  j <- 1
  for (i in 1:n.terms){
     for (kk in 1:length(G$smooth[[i]]$S)){
       start <- G$off[j]
       finish <- start + ncol(G$smooth[[i]]$S[[kk]])-1
         # matrix norm of the Hessian elements penalized by S[[kk]]...
       Hi.norm <- sum(H[start:finish,start:finish]*H[start:finish,start:finish]) 
       Si.norm <- sum(G$smooth[[i]]$S[[kk]]*G$smooth[[i]]$S[[kk]])
       def.sp[j] <- (Hi.norm/Si.norm)^0.5
       j <- j+1
     }
  }
 ## Create again new environments with `start' initially empty...
  ee <- new.env()
  assign("start",rep(0,0),envir=ee)
  eb <- new.env()
  assign("dbeta.start",rep(0,0),envir=eb)
  esp <- new.env()
  assign("sp.last",rep(0,0),envir=esp)

##as.numeric(def.sp)
  def.sp
}


##------------------------------------------------------------------------
## function to get list of penalty matrices and 
## vector of parameter identifications .....

penalty_pident <- function(object)
## function to get the list of penalties and vector of model patameters 
## identifications from the gam() setting...
{ n.terms <- length(object$smooth)  # number of terms in the model
  q <- ncol(object$X)          # total number of parameters
  cons.terms <- rep(0,n.terms) # define whether each term is constrained or not
  for (i in 1:n.terms){
        if (!is.null(object$smooth[[i]]$p.ident))
             cons.terms[i] <- 1  
  }
  p.ident <- rep(0,q) # initialize vector of parameter identifications
                      # with `1' - for a parameter to be exponentiated, `0' - otehrwise
  off.terms <- rep(0,n.terms) # starting points for each term
  off <- object$off
  if (n.terms ==length(off)){
        off.terms <- off}
  else 
      { off.terms[1] <- off[1]
       k <- 1
       l <- 1
       while (l<length(off)){
              if (off[l]!=off[l+1]){
                     off.terms[k+1] <- off[l+1] 
                     k <- k+1; l <- l+1 
              } 
              else l <- l+1
       }
     #  if (off[l]!=off[l-1]){
     #                off.terms[k] <- off[l] }
  }
  for (i in 1:n.terms){
        if (cons.terms[i]==1) 
            #  p.ident[off.terms[i]:(off.terms[i+1]-1)] <- 1
              p.ident[off.terms[i]:(off.terms[i]+ncol(object$smooth[[i]]$S[[1]])-1)] <- 
                 object$smooth[[i]]$p.ident
  }
#  if (cons.terms[n.terms]==1) 
#              p.ident[off.terms[n.terms]:q] <- 1
    
  ## getting the list of penalty matrices in terms of the full model
     ## coefficients vector ...
  S <- list()
  j <- 1
  for(i in 1:n.terms)
        { for (kk in 1:length(object$smooth[[i]]$S))
                { S[[j]] <- matrix(0,q,q) # initialize penalty matrix
                 S[[j]][off.terms[i]:(off.terms[i]+ncol(object$smooth[[i]]$S[[kk]])-1),
                 off.terms[i]:(off.terms[i]+ncol(object$smooth[[i]]$S[[kk]])-1)] <-object$smooth[[i]]$S[[kk]]
                 j <- j+1       
          }
  }

 object$S <- S 
 object$p.ident <- p.ident
 object
}




#----------------------------------------------------------------------------
# ***************************************************************************
# Function to fit mono-GAM based on Full Newton method ******************

scam.fit <- function(G,sp, SVD=TRUE,ee,eb,esp, maxit=200,epsilon=1e-8,
                start=NULL, etastart=NULL, mustart=NULL)

# y,X,S,sp,q0,q.f,p.ident,n.terms,family,SVD=TRUE,ee,eb,esp,
 #        maxit=200,epsilon=1e-8,offset,
 #        intercept,weights)   
# function to fit mono-GAM with several terms
# y - values of the response variable
# X - model matrix
# S - list of penalties matrices 
# sp- vector of smoothing parameters
# q0 - number of parameters of the strictly parametric model 
# q.f -vector of the basis dimensions for each smooth term
# n.terms -number of smooth terms
# p.ident - vector of 0's & 1's for the model parameters identification: 
# 1's stand for the parameters which will be exponentiated, 0's - otherwise
# SVD - if TRUE then svd is applied to the augmented `model' matrix, 
# if SVD=FALSE then qr decomposition will be used
# weights - values of constant weights omega_i
{ 
  y <- G$y;  X <- G$X;  S <- G$S
  q0 <- G$q0; q.f <- G$q.f
  p.ident <- G$p.ident; n.terms <- G$n.terms
  family <- G$family; intercept <- G$intercept; offset <- G$offset;
  weights <- G$weights;  
  n <- nobs <- NROW(y)
  q <- ncol(X)
  dg <- fix.family.link(family)
  dv <- fix.family.var(family)

## ADDED CODE........
  nvars <- ncol(X)
  EMPTY <- nvars == 0
  variance <- family$variance
  linkinv <- family$linkinv
  if (!is.function(variance) || !is.function(linkinv)) 
          stop("'family' argument seems not to be a valid family object")
  dev.resids <- family$dev.resids
  aic <- family$aic
  mu.eta <- family$mu.eta
  mu.eta <- family$mu.eta
  if (!is.function(variance) || !is.function(linkinv)) 
        stop("illegal `family' argument")
  valideta <- family$valideta
  if (is.null(valideta)) 
        valideta <- function(eta) TRUE
  validmu <- family$validmu
  if (is.null(validmu)) 
        validmu <- function(mu) TRUE
  if (is.null(mustart)) {
        eval(family$initialize)
  }
  else {
        mukeep <- mustart
        eval(family$initialize)
        mustart <- mukeep
  }
  if (EMPTY) {
      eta <- rep.int(0, nobs) + offset
      if (!valideta(eta)) 
            stop("Invalid linear predictor values in empty model")
      mu <- linkinv(eta)
      if (!validmu(mu)) 
          stop("Invalid fitted means in empty model")
      dev <- sum(dev.resids(y, mu, weights))
      w <- (weights * mu.eta(eta)^2)/variance(mu) ## incorrect for Newton
      residuals <- (y - mu)/mu.eta(eta)
      good <- rep(TRUE, length(residuals))
      boundary <- conv <- TRUE
      coef <- numeric(0)
      iter <- 0
      V <- variance(mu)
      alpha <- dev
      trA <- 0
      GCV <- nobs * alpha/(nobs - gamma * trA)^2
      UBRE <- alpha/nobs - scale + 2 * gamma/n * trA
      scale.est <- alpha/(nobs - trA)
      aic.model <- aic(y, n, mu, weights, dev) +  2 * trA
  } ### end if (EMPTY)
  else {
      eta <- if (!is.null(etastart)) 
          etastart
    ### NB: initialization of parameters differs from gam() and glm() (see below)
    #  else if (!is.null(start)) 
    #      if (length(start) != nvars) 
    #            stop("Length of start should equal ", nvars, 
    #              " and correspond to initial coefs for ", deparse(xnames))
    #      else {
    #            coefold <- start
    #            offset + as.vector(if (NCOL(x) == 1) 
    #               x * start
    #            else x %*% start)
    #      }
      else family$linkfun(mustart)
      mu <- linkinv(eta)
      if (!(validmu(mu) && valideta(eta))) 
            stop("Can't find valid starting values: please specify some")
      
 ## end of added code except the loop for "if  (!EMPTY)" .....
      S.t <- matrix(0,q,q) # define the total sum of the penalties times sp
      n.pen <- length(S) # number of penalties 
      if (length(sp)!=n.pen) stop (paste("length of sp has to be equal to", n.pen))
      for (j in 1:n.pen)
      S.t<-S.t+sp[j]*S[[j]]
      # get sqrt of total penalty matrix...
      rS <- t(mroot(S.t,method="svd"))
      count<-0 # count the number of the parameters to be exponentiated
      for (i in 1:q)
            {if (p.ident[i]==1) count<-count+1}
      ia<-array(0, dim=c(count,2)) # define an index array for the monotone parameters
      iv<-array(0, dim=c(count,1)) # define an index vector for the monotone parameters
      k<-1
      for (i in 1:q)
           {if (p.ident[i]==1) {ia[k,]<-i; iv[k]<-i; k<-k+1}}
      #----------------------------------------
      ## Initialization of parameters start here 
      beta0 <- get("start",envir=ee)
      dbeta0<-get("dbeta.start",envir=eb)
      sp.old<-get("sp.last",envir=esp)
      if (length(beta0)==0) {
          # list argument to pcls for initializing model coefficients
          M<-list(X=X,p=rep(0.1,q),C=matrix(0,0,0),sp=sp,y=eta-offset,w=y*0+1) ###?????????
          M$Ain<-matrix(0,q,q); diag(M$Ain)<-rep(1,q);
          M$bin<-rep(-1e+12,q); M$bin[iv]<-1e-12
          M$off<-rep(0,n.pen); M$S<-list()
          for (j in 1:n.pen) {M$S[[j]]<-matrix(0,q,q);M$S[[j]]<- S[[j]]}
          beta.t<-pcls(M)      # initialize model coefficients (beta tilde)
          beta<-beta.t         # initialize beta
          beta[iv]<-log(beta.t[iv]) # values of beta of the monotone terms
      }
      else {
          # beta<-beta0+dbeta0*(log(sp)-log(sp.old))
          beta<-beta0
          beta.t<-beta               # current beta tilde
          beta.t[iv]<-exp(beta[iv])  # values of beta tilde of the monotone term
      }
      ## Initialization of parameters finishes here 
      #-------------------------------------------
      eta <- X%*%beta.t + offset     # define initial linear predictor
      mu <- linkinv(eta)  # define initial fitted model
      # end of the initialization ------------------------------------------------
      # --------------------------------------------------------------------------
      dev <- sum(dev.resids(y,mu,weights)) # define initial norm
      pdev<- dev + sum((rS%*%beta)^2) # define initial norm + penalty
      old.pdev <- pdev       # initialize convergence control for norm + penalty
      pdev.plot<-0     # define initial norm + penalty for plotting it 
      C <- diag(1,q,q)     # define matrix C
      C1 <- matrix(0,q,q)  # define matrix C1
      E <- matrix(0,q,q)   # define diagonal matrix E- second term of the Hessian
      C[ia] <- beta.t[iv]  # diagonal elements of C corresponding to the monotone terms
      C1[ia] <- beta.t[iv] # diagonal elements of C1 corresponding to the monotone terms
      X1 <- X%*%C 
      g.deriv <- 1/mu.eta(eta)            # diag(G)
      w1 <- weights/(variance(mu)*g.deriv^2)    # diag(W1)
      Dp.g <- -t(X1)%*%(w1*g.deriv*(y-mu))+S.t%*%beta # the gradient vector of the penalized deviance
      Dp.gnorm <- max(abs(Dp.g)) # set convergence on the max value of the Dp.g
      old.beta <- beta
      conv <- FALSE
     # --------------------------------------------------------------------------
     # MAIN ITERATIONS START HERE ----------------------------------------------
     ## while (abs(pnorm-old.pnorm)>1e-10*pnorm)  # repeat un-converged
     ##  while (Dp.gnorm>1e-8)  # repeat un-converged
     for (iter in 1:maxit)  {
         ## ADDED CODE...... 
         good <- weights > 0
         var.val <- variance(mu)
         varmu <- var.val[good]
         if (any(is.na(varmu))) 
              stop("NAs in V(mu)")
         if (any(varmu == 0)) 
              stop("0s in V(mu)")
         mu.eta.val <- mu.eta(eta)
         if (any(is.na(mu.eta.val[good]))) 
              stop("NAs in d(mu)/d(eta)")
         good <- (weights > 0) & (mu.eta.val != 0)
         if (all(!good)) {
              conv <- FALSE
              warning("No observations informative at iteration ", 
                  iter)
              break
         }
         ## END of ADDED CODE...... 
         C[ia] <- beta.t[iv]  # diagonal elements of C corresponding to the monotone terms
         C1[ia] <- beta.t[iv] # diagonal elements of C1 corresponding to the monotone terms
         X1 <- X%*%C 
         g.deriv<-1/mu.eta(eta)            # diag(G)
         w1 <- weights/(variance(mu)*g.deriv^2)    # diag(W1)
         alpha<-1+(y-mu)*(dv$dvar(mu)/variance(mu)+dg$d2link(mu)/g.deriv) # alpha elements of W
         w<-w1*alpha          # diag(W)
         diag(E)<-t(X%*%C1)%*%(w1*g.deriv*(y-mu)) # diagonal elements of E
         abs.w<-abs(w)      # absolute values of the diag(W)
         I.minus<-rep(0,nobs)  # define diagonal elements of the matrix I^{-}
         z1<-g.deriv*(y-mu)/alpha  # the first term of the pseudodata
         for (i in 1:nobs) {
             if (w[i]<0) {
                  I.minus[i]<-1; z1[i]<--z1[i]}
         }
         wX1<-sqrt(abs.w)[1:nobs]*X1
         wX11<-rbind(wX1,rS) # augmented model matrix 
         if (SVD) {
             # svd.s<-svd(wX11)
             svd.s <- try(svd(wX11),silent = TRUE) 
             if (inherits(svd.s, "try-error"))
                   svd.s <- try(svd(wX11,LINPACK = TRUE),silent = TRUE) 
             d.inv<-rep(0,q) # initial vector of inverse singular values
             for (i in 1:q) {       # corrected vector of inverse singular values 
                if (svd.s$d[i] < max(svd.s$d)*sqrt(.Machine$double.eps)) {
                       d.inv[i]<-0}   
                else d.inv[i]<-1/svd.s$d[i]
             }
             good <- d.inv > 0
             Q <- svd.s$u[,good]
             R <- svd.s$d[good]*t(svd.s$v[,good])
             R.inv<-svd.s$v[,good]%*%diag(d.inv[good])# inverse of R 
             tR.inv<-d.inv[good]*t(svd.s$v[,good])    # inverse of transpose of R
         } ### end if (SVD)
         else {
             qrst<-qr(wX11)        # QR decomposition of the augmented model matrix
             R<-qr.R(qrst)         # matrix R from QR decomposition
             Q<-qr.Q(qrst)         # matrix Q from QR decomposition
             tR.inv<-backsolve(t(R),diag(1,q),upper.tri=FALSE)  # inverse of transpose of R
             R.inv<-backsolve(R,diag(1,q))      # inverse of R 
         }
         QtQRER <- 2*t(Q[1:nobs,])%*%(I.minus*Q[1:nobs,])+tR.inv%*%(diag(E)*R.inv) # term in the Hessian 
         ei <- eigen(QtQRER,symmetric=TRUE)    # eigen-decomposition of QtQRER
         d <- ei$values        # vector of eigenvalues
         ok1 <- FALSE
         for (i in 1:length(d)) {     # checking positive semi-definiteness      
               if (d[i]>1) ok1 <- TRUE}
         if (ok1==TRUE) {# Fisher step in case of not positive semi-definiteness 
                  # of penalized loglikelihood
                  # set alpha =1
             eta.t<-X1%*%beta     # eta tilde for pseudodata 
             wX1<-sqrt(w1)[1:nobs]*X1
             wX11<-rbind(wX1,rS)
             z<-g.deriv*(y-mu)+eta.t      # pseudodata
             wz<-sqrt(w1)*z               # weighted pseudodata
             wz.aug<-c(wz,rep(0,nrow(rS)))   # augmented pseudodata
             s1 <- try(svd(wX11),silent=TRUE)   # SVD of the augmented model matrix
             if (inherits(s1, "try-error"))
                   s1 <- try(svd(wX11,LINPACK = TRUE),silent=TRUE)
             d.inv1 <- rep(0,q)  # initial vector of inverse singular values  
             for (i in 1:q) { # corrected vector of sv and inverse sv 
                 if (s1$d[i] < max(s1$d)*(.Machine$double.eps)^(1/2))
                         d.inv1[i] <- 0   
                 else d.inv1[i] <- 1/s1$d[i]
             }  
             good1 <- d.inv1 >0 
             ## old.beta<-beta  # store for 'step reduction'
             beta<-s1$v[,good1]%*%(d.inv1[good1]*t(s1$u[,good1]))%*%wz.aug 
         }  ### end of if (ok1) - Fisher step 
         else  {       ##  full Newton step
             Id.inv.r<-1/(1-d)^(1/2)   # vector of inverse values of (1-sv)^1/2
             for (i in 1:length(d)){
                 if ((1-d[i]) < .Machine$double.eps) Id.inv.r[i] <- 0
             }
             P<-R.inv%*%ei$vectors%*%diag(Id.inv.r) # define matrix P
             K<-Q[1:nobs,]%*%(ei$vectors%*%diag(Id.inv.r))  # define matrix K 
             wz1<-sqrt(abs.w)*z1     # the first term of the weighted pseudodata
             ## old.beta<-beta  # store for 'step reduction'
             beta<-c(old.beta)+P%*%(t(K)%*%wz1)-P%*%(t(P)%*%(S.t%*%old.beta))
         }  ###  end of if (!ok1) - Newton step
         delta<-beta-c(old.beta)         # trial step
         step<-1                      # initial trial step length
         beta<-c(old.beta)+step*delta    # current parameter estimates
         beta.t<-beta                 # current beta tilde
         beta.t[iv]<-exp(beta[iv])  # values of beta tilde of the monotone term
         eta<-X%*%beta.t + offset     # linear predictor
         mu<- linkinv(eta)          # fitted values
         dev<-sum(dev.resids(y,mu,weights)) # deviance of the working model
         ## old.pdev<-pdev                 # store for convergence test
         pdev<- dev + sum((rS%*%beta)^2) # deviance + penalty of the working model
               
         ## `step reduction' approach starts here ------------------------------------
         ii <- 1 
         div.thresh <- 10*(0.1 +abs(old.pdev))*.Machine$double.eps^0.5
         ##  while (!(valideta(eta) && validmu(mu))){ # `step reduction'
         while (is.na(pdev) || (pdev-old.pdev) > div.thresh) { # 'step reduction' approach
             if (ii > 200) 
               stop ("step reduction failed")
             ii <- ii+1
             step<-step/2                # decrease step length 
             beta<-c(old.beta)+step*delta   # update current parameter estimates
             beta.t<-beta                # update current beta tilde
             beta.t[iv]<-exp(beta[iv])   # values of beta tilde of the monotone term
             eta<-X%*%beta.t + offset    # linear predictor  
             mu<- linkinv(eta)         # fitted values
             dev <- sum(dev.resids(y,mu,weights))# update deviance of the working model
             pdev <- dev+sum((rS%*%beta)^2) # update deviance + penalty of the working model
         }
         ## `step reduction' finishes here -------------------------------------------
     
         Dp.g <- -t(X1)%*%(w1*g.deriv*(y-mu))+S.t%*%beta # the gradient vector of the penalized deviance
         Dp.gnorm<-max(abs(Dp.g)) 
         pdev.plot[iter]<-pdev      # store deviance+penalty of the working model for plotting
          
         ## cheking convergence .......
         if (abs(pdev - old.pdev)/(0.1 + abs(pdev)) < epsilon) {
             if (max(abs(beta - c(old.beta))) > epsilon * 
                          max(abs(beta + c(old.beta)))/2) {
                old.beta <- beta
                old.pdev <- pdev
             }
             else {
                conv <- TRUE
                beta <- beta
                break
             }
         }
         else {
             old.pdev <- pdev
             old.beta <- beta
         }
     } ### main iteration procedure is completed here --------------------------------
     # ---------------------------------------------------------------------------
  
     # ___________________________________________________________________________

     # ---------------------------------------------------------------------------
     # define matrices at their converged values from the full Newtom method------
  
     dev <- sum(dev.resids(y,mu,weights))
     beta.t <- beta                # estimates of beta tilde
     beta.t[iv]<- exp(beta[iv])   # values of beta tilde of the monotone term
     eta <- X%*%beta.t + offset      # linear predictor  
     mu <- linkinv(eta)         # fitted values
     C[ia] <- beta.t[iv]  # diagonal elements of C corresponding to the monotone terms
     C1[ia] <- beta.t[iv] # diagonal elements of C1 corresponding to the monotone terms
     X1 <- X%*%C 
     g.deriv <- 1/ mu.eta(eta)        # diag(G)
     w1 <- weights/(variance(mu)*g.deriv^2)    # diag(W1)
     alpha <- 1+(y-mu)*(dv$dvar(mu)/variance(mu)+dg$d2link(mu)/g.deriv) # alpha elements of W
     w <- w1*alpha          # diag(W)
     diag(E) <- t(X%*%C1)%*%(w1*g.deriv*(y-mu)) # diagonal elements of E
     abs.w <- abs(w)      # absolute values of the diag(W)
     I.minus <- rep(0,nobs)  # define diagonal elements of the matrix I^{-}
     for (i in 1:nobs) {
         if (w[i] < 0) I.minus[i]<-1
     }
     wX1 <- sqrt(abs.w)[1:nobs]*X1
     wX11 <- rbind(wX1,rS) # augmented model matrix 
     if (SVD) {
        # svd.s<-svd(wX11)
        svd.s <- try(svd(wX11),silent = TRUE) 
        if (inherits(svd.s, "try-error"))
                 svd.s <- try(svd(wX11,LINPACK = TRUE),silent = TRUE) 
        d.inv <- rep(0,q)        # initial vector of inverse singular values
        illcond <- FALSE
        for (i in 1:q) {  # corrected vector of inverse singular values 
            if (svd.s$d[i] < max(svd.s$d)*sqrt(.Machine$double.eps))
                 {d.inv[i] <- 0; illcond <- TRUE}   
            else d.inv[i ]<- 1/svd.s$d[i]
        }
        good <- d.inv > 0
        R <- svd.s$d[good]*t(svd.s$v[,good])
        Q <- svd.s$u[,good]
        R.inv <- svd.s$v[,good]%*%diag(d.inv[good])# inverse of R 
        tR.inv <- d.inv[good]*t(svd.s$v[,good])    # inverse of transpose of R
     }  ### end of if (SVD)
     else  {
        qrst<-qr(wX11)        # QR decomposition of the augmented model matrix
        R<-qr.R(qrst)         # matrix R from QR decomposition
        Q<-qr.Q(qrst)         # matrix Q from QR decomposition
        tR.inv<-backsolve(t(R),diag(1,q),upper.tri=FALSE)  # inverse of transpose of R
        R.inv<-backsolve(R,diag(1,q))      # inverse of R 
     }  ### end of if (!SVD)
     QtQRER<-2*t(Q[1:nobs,])%*%(I.minus*Q[1:nobs,])+tR.inv%*%(diag(E)*R.inv) # term in the Hessian 
     ei<-eigen(QtQRER,symmetric=TRUE)    # eigen-decomposition of QtQRER
     d<-ei$values        # vector of eigenvalues
     ok1 <- FALSE
     for (i in 1:length(d))  {    # checking positive semi-definiteness      
         if (d[i]>1) ok1 <- TRUE}
     if (ok1) { ## Fisher step in case of not positive semi-definiteness 
                ## of penalized loglikelihood
                ## set alpha =1 ...
         eta.t<-X1%*%beta  # eta tilde for pseudodata 
         wX1<-sqrt(w1)[1:nobs]*X1
         wX11<-rbind(wX1,rS)
         s1 <- try(svd(wX11),silent=TRUE)   # SVD of the augmented model matrix
         if (inherits(s1, "try-error"))
                 s1 <- try(svd(wX11,LINPACK = TRUE),silent=TRUE)
         d.inv1<-c(0,q)  # initial vector of inverse singular values  
         for (i in 1:q) {  # corrected vector of sv and inverse sv 
             if (s1$d[i] < max(s1$d)*(.Machine$double.eps)^(1/2))
                     {d.inv1[i]<-0}   
             else d.inv1[i]<-1/s1$d[i]
         }
         good1 <- d.inv1 > 0
         P <- s1$v[,good1]%*%diag(d.inv1[good1])  # define matrix P
         K <- s1$u[1:nobs,good1]   # define matrix K
     } ## end of if (ok1)
     else  {   ## full Newton step
         Id.inv.r<-1/(1-d)^(1/2)   # vector of inverse values of (1-sv)^1/2
         for (i in 1:length(d)){
              if ((1-d[i]) < .Machine$double.eps) Id.inv.r[i] <- 0
         }
         P <- R.inv%*%ei$vectors%*%diag(Id.inv.r) # define matrix P
         K <- Q[1:nobs,]%*%(ei$vectors%*%diag(Id.inv.r))  # define matrix K 
     }  ## end of if (!ok1)
     # end of calculation of the matrices at their converged values --------------
     # ---------------------------------------------------------------------------
     Dp.g <- -t(X1)%*%(w1*g.deriv*(y-mu))+S.t%*%beta # the gradient vector of the penalized deviance
     Dp.gnorm<-max(abs(Dp.g)) 

     # calculating tr(A) ---------------------------------------------------------
     I.plus<-rep(1,nobs)   # define diagonal elements of the matrix I^{+}
     for (i in 1:nobs) {
          if (w[i]<0) I.plus[i]<--1}
     L<-c(1/alpha)    # define diagonal elements of L=diag(1/alpha)
     PKt <- P%*%t(K)   
     edf <- rowSums(PKt*t(L*I.plus*wX1)) # edf for each model parameter
     trA <- sum(edf)
     # ---------------------------------------------------------------------------
     scale.est <- dev/(nobs-trA) #  scale estimate...
     residuals <- rep.int(NA, nobs)
     residuals <- (y-mu)*g.deriv
     #  good <- weights > 0
     #  mug <- mu[good];yg <- y[good]
     #  weg <- weights[good]
     #  residuals[good] <- (yg - mug)*g.deriv
     # ---------------------------------------------------------------------------
     ## calculation of the derivatives of beta by the Implicit Function Theorem starts here
     dbeta.rho<-matrix(0,q,n.pen) # define matrix of the parameters derivatives
     PPt<-P%*%t(P)  # inverse matrix of the Hessian
     for (j in 1:n.pen) {
        dbeta.rho[,j]<--sp[j]*PPt%*%(S[[j]]%*%beta) # derivative of beta wrt rho[j]
     }
     # end of calculating the parameters derivatives
     aic.model <- aic(y, n, mu, weights, dev) +  2 * sum(edf)
     assign("start",beta,envir=ee)
     assign("dbeta.start",dbeta.rho,envir=eb)
     assign("sp.last",sp,envir=esp)
  } ### end if (!EMPTY) 

 list(L=L,C1=C1,E=E,iter=iter, old.beta=old.beta, step=step,trA=trA,gcv=dev*nobs/(nobs-trA)^2,
      sp=sp, mu=mu,X=X, X1=X1,beta=beta,beta.t=beta.t,iv=iv,S=S,S.t=S.t,rS=rS,
      d.inv=d.inv, P=P,K=K,dlink.mu=1/mu.eta(eta),Var=variance(mu), abs.w=abs.w,
      link=family$linkfun(mu),w=w,w1=w1,d2link.mu=dg$d2link(mu),wX1=wX1,I.plus=I.plus,
      dvar.mu=dv$dvar(mu),d2var.mu=dv$d2var(mu),deviance=dev,scale.est=scale.est,
      ok1=ok1,alpha=alpha,d3link.mu=dg$d3link(mu),eta=eta, edf=edf,iter=iter,
      svd.d=svd.s$d,Dp.gnorm=Dp.gnorm, Dp.g=Dp.g,d=d, conv=conv, illcond=illcond,
      residuals=residuals,z=g.deriv*(y-mu)+X1%*%beta,dbeta.rho=dbeta.rho, aic=aic.model)
}



## ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


scam.fit.post <- function(y,X,object,sig2,offset=rep(0, nrow(X)),
                    intercept=TRUE,weights=rep.int(1,nrow(X)),scale.known)
 {
## get null deviance and covariance matrices after a scam fit. 
  n <- nrow(X) # number of observations
  linkinv <- object$family$linkinv
  wtdmu <- if (intercept) 
        sum(weights * y)/sum(weights)
  else linkinv(offset)
  dev.resids <- object$family$dev.resids
  nulldev <- sum(dev.resids(y, wtdmu, weights)) # deviance for single parameter model
  
  # calculating the approximate covariance matrices 
  # (dealing with the expected Hessian of the log likelihood) ...
#  if (sig2 > 0) scale.known <- TRUE else scale.known <- FALSE
  if (!scale.known) sig2 <- object$scale.est
  
  ## get the inverse of the expected Hessian...
  wX1 <- sqrt(object$w1)[1:n]*object$X1
  wX11 <- rbind(wX1,object$rS)
  s1 <- try(svd(wX11),silent=TRUE)   # SVD of the augmented model matrix
  if (inherits(s1, "try-error"))
        s1 <- try(svd(wX11,LINPACK = TRUE),silent=TRUE)
  q <- ncol(X)
  d.inv1 <- rep(0,q)  # initial vector of inverse singular values  
  for (i in 1:q) {  # corrected vector of sv and inverse sv 
       if (s1$d[i] < max(s1$d)*(.Machine$double.eps)^(1/2))
                     {d.inv1[i] <- 0}   
       else d.inv1[i] <- 1/s1$d[i]
  }
  good <- d.inv1 > 0
  P <- s1$v[,good]%*%diag(d.inv1[good])  
  K <- s1$u[1:n,good]   
  PKt <- P%*%t(K)
  Vb <- P%*%t(P)*sig2 ## object$scale.est  # Bayesian posterior covariance matrix for the parameters 
  Ve <- PKt%*%t(PKt)*sig2 # covariance matrix of the parametesr estimators 
  
  # Delta method to get covariance matrix for the reparametrized parameters...
    df.p <- diag(rep(1,q)) # diag matrix of parameters derivatives
    diag(df.p)[object$iv]<- object$beta.t[object$iv] 
    Vb.t <- (diag(df.p)*Vb)%*%df.p
    Ve.t <- (diag(df.p)*Ve)%*%df.p
 
#  PKt0 <- object$P%*%t(object$K)   
#  edf <- rowSums(PKt0*t(object$L*object$I.plus*object$wX1))
#  edf <- rowSums(PKt*t(sqrt(object$weights)*X))
#  hat <- rowSums(object$K*object$K)
  list (nulldev=nulldev, Vb=Vb,Vb.t=Vb.t,Ve=Ve,Ve.t=Ve.t,sig2=sig2)
}


