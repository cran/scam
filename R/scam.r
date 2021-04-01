
#############################################################
## the wrapper overall function to fit scam...             ##
#############################################################

scam <- function(formula, family=gaussian(), data=list(), gamma=1, sp=NULL, 
                 weights=NULL, offset=NULL, optimizer="bfgs", optim.method=c("Nelder-Mead","fd"),
                 scale=0, knots=NULL, not.exp=FALSE, start=NULL, etastart, mustart, control=list(),
                 AR1.rho=0, AR.start=NULL,drop.unused.levels=TRUE) ##,devtol.fit=1e-8, steptol.fit=1e-8, check.analytical=FALSE, del=1e-4)
{  ## scale - scale parameter of the exponential deistribution as in gam(mgcv)
   ## optimizer - numerical optimization method to use to optimize the smoothing parameter estimation criterion: "bfgs", "optim", "nlm", "nlm.fd", "efs"
   ## optim.method - if optimizer=="optim" then the first argument of optim.method     specifies the method, and the second can be either "fd" for finite-difference approximation of the gradient or "grad" - to use analytical gradient of gcv/ubre
   ## not.exp - if TRUE then notExp() function will be used in place of exp in positivity ensuring beta parameters re-parameterization
   
   control <- do.call("scam.control",control)
   ## Setting from mgcv(gam).......
   ## create model frame..... 
   gp <- interpret.gam(formula) # interpret the formula 
   cl <- match.call() # call needed in gam object for update to work
   mf <- match.call(expand.dots=FALSE)
   mf$formula <- gp$fake.formula 
   mf$family <- mf$control<-mf$scale<-mf$knots<-mf$sp<-mf$min.sp<-mf$H<-mf$select <- mf$drop.intercept <-
                mf$gamma<-mf$method<-mf$fit<-mf$paraPen<-mf$G<-mf$optimizer <- mf$optim.method <- mf$not.exp <- mf$in.out <- mf$AR1.rho <- mf$devtol.fit <- mf$steptol.fit <- mf$del <- mf$...<-NULL
   mf$drop.unused.levels <- drop.unused.levels
   mf[[1]] <- quote(stats::model.frame) ## as.name("model.frame")
   pmf <- mf
   mf <- eval(mf, parent.frame()) # the model frame now contains all the data 
   if (nrow(mf)<2) stop("Not enough (non-NA) data to do anything meaningful")
   terms <- attr(mf,"terms") 

   ## summarize the *raw* input variables
   ## note can't use get_all_vars here -- buggy with matrices
   vars <- all.vars1(gp$fake.formula[-2]) ## drop response here
   inp <- parse(text = paste("list(", paste(vars, collapse = ","),")"))

   ## allow a bit of extra flexibility in what `data' is allowed to be (as model.frame actually does)
   if (!is.list(data)&&!is.data.frame(data)) data <- as.data.frame(data) 

   dl <- eval(inp, data, parent.frame())
   names(dl) <- vars ## list of all variables needed
   var.summary <- variable.summary(gp$pf,dl,nrow(mf)) ## summarize the input data
   rm(dl) ## save space    

   ## pterms are terms objects for the parametric model components used in 
   ## model setup - don't try obtaining by evaluating pf in mf - doesn't
   ## work in general (e.g. with offset)...

   if (is.list(formula)) { ## then there are several linear predictors
      environment(formula) <- environment(formula[[1]]) ## e.g. termplots needs this
      pterms <- list()
      tlab <- rep("",0)
      for (i in 1:length(formula)) {
        pmf$formula <- gp[[i]]$pf 
        pterms[[i]] <- attr(eval(pmf, parent.frame()),"terms")
        tlabi <- attr(pterms[[i]],"term.labels")
        if (i>1&&length(tlabi)>0) tlabi <- paste(tlabi,i-1,sep=".")
        tlab <- c(tlab,tlabi)
     }
     attr(pterms,"term.labels") <- tlab ## labels for all parametric terms, distinguished by predictor
   } else { ## single linear predictor case
      pmf$formula <- gp$pf
      pmf <- eval(pmf, parent.frame()) # pmf contains all data for parametric part
      pterms <- attr(pmf,"terms") ## pmf only used for this
   }

   if (is.character(family)) family <- eval(parse(text=family))
   if (is.function(family)) family <- family()
   if (is.null(family$family)) stop("family not recognized")
  
   if (family$family[1]=="gaussian" && family$link=="identity") am <- TRUE
   else am <- FALSE

   if (AR1.rho!=0&&!is.null(mf$"(AR.start)")) if (!is.logical(mf$"(AR.start)")) stop("AR.start must be logical")
   if (AR1.rho!=0 && !am) stop("residual autocorrelation, AR1, is currently available only for the Gaussian identity link model.")
    
   if (!control$keepData) rm(data) ## save space

   ## check whether family requires intercept to be dropped...
   # drop.intercept <- if (is.null(family$drop.intercept) || !family$drop.intercept) FALSE else TRUE
   # drop.intercept <- as.logical(family$drop.intercept)
  ### if (is.null(family$drop.intercept)) { ## family does not provide information
  ###    lengthf <- if (is.list(formula)) length(formula) else 1
  ###    if (is.null(drop.intercept)) drop.intercept <- rep(FALSE, lengthf) else {
  ###      drop.intercept <- rep(drop.intercept,length=lengthf) ## force drop.intercept to correct length
  ###	if (sum(drop.intercept)) family$drop.intercept <- drop.intercept ## ensure prediction works
  ###    }
  ### } else drop.intercept <- as.logical(family$drop.intercept) ## family overrides argument
    
 ### if (inherits(family,"general.family")&&!is.null(family$presetup)) eval(family$presetup)

 ###   gsname <- if (is.list(formula)) "gam.setup.list" else "gam.setup" 
 ###   G <- do.call(gsname,list(formula=gp,pterms=pterms,
 ###                data=mf,knots=knots,sp=sp, min.sp=min.sp,
 ###                H=H,absorb.cons=TRUE,sparse.cons=0,select=select,
 ###                idLinksBases=control$idLinksBases,scale.penalty=control$scalePenalty,
 ###                paraPen=paraPen,drop.intercept=drop.intercept))
   G <- do.call("gam.setup",list(formula=gp,pterms=pterms,
                 data=mf,knots=knots,sp=sp, absorb.cons=TRUE,sparse.cons=0))  
    
   G$var.summary <- var.summary
   G$family <- family
   
   if ((is.list(formula)&&(is.null(family$nlp)||family$nlp!=gp$nlp))||
        (!is.list(formula)&&!is.null(family$npl)&&(family$npl>1))) stop("incorrect number of linear predictors for family")
       
   G$terms<-terms;
   G$mf<-mf;G$cl<-cl;
   G$am <- am
   G$AR1.rho <- AR1.rho; G$AR.start <- AR.start

   if (is.null(G$offset)) G$offset<-rep(0,G$n)
     
   G$min.edf <- G$nsdf ## -dim(G$C)[1]
   if (G$m) for (i in 1:G$m) G$min.edf<-G$min.edf+G$smooth[[i]]$null.space.dim

   G$formula <- formula
   G$pred.formula <- gp$pred.formula
   environment(G$formula)<-environment(formula)
  
  if (ncol(G$X)>nrow(G$X)) stop("Model has more coefficients than data") 

###  G <- gam(formula, family,data=data, knots=knots,fit=FALSE) 
###   n.terms <- length(G$smooth)  ## number of smooth terms in the model
###   n <- nrow(G$X)
##   intercept <- G$intercept ## TRUE or FALSE
   ## now need to set 'offset' as the above G wouldn't take in 'offset' that is outside of formula..
###   gp <- interpret.gam(formula) # interpret the formula 
###   cl <- match.call() # call needed in gam object for update to work
###   mf <- match.call(expand.dots=FALSE)
###   mf$formula <- gp$fake.formula 
###   mf$family <- mf$control<-mf$scale<-mf$knots<-mf$sp<-mf$min.sp<-mf$H<-mf$select <-
###                 mf$gamma<-mf$method<-mf$fit<-mf$paraPen<-mf$G<-mf$optimizer <- mf$optim.method <- mf$not.exp <- mf$in.out <- 
###                  mf$devtol.fit <- mf$steptol.fit <- mf$del <- mf$...<-NULL
###   mf[[1]] <- quote(stats::model.frame) ## as.name("model.frame")
###   pmf <- mf
###   mf <- eval(mf, parent.frame()) # the model frame now contains all the data 
###   if (nrow(mf)<2) stop("Not enough (non-NA) data to do anything meaningful")
   # terms <- attr(mf,"terms") 
  
   n.terms <- length(G$smooth)  ## number of smooth terms in the model
   n <- nrow(G$X)
   intercept <- G$intercept

   G$offset <- as.vector(model.offset(mf)) 
   if (is.null(G$offset)) 
        G$offset <- rep.int(0, n)
   ## done offset 
   
 ###  if (is.null(weights)) 
 ###        weights <- rep.int(1,n)
   weights <- G$w
   

   fam.name <- G$family[1]
   if (scale == 0) 
      {  if (fam.name == "binomial" || fam.name == "poisson") 
               sig2 <- 1
         else sig2 <- -1
      }
   else { sig2 <- scale }
   if (sig2 > 0) scale.known <- TRUE else scale.known <- FALSE
   ## get penalty matrices and 
   ## vector of identifications for exponentiated model coefficients...
   Q <- penalty_pident(G)
   ## checking sp...
   if (!is.null(sp)) {
         neg <- FALSE
         if (length(sp)!= length(G$off)) {
              warning("Supplied smoothing parameter vector is too short - ignored.")
              sp <- NULL       
         } else if (sum(is.na(sp))) {
               warning("NA's in supplied smoothing parameter vector - ignoring.")
               sp <- NULL
           } else {
                good <- sp < 0
                if (sum(good) > 0) { ## cheking negative values..
                   warning("Supplied smoothing parameter vector has negative values - ignored.")
                   neg <- TRUE
                }
             }                
         if (neg) sp <- NULL
     }
   ## Create new environments with `start' initially empty...
 #  ee <- new.env() ## parent = .GlobalEnv
   env <- new.env() 
   assign("start",rep(0,0),envir=env)
   assign("dbeta.start",rep(0,0),envir=env)
   assign("sp.last",rep(0,0),envir=env)
   
   q.f <- rep(0,n.terms)
   if (n.terms >0) for (i in 1:n.terms) 
                      q.f[i] <- ncol(G$smooth[[i]]$S[[1]]) + 1 
   G$S <- Q$S
   G$q.f <- q.f
   G$q0 <- G$off[1]-1  ## number of the parameters of the strictly parametric model
   G$p.ident <- Q$p.ident  ## vector of 0's & 1's for the model parameters identification: 
   G$n.terms <- n.terms   ## number of the smooth terms in the SCAM
  # G$intercept <- intercept
   G$weights <- weights
   G$sig2 <- sig2
   G$scale.known <- scale.known
   G$not.exp <- not.exp

  ### if (!control$keepData) rm(data) ## save space

   object <- list() 
   if (is.null(sp)) { 
       ## get initial estimates of the smoothing parameter...
         start <- etastart <- mustart <- NULL
         y <- G$y; family <- G$family
         nobs <- NROW(y)
         eval(family$initialize)
         G$y <- y  ## needed to set factor response values as numeric
         def.sp <- initial.sp.scam(G, Q,q.f=q.f,n.terms=n.terms,family=family,
              intercept=intercept,offset=G$offset, env=env,
              weights=weights, control=control) 
         rho <- log(def.sp+1e-4) ## get initial log(sp) ...
         ## minimize GCV/UBRE by optimizer....
         ptm <- proc.time()
         re <- estimate.scam(G=G,optimizer=optimizer,optim.method=optim.method,
               rho=rho, gamma=gamma, env=env,control=control)  
         CPU.time <- proc.time()-ptm
         best <- re
         object$gcv.ubre <- re$gcv.ubre
         object$dgcv.ubre <- re$dgcv.ubre
     #    object$aic <- re$aic
         best$p.ident <- Q$p.ident
         best$S <- Q$S
         object$optimizer <- optimizer
         object$edf1 <- re$edf1
         object$termcode <- re$termcode
         if (optimizer == "bfgs")
            {   object$check.grad <- re$check.grad
                object$dgcv.ubre.check <- re$dgcv.ubre.check
            }
   } else {   ## no GCV minimization if sp is given...
            best <- scam.fit(G=G, sp=sp,gamma=gamma,env=env, control=control) 
            object$optimizer <- "NA"           
      }
   ## post-fitting values...
   best$n.smooth <- object$n.smooth <- n.terms
   best$formula <- object$formula <- formula
   best$family <- object$family <- G$family
   best$smooth <- object$smooth <- G$smooth
   best$model <- object$model <- G$mf

   object$R <- best$R
   if (is.null(object$R)){
         rr <- scam.fit(G=G, sp=best$sp, gamma=gamma,env=env, control=control) 
         object$R <- rr$R } ## not sure if it's needed?
  
  # object$df.residual <- nrow(best$X) - sum(best$edf)

   object$sp <- best$sp
   names(object$sp) <- names(G$sp)
   if (sum(is.na(names(object$sp)))!=0){  ## create names for sp if NA's from G
      if (n.terms >0) for (i in 1:n.terms) names(object$sp)[i] <- object$smooth[[i]]$label
   }
  # object$deviance <- best$deviance
  # object$residuals <- best$residuals
#   object$X <- best$X

   object$conv <- best$conv # whether or not the inner full Newton method converged
   post <- scam.fit.post(G, object=best) ##,sig2=sig2,offset = G$offset,
                 ##  intercept=G$intercept, weights=weights, scale.known=scale.known,not.exp=not.exp) 

   object$edf <- best$edf  # post$edf 
   object$edf1 <- post$edf1
   object$trA <- best$trA  # post$trA
   names(object$edf) <- G$term.names
   names(object$edf1) <- G$term.names
   object$aic <- post$aic  # best$aic
   object$null.deviance <- post$null.dev
   object$deviance <- post$deviance
   object$residuals <- post$residuals
   object$df.residual <- nrow(G$X) - sum(post$edf)
   object$rank <- post$rank

   object$var.summary <- G$var.summary 
   object$cmX <- G$cmX ## column means of model matrix --- useful for CIs
   object$model<-G$mf # store the model frame
   
   object$full.sp <- G$full.sp  ## incorrect !!!
   if (!is.null(object$full.sp))   names(object$full.sp) <- names(G$full.sp)

   object$na.action <- attr(G$mf,"na.action") # how to deal with NA's
   object$df.null <- post$df.null
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
   object$contrasts <- G$contrasts
   object$xlevels <- G$xlevels
   object$nsdf <- G$nsdf
   object$y <- G$y
 #  object$data <- G$mf
   if (control$keepData) object$data <- data 
   object$control <- control
   object$offset <- G$offset
   object$not.exp <- G$not.exp
 #  object$scale.known <- scale.known # to be passed in the summary function
   object$scale.estimated <- !scale.known # to be passed in the summary function
   object$prior.weights <-weights # prior weights on observations
   object$weights <- best$w  # final weights used in full Newton iteration
   object$fitted.values <- post$mu  #best$mu
   object$linear.predictors <- post$eta #best$eta
 #  cl<-match.call() # call needed in gam object for update to work
   object$call <- cl 
   object$p.ident <- Q$p.ident
   object$intercept <- G$intercept
   object$min.edf <- G$min.edf ## Minimum possible degrees of freedom for whole model
   object$gamma <- gamma
   object$iter <- best$iter  # number of iterations of the Full Newton
   if (is.null(sp)) 
         object$CPU.time <- CPU.time
   else 
         object$CPU.time <- NULL
   object$AR1.rho <- AR1.rho

   ## get the optimizer info (smoothing parameter selection).....
  if (is.null(sp))
     {   if (optimizer == "bfgs")
            {   ## get the bfgs info in case of sp selection... 
                object$bfgs.info <- list()
                object$bfgs.info$conv <- re$conv.bfgs  
                object$bfgs.info$iter <- re$iterations 
                object$bfgs.info$grad <- re$dgcv.ubre
                object$bfgs.info$score.hist <- re$score.hist
            }
         else if (optimizer == "nlm.fd" || optimizer == "nlm")
                 {   object$nlm.info <- list()
                     object$nlm.info$conv <- re$conv 
                     object$nlm.info$iter <- re$iterations 
                     object$nlm.info$grad <- re$dgcv.ubre
                 }
         else if (optimizer=="optim")
                 {   object$optim.info <- list()
                     object$optim.info$conv <- re$conv 
                     object$optim.info$iter <- re$iterations 
                     object$optim.method <- re$optim.method  
                 }
         else if (optimizer=="efs")
                 {   object$efs.info <- list()
                     object$efs.info$conv <- re$conv 
                     object$efs.info$iter <- re$iterations 
                     object$efs.info$score.hist <- re$score.hist
                 }        
      }
   if (scale.known)
         object$method <- "UBRE"
   else  object$method <- "GCV" 
   if (G$nsdf > 0) 
         term.names <- colnames(G$X)[1:G$nsdf]
   else term.names <- array("", 0)
   if (n.terms >0) 
        for (i in 1:n.terms) 
            {   k <- 1  
                for (j in G$smooth[[i]]$first.para:G$smooth[[i]]$last.para) 
                   {   term.names[j] <- paste(G$smooth[[i]]$label, ".", 
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
   class(object) <- c("scam","glm","lm")  
   rm(G)
   dev <- object$deviance
   if (AR1.rho!=0){
      object$std.rsd <- AR.resid(object$residuals,AR1.rho,object$model$"(AR.start)") ##standardised residuals for AR1 model (approximately uncorrelated under correct model)
      dev <- sum(object$std.rsd^2)
      object$deviance <- sum(object$residuals^2) ## ss of the working residuals of the fitted model
   }
   object$aic <- object$family$aic(object$y,1,object$fitted.values,object$prior.weights,dev)
   object$aic <- object$aic -
                2 * (length(object$y) - sum(sum(object$model[["(AR.start)"]])))*log(1/sqrt(1-AR1.rho^2)) + ## correction for AR
                2*sum(object$edf)
   object
}  ## end scam


##############################################################
## control function for scam (similar to gam.control(mgcv)) ##
##############################################################

scam.control <- function (maxit = 200, maxHalf=30, devtol.fit=1e-7, steptol.fit=1e-7,
                          keepData=FALSE,efs.lspmax=15,efs.tol=.1, nlm=list(),optim=list(),bfgs=list(),
                          trace =FALSE, print.warn=FALSE) 
# Control structure for a scam. 
# devtol.fit is the tolerance to use in the scam.fit call within each IRLS. 
# check.analytical - logical whether the analytical gradient of GCV/UBRE should be checked for bfgs method
# del - increment for finite differences when checking analytical gradients for bfgs method
# maxHalf is the number of step halvings to employ in bfgs_gcv.ubre search, before giving up on a search direction. 
# trace turns on or off some de-bugging information.
# print.warn =FALSE - when set to 'FALSE' turns off printing warning messages for step halving under non-finite exponentiated coefficients,  non-finite deviance and/or if mu or eta are out of bounds.
{   if (!is.numeric(devtol.fit) || devtol.fit <= 0) 
        stop("value of devtol.fit must be > 0")
    if (!is.numeric(steptol.fit) || devtol.fit <= 0) 
        stop("value of steptol.fit must be > 0")
    if (!is.numeric(maxit) || maxit <= 0) 
        stop("maximum number of iterations must be > 0")
    if (!is.numeric(maxHalf) || maxHalf <= 0) 
        stop("maximum number of step halving must be > 0")
    if (!is.logical(trace)) stop("trace must be logical")
    if (!is.logical(print.warn)) stop("print.warn must be logical")

    # work through nlm defaults
    if (is.null(nlm$ndigit)||nlm$ndigit<2) nlm$ndigit <- max(2,ceiling(-log10(1e-7)))
    nlm$ndigit <- round(nlm$ndigit)
    ndigit <- floor(-log10(.Machine$double.eps))
    if (nlm$ndigit>ndigit) nlm$ndigit <- ndigit
    if (is.null(nlm$gradtol)) nlm$gradtol <- 1e-6
    nlm$gradtol <- abs(nlm$gradtol)
    ## note that nlm will stop after hitting stepmax 5 consecutive times
    ## hence should not be set too small ... 
    if (is.null(nlm$stepmax)||nlm$stepmax==0) nlm$stepmax <- 2
    nlm$stepmax <- abs(nlm$stepmax)
    if (is.null(nlm$steptol)) nlm$steptol <- 1e-4
    nlm$steptol <- abs(nlm$steptol)
    if (is.null(nlm$iterlim)) nlm$iterlim <- 200
    nlm$iterlim <- abs(nlm$iterlim)
    ## Should be reset for a while anytime derivative code altered...
    if (is.null(nlm$check.analyticals)) nlm$check.analyticals <- FALSE
    nlm$check.analyticals <- as.logical(nlm$check.analyticals) 

    # and bfgs defaults
    if (is.null(bfgs$check.analytical)) bfgs$check.analytical <- FALSE
    if (is.null(bfgs$del)) bfgs$del <- 1e-4
    if (is.null(bfgs$steptol.bfgs)) bfgs$steptol.bfgs <- 1e-7

    if (is.null(bfgs$gradtol.bfgs)) bfgs$gradtol.bfgs <- 1e-06
    if (is.null(bfgs$maxNstep)) bfgs$maxNstep <- 5 ## the maximum allowable step length
    if (is.null(bfgs$maxHalf)) bfgs$maxHalf <- maxHalf ## the maximum number of step halving in "backtracking"
 
    # and optim defaults
    if (is.null(optim$factr)) optim$factr <- 1e7
    optim$factr <- abs(optim$factr)
    if (efs.tol<=0) efs.tol <- .1
    
    list(maxit=maxit, devtol.fit=devtol.fit, steptol.fit=steptol.fit, keepData=as.logical(keepData[1]), nlm=nlm, optim=optim,bfgs=bfgs,efs.lspmax=efs.lspmax,efs.tol=efs.tol,trace = trace, print.warn=print.warn)    
} ## end scam.control


#################################################################
## function to get initial estimates of smoothing parameters...##
#################################################################

initial.sp.scam <- function(G,Q,q.f,n.terms,family,intercept,offset, env= env,
                      weights, control=control) 
{  ## function to get initial estimates of smoothing parameters
   control$devtol.fit <- 1e-4
   control$steptol.fit <- 1e-4
   ## step 1: set sp=rep(0.5,p) and estimate hessian... 
   b <- scam.fit(G=G,sp=rep(0.05,length(G$off)), env=env, control=control) 
   H <- crossprod(b$wX1) - b$E
   ## step 2:...
   n.p <- length(Q$S) ## number of penalty matrices
   def.sp <- array(0,n.p) ## initialize the initial sp values
   j <- 1
   if (n.terms >0)  for (i in 1:n.terms){
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
   env <- new.env()
   assign("start",rep(0,0),envir=env)
   assign("dbeta.start",rep(0,0),envir=env)
   assign("sp.last",rep(0,0),envir=env)
   def.sp
}


#########################################################
## function to get list of penalty matrices and        ## 
## vector of parameter identifications .....           ##
#########################################################


penalty_pident <- function(object)
{  ## function to get the list of penalties and vector of model parameters 
   ## identifications from the gam() setting...
   n.terms <- length(object$smooth)  # number of terms in the model
   q <- ncol(object$X)          # total number of parameters
   cons.terms <- rep(0,n.terms) # define whether each term is constrained or not
   if (n.terms>0) for (i in 1:n.terms){
                     if (!is.null(object$smooth[[i]]$p.ident))
                     cons.terms[i] <- 1  
                  }
   p.ident <- rep(FALSE,q) # initialize vector of parameter identifications
                      # with `TRUE' - for a parameter to be exponentiated, `FALSE' - otherwise
   off.terms <- rep(0,n.terms) # starting points for each term
   off <- object$off
   if (n.terms ==length(off))
          off.terms <- off
   else 
      {   off.terms[1] <- off[1]
          k <- 1
          l <- 1
          while (l<length(off))
              {   if (off[l]!=off[l+1])
                     {   off.terms[k+1] <- off[l+1] 
                         k <- k+1; l <- l+1 
                     } 
                  else l <- l+1
              }
     
      }
   if (n.terms>0) for (i in 1:n.terms){
                    if (cons.terms[i]==1) 
                       # if (inherits(object$smooth[[i]], c("mipc.smooth")))
                        #    p.ident[off.terms[i]:(off.terms[i]+ncol(object$smooth[[i]]$S[[1]])-1)] <- 
                         #             object$smooth[[i]]$p.ident[2:length(object$smooth[[i]]$p.ident)]
                       # else 
                            p.ident[off.terms[i]:(off.terms[i]+ncol(object$smooth[[i]]$S[[1]])-1)] <- 
                                                                       object$smooth[[i]]$p.ident
                   }
   ## getting the list of penalty matrices in terms of the full model vector of coefficients...
   S <- list()
   j <- 1
   if (n.terms>0) for(i in 1:n.terms){
                      for (kk in 1:length(object$smooth[[i]]$S)){
                            S[[j]] <- matrix(0,q,q) # initialize penalty matrix
                            S[[j]][off.terms[i]:(off.terms[i]+ncol(object$smooth[[i]]$S[[kk]])-1),
                   off.terms[i]:(off.terms[i]+ncol(object$smooth[[i]]$S[[kk]])-1)] <- object$smooth[[i]]$S[[kk]]
                            j <- j+1       
                      }
                   }
   object$S <- S 
   object$p.ident <- p.ident
   object
} ## end penalty_pident





#############################################################
## Function to fit SCAM based on Full Newton method        ##     
#############################################################

scam.fit <- function(G,sp, gamma=1, start=NULL, etastart=NULL, mustart=NULL, env=env, 
              null.coef=rep(0,ncol(G$X)), control=scam.control()) 
##  maxit=200, devtol.fit=1e-8, steptol.fit=1e-8, trace=FALSE, print.warn=FALSE
   ## G - list of items from gam(...,fit=FALSE) needed to fit a scam
   ## sp- vector of smoothing parameters
   ## not.exp - if TRUE then notExp() function will be used in place of exp in positivity ensuring beta parameters re-parameterization
   ## null.coef - coefficients for a null model, in order to be able to check for immediate divergence. null.coef give some sort of upper bound on deviance. This allows immediate divergence problems to be controlled (from mgcv). here some of coeff need to be exponentiated. 
   ## control - control list; it includes:
      ##* maxit - a positive scalar which gives the maximum number of iterations for Newton's method
      ##* devtol.fit - a scalar giving the tolerance at which the relative penalized deviance is considered to be close enougth to 0 to terminate the algorithm
      ##* steptol.fit - a scalar giving the tolerance at which the scaled distance between two successive iterates is considered close enough to zero to terminate the algorithm 
      ##* trace turns on or off some de-bugging information.
      ##* print.warn =FALSE - when set to 'FALSE' turns off printing warning messages for step halving under non-finite exponentiated coefficients,  non-finite deviance and/or if mu or eta are out of bounds. 
{ y <- G$y;  X <- G$X;  S <- G$S; not.exp <- G$not.exp;
  AR1.rho <- G$AR1.rho
  attr(X,"dimnames") <- NULL
  q0 <- G$q0; q.f <- G$q.f
  p.ident <- G$p.ident; n.terms <- G$n.terms
  family <- G$family; intercept <- G$intercept; offset <- G$offset;
  weights <- G$weights;  
  n <- nobs <- NROW(y)
  q <- ncol(X)
  dg <- fix.family.link(family)
  dv <- fix.family.var(family)
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
  
   if (AR1.rho!=0) {
        ld <- 1/sqrt(1-AR1.rho^2) ## leading diagonal of root inverse correlation
        sd <- -AR1.rho*ld         ## sub diagonal
        row <- c(1,rep(1:nobs,rep(2,nobs))[-c(1,2*nobs)])
        weight.r <- c(1,rep(c(sd,ld),nobs-1))
        end <- c(1,1:(nobs-1)*2+1) 
        if (!is.null(G$mf$"(AR.start)")) { ## need to correct the start of new AR sections...
                 ii <- which(G$mf$"(AR.start)"==TRUE)
                 if (length(ii)>0) {
                    if (ii[1]==1) ii <- ii[-1] ## first observation does not need any correction
                    weight.r[ii*2-2] <- 0 ## zero sub diagonal
                    weight.r[ii*2-1] <- 1 ## set leading diagonal to 1
                 }
        }
        ## apply transform...
        X <- rwMatrix(end,row,weight.r,X)
        y <- rwMatrix(end,row,weight.r,y)    
      }   



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

   ## Added code...
   if (family$family=="gaussian"&&family$link=="identity") strictly.additive <- TRUE else
      strictly.additive <- FALSE

   ## end of added code

  if (EMPTY) {
      eta <- rep.int(0, nobs) + offset
      if (!valideta(eta)) 
            stop("Invalid linear predictor values in empty model")
      mu <- linkinv(eta)

      if (AR1.rho!=0) {
        mu <- rwMatrix(end,row,weight.r,mu)  
      }   

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
      UBRE <- alpha/nobs + 2 * gamma* trA*scale/n - scale 
      scale.est <- alpha/(nobs - trA)
      aic.model <- aic(y, n, mu, weights, dev) +  2 * trA
  } ### end if (EMPTY)
  else {
      eta <- if (!is.null(etastart)) 
          etastart
            else family$linkfun(mustart)
      mu <- as.numeric(linkinv(eta))
    # if (!(validmu(mu) && valideta(eta))) 
    #        stop("Can't find valid starting values: please specify some")
      S.t <- matrix(0,q,q) # define the total sum of the penalties times sp
      n.pen <- length(S) # number of penalties 
      if (length(sp)!=n.pen) stop (paste("length of sp has to be equal to", n.pen))
      if (n.pen>0) for (j in 1:n.pen) S.t <- S.t + sp[j]*S[[j]]
      # get sqrt of total penalty matrix...
      er <- eigen(S.t,symmetric=TRUE);  er$values[er$values<0] <- 0
      rS <- crossprod(sqrt(sqrt(er$values))*t(er$vectors))
     # ii <- p.ident==1 ## set to TRUE/FALSE
     ## count <- sum(p.ident) 
     ## iv <- array(0, dim=c(count,1)) # define an index vector for the coeff-s to be exponentiated
      iv <- (1:q)[p.ident] # define an index vector for the coeff-s to be exponentiated
      ##############################################
      ## Initialization of parameters start here... 
      beta0 <- get("start",envir=env)
      dbeta0 <- get("dbeta.start",envir=env)
      sp.old <- get("sp.last",envir=env)
      if (length(beta0)==0) {
          # list argument to pcls for initializing model coefficients
          M <- list(X=X,p=rep(0.1,q),C=matrix(0,0,0),sp=sp,y=eta-offset,w=y*0+1) 
          M$Ain <- matrix(0,q,q); diag(M$Ain) <- rep(1,q);
          M$bin <- rep(-1e+12,q); M$bin[iv] <- 1e-12
          M$off <- rep(0,n.pen); M$S <- list()
          if (n.pen>0) for (j in 1:n.pen) {M$S[[j]] <- matrix(0,q,q); M$S[[j]] <- S[[j]]}
          beta.t <- pcls(M)      # initialize model coefficients (re-parameterized beta)
          beta <- beta.t         # initialize beta
          beta[iv] <- log(beta.t[iv]) # values of beta of the constrained terms
      }
      else {
          beta <- beta0
          beta.t <- beta    ## current beta tilde
          beta.t[iv] <- if (!not.exp) exp(beta[iv]) else notExp(beta[iv]) # values of re-para beta of the constrained term
      }
      ## Initialization of parameters finishes here 
      #################################################
      eta <- as.numeric(X%*%beta.t + as.numeric(offset))  # define initial linear predictor
      mu <- linkinv(eta)  # define initial fitted model
    ##  dev <- sum(dev.resids(y,mu,weights)) # define initial norm/deviance
    ##  pdev <- dev + sum((rS%*%beta)^2) # define initial penalized deviance 
    ## old.pdev <- pdev       # initialize convergence control for penalized deviance

      ##################################################################
      ## added code here made on May 6, 2020 for scam version 1-2-6,
      ## following gam.fit3(mgcv_1.8-31))....
      ################################################################### 
        ## shrink towards null.coef if immediately invalid 
      ##  betaold <- null.coef
      ##  etaold <- null.eta <- as.numeric(X%*%null.coef + as.numeric(offset))
      ##  old.pdev <- sum(dev.resids(y, linkinv(null.eta), weights)) + sum((rS%*%null.coef)^2) 
      ##  null.eta <- as.numeric(X%*%null.coef + as.numeric(offset))    
        ii <- 0
        while (!(validmu(mu) && valideta(eta))) { ## shrink towards null.coef if immediately invalid
          ii <- ii + 1
          if (ii>20) stop("Can't find valid starting values: please specify some")
          beta <- beta * .9 + null.coef * .1
          beta.t <- beta    ## current beta tilde
          beta.t[iv] <- if (!not.exp) exp(beta[iv]) else notExp(beta[iv])
          eta <- as.numeric(X%*%beta.t + offset)   
          mu <- linkinv(eta)
        }
        betaold <- null.coef <- beta
        betaold.t <- beta   
        betaold.t[iv] <- if (!not.exp) exp(betaold[iv]) else notExp(betaold[iv])
        etaold <- as.numeric(X%*%betaold.t + offset) 
        old.pdev <- sum(dev.resids(y,linkinv(etaold),weights)) + sum((rS%*%betaold)^2)
      ##################################################################

      pdev.plot <- 0     # define initial pen dev for plotting it 
      E <- matrix(0,q,q)   # define diagonal matrix E- second term of the Hessian
      Cdiag <- rep(1,q); C1diag <- rep(0,q)
    ##  if (!not.exp) {
    ##      Cdiag[iv] <- C1diag[iv] <- beta.t[iv]
    ##  } else {
    ##      Cdiag[iv] <- DnotExp(beta[iv]); C1diag[iv] <- D2notExp(beta[iv])
    ##  }
    ##  tX1 <- Cdiag*t(X)
    ##  g.deriv <- 1/mu.eta(eta)        # diag(G)
    ##  w1 <- weights/(variance(mu)*g.deriv^2)    # diag(W1)
    ##  Dp.g <- - drop(tX1%*%(w1*g.deriv*(y-mu))) + S.t%*%beta
                     # the gradient vector of the penalized deviance
    ##  Dp.gnorm <- max(abs(Dp.g)) # set convergence on the max value of the Dp.g
     ## betaold <- beta
      boundary <- conv <- FALSE
      old.warn <- getOption("warn")
      if (!control$print.warn) curr.warn <- -1 else curr.warn <- old.warn

     ###################################################
     ## MAIN ITERATIONS START HERE...
     #cat("\nscam.fit iter start")
     for (iter in 1:control$maxit)  {
         #cat(".")
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
       
         if (!not.exp) {
             Cdiag[iv] <- C1diag[iv] <- beta.t[iv]
         } else {
             Cdiag[iv] <- DnotExp(beta[iv]); C1diag[iv] <- D2notExp(beta[iv])
           }

         tX1 <- Cdiag*t(X)
         g.deriv <- 1/mu.eta(eta)  # diag(G)
         w1 <- weights/(variance(mu)*g.deriv^2)    # diag(W1) - Fisher weights
         y.mu <- y - mu
         alpha <- 1+ y.mu*(dv$dvar(mu)/variance(mu)+dg$d2link(mu)/g.deriv) # alpha elements of W
         w <- w1*alpha          # diag(W) - Newton weights
         diag(E) <- drop((C1diag*t(X))%*%(w1*g.deriv*y.mu))
         abs.w <- abs(w)      # absolute values of the diag(W)
         I.minus <- rep(0,nobs)  # define diagonal elements of the matrix I^{-}
         z1 <- g.deriv*y.mu/alpha  # the first term of the pseudodata
         iin <- w < 0;  I.minus[iin] <- 1;z1[iin] <- -z1[iin]
         wX11 <- rbind(sqrt(abs.w)[1:nobs]*t(tX1),rS)
           illcond <- FALSE
           Q <- qr(wX11,LAPACK=TRUE) 
           R <- qr.R(Q)
           rp <- 1:ncol(R)
           rp[Q$pivot] <- rp ## reverse pivot X=Q%*%R[,rp]
           if (Rrank(R)==ncol(R)) { ## no need to truncate, can just use QR
             R.inv <- backsolve(R,diag(ncol(R)))[rp,] ## inverse of unpivoted R
             tR.inv <- t(R.inv)
           } else { ## need SVD step
             R <- R[,rp] ## unpivoted R
             svd.r <- svd(R)
             d.inv <- rep(0,q)  # initial vector of inverse singular values
             good1 <- svd.r$d >= max(svd.r$d)*.Machine$double.eps^.5
             d.inv[good1] <- 1/svd.r$d[good1]
             if (sum(!good1)>0) illcond <- TRUE
             R <- svd.r$d*t(svd.r$v)
             Q <- qr.qy(Q,rbind(svd.r$u,matrix(0,nobs,q))) ## this is inefficient don't need to extract Q really
             tR.inv <- d.inv*t(svd.r$v)    # inverse of transpose of R
             R.inv <- t(tR.inv)
           }
           
          QtQRER <- tR.inv%*%(diag(E)*R.inv)
          if (sum(I.minus)>0) {
             if (is.qr(Q)) { 
            # QtQRER <- QtQRER + 2*tcrossprod(qr.qty(Q,diag(nrow(wX11))[,(1:nobs)[as.logical(I.minus)]]))
             QtQRER <- QtQRER + 2*crossprod(I.minus*qr.Q(Q)[1:nobs,])  
             } else {
             QtQRER <- QtQRER + 2*crossprod(I.minus*Q[1:nobs,])
             }
         }
         ei <- eigen(QtQRER,symmetric=TRUE)    # eigen-decomposition of QtQRER
         d <- ei$values        # vector of eigenvalues
         ok1 <- sum(d>1) > 0 # checking positive semi-definiteness 
         if (ok1 == TRUE) {# Fisher step in case of not +ve semi-definiteness of penalized loglikelihood
                         # set alpha =1
             eta.t <- drop(t(beta)%*%tX1)     # eta tilde for pseudodata 
             wX11 <- rbind(sqrt(w1)[1:nobs]*t(tX1),rS)   
             z<-g.deriv*y.mu+eta.t      # pseudodata
             wz<-w1^.5*z               # weighted pseudodata
             wz.aug<-c(wz,rep(0,nrow(rS)))   # augmented pseudodata
             Q <- qr(wX11,LAPACK=TRUE) 
             R <- qr.R(Q)
             rp <- 1:ncol(R)
             rp[Q$pivot] <- rp ## reverse pivot X=Q%*%R[,rp]
             if (Rrank(R)==ncol(R)) { ## no need to truncate, can just use QR
               beta <- backsolve(R,qr.qty(Q,wz.aug)[1:q])[rp]
             } else { ## need SVD step
               R <- R[,rp] ## unpivoted R
               s1 <- svd(R)
               d.inv1 <- rep(0,q)
               good1 <- s1$d >= max(s1$d)*.Machine$double.eps^.5
               d.inv1[good1] <- 1/s1$d[good1]
               beta <- s1$v%*%((d.inv1*t(s1$u))%*%qr.qty(Q,wz.aug)[1:q])
             }
         }  ### end of if (ok1) - Fisher step 
         else  {       ##  full Newton step
             Id.inv.r<-1/(1-d)^.5   # vector of inverse values of (1-sv)^.5
             iin <- (1-d) < .Machine$double.eps
             Id.inv.r[iin] <- 0
             eidrop <- t(Id.inv.r*t(ei$vectors))
             wz1<-sqrt(abs.w)*z1  # the first term of the weighted pseudodata
             if (is.qr(Q)) {
               beta <- R.inv%*%(eidrop%*%(t(eidrop)%*%qr.qty(Q,c(wz1,rep(0,nrow(rS))))[1:nrow(eidrop)]))
             } else {
               beta <- R.inv%*%(eidrop%*%(t(eidrop)%*%(t(Q[1:nobs,])%*%wz1)[1:nrow(eidrop)]))
             }
             beta <- betaold + 
                     drop(beta - R.inv%*%(eidrop%*%(t(eidrop)%*%(tR.inv%*%(S.t%*%betaold)))))
         }  ###  end of if (!ok1) - Newton step


      ##   delta <- beta-c(betaold)         # trial step
      ##   step <- 1                      # initial trial step length
      ##   beta <- c(betaold)+step*delta    # current parameter estimates
         beta.t <- beta                 # current reparameterized beta
         beta.t[iv] <- if (!not.exp) exp(beta[iv]) else notExp(beta[iv])  # values of re-para beta of the shape constrained term
       ##  eta <- drop(as.numeric(X%*%beta.t + offset))  # linear predictor
       ##  mu <- linkinv(eta)          # fitted values
       ##  dev <- sum(dev.resids(y,mu,weights)) # deviance of the working model
       ##  pdev <- dev + sum((rS%*%beta)^2) # deviance + penalty of the working model

         ##################################################################
         ## added code here made on May 6, 2020 for scam version 1-2-6,
         ## following gam.fit3(mgcv_1.8-31))....
         ###################################################################   
         if (any(!is.finite(beta))) {
                conv <- FALSE
                warning(gettextf("Non-finite coefficients at iteration %d", 
                  iter))
                break
         }  

         options(warn = curr.warn) ## turn of/on printing warning messages under non-finite deviance resulted in non-finite exponentiated coefficients...      
         if (any(!is.finite(beta.t))) {
                conv <- FALSE
                warning(gettextf("Non-finite exponentiated coefficients at iteration %d", 
                  iter))
             #   break ## takes out of the main iteration loop 'iter in 1:maxit'
         }     

         eta <- drop(as.numeric(X%*%beta.t + offset))  # linear predictor
         mu <- linkinv(eta)          # fitted values
         dev <- sum(dev.resids(y,mu,weights)) # deviance of the working model
         penalty <- sum((rS%*%beta)^2)

         if (control$trace) 
             message(gettextf("Deviance = %s Iterations - %d", dev, iter, domain = "R-scam"))
          boundary <- FALSE
         
         ## step halve under non-finite deviance...
         if (!is.finite(dev)) {
                if (is.null(betaold)) {
                  if (is.null(null.coef)) 
                  stop("no valid set of coefficients has been found:please supply starting values", 
                    call. = FALSE)
                  ## Try to find feasible coefficients from the null.coef and null.eta
                  betaold <- null.coef                
                }
                warning("Step size truncated due to divergence", 
                  call. = FALSE)
                ii <- 0
             ##   step <- 1
                while (!is.finite(dev)) {
                  if (ii > control$maxit) 
                    stop("inner loop 1; can't correct step size")
                  ii <- ii + 1
               ## beta <- (beta + betaold)/2 
               ## eta <- (eta + etaold)/2               
               ## mu <- linkinv(eta)
               ## step <- step/2      ## decrease step length 
               ## beta <- c(betaold)+step*delta   ## update current parameter estimates
                  beta <- (beta + c(betaold))/2 
                  beta.t <- beta       ## update current re-para beta
                  beta.t[iv] <- if (!not.exp) exp(beta[iv]) else notExp(beta[iv])   
                  eta <- as.numeric(X%*%beta.t + offset)   ## linear predictor  
                  mu <- linkinv(eta)   ## fitted values
                  dev <- sum(dev.resids(y, mu, weights))
                }
                boundary <- TRUE 
                penalty <- sum((rS%*%beta)^2) ## t(beta)%*%S.t%*%beta # reset penalty too
                if (control$trace) 
                  cat("Step halved: new deviance =", dev, "\n")
            } ## end of infinite deviance correction

          ## now step halve if mu or eta are out of bounds... 
          if (!(valideta(eta) && validmu(mu))) {
                warning("Step size truncated: out of bounds", 
                  call. = FALSE)
                ii <- 0
               ## step <- 1
                while (!(valideta(eta) && validmu(mu))) {
                  if (ii > control$maxit) 
                    stop("inner loop 2; can't correct step size")
                  ii <- ii + 1
               ## beta <- (beta + betaold)/2 
               ## eta <- (eta + etaold)/2               
               ## mu <- linkinv(eta)
               ## step <- step/2      ## decrease step length 
               ## beta <- c(betaold)+step*delta   ## update current parameter estimates
                  beta <- (beta + c(betaold))/2 ## update current parameter estimates
                  beta.t <- beta       ## update current re-para beta
                  beta.t[iv] <- if (!not.exp) exp(beta[iv]) else notExp(beta[iv])   
                  eta <- as.numeric(X%*%beta.t + offset)   ## linear predictor  
                  mu <- linkinv(eta)   ## fitted values
                }
                boundary <- TRUE
                penalty <- sum((rS%*%beta)^2) ## t(beta)%*%S.t%*%beta
                dev <- sum(dev.resids(y, mu, weights))
                if (control$trace) 
                  cat("Step halved: new deviance =", dev, "\n")
          } ## end of invalid mu/eta handling

         options(warn = old.warn) ## return to the default option for printing warning messages (which is '0')

         ## now check for divergence of penalized deviance....
         pdev <- dev +  penalty  ## the penalized deviance   

         if (control$trace) 
                message(gettextf("penalized deviance = %s", pdev, domain = "R-scam"))

         ######################################################
         ## end of added coding
         ######################################################        

         div.thresh <- 10*(.1 +abs(old.pdev))*.Machine$double.eps^.5

         ## `step reduction' approach starts here...
       ##  ii <- 1  
       ##  while ( (pdev-old.pdev) > div.thresh) { # 'step reduction' approach
       ##      if (ii > maxHalf.fit) 
       ##          break ## stop ("step reduction failed")
       ##      ii <- ii+1
       ##      step <- step/2         # decrease step length 
       ##      beta <- c(old.beta)+step*delta   # update current parameter estimates
       ##      beta.t <- beta                # update current re-para beta
       ##      beta.t[iv] <- if (!not.exp) exp(beta[iv]) else notExp(beta[iv])   
       ##      eta <- as.numeric(X%*%beta.t + offset)    # linear predictor  
       ##      mu <- linkinv(eta)         # fitted values
       ##      dev <- sum(dev.resids(y,mu,weights)) # update deviance of the working model
       ##      pdev <- dev+sum((rS%*%beta)^2) # update pen deviance of the working model
       ##  }
         ## `step reduction' finishes here 
     
           ## ... threshold for judging divergence --- too tight and near
           ## perfect convergence can cause a failure here

          if (pdev-old.pdev>div.thresh) { ## solution diverging
             ii <- 0 ## step halving counter
           ##  step <- 1
           ##  if (iter==1) { ## immediate divergence, need to shrink towards zero 
           ##     betaold <- null.coef
           ##  }
             while (pdev -old.pdev > div.thresh){ ## step halve until pdev <= old.pdev
                if (ii > 100) 
                   stop("inner loop 3; can't correct step size")
                ii <- ii + 1
              ## beta <- (beta + betaold)/2 
              ## eta <- (eta + etaold)/2               
              ## mu <- linkinv(eta)
              ##  step <- step/2      ## decrease step length 
              ## beta <- c(betaold)+step*delta   ## update current parameter estimates
                beta <- (beta + c(betaold))/2
                beta.t <- beta       ## update current re-para beta
                beta.t[iv] <- if (!not.exp) exp(beta[iv]) else notExp(beta[iv])   
                eta <- as.numeric(X%*%beta.t + offset)   ## linear predictor  
                mu <- linkinv(eta)   ## fitted values
                dev <- sum(dev.resids(y, mu, weights)) ## update deviance of the working model
                pdev <- dev + sum((rS%*%beta)^2) ## update the penalized deviance
                if (control$trace) 
                  message(gettextf("Step halved: new penalized deviance = %g", pdev, "\n"))
              }
          } ## end of pdev divergence

         Dp.g <- -drop(tX1%*%(w1*g.deriv*(y-mu)))+S.t%*%beta # the gradient vector of the penalized deviance
         Dp.gnorm <- max(abs(Dp.g)) 
         pdev.plot[iter] <- pdev      # store penilized deviance of the working model for plotting
          
      ##   if (strictly.additive) { conv <- TRUE; break;}

         ## checking convergence adequately...
         if (abs(pdev - old.pdev)/(.1 + abs(pdev)) < control$devtol.fit) {
              if (Dp.gnorm > control$devtol.fit*max(abs(beta+c(betaold)))/2) {
           ## if (max(abs(beta - c(betaold))) > control$steptol.fit*max(abs(beta + c(betaold)))/2) {
                old.pdev <- pdev ## not converged quite enough
                betaold <- beta                   
              }
              else { ## converged
                 conv <- TRUE
                 beta <- beta              
                 break
              }
         }
         else { ## not converged
             old.pdev <- pdev
             betaold <- beta              
         }
     } ## main iteration procedure is completed here.
     #######################################################
     ## at this stage the model has been fully estimated


     ## now define matrices at their converged values from the full Newton method...
  
   #  beta.t <- beta   ## estimates of re-para beta
   #  beta.t[iv] <- if (!not.exp) exp(beta[iv]) else notExp(beta[iv])   
   #  eta <- as.numeric(X%*%beta.t + offset)      ## linear predictor  
   #  mu <- linkinv(eta)   ## fitted values
   #  dev <- sum(dev.resids(y,mu,weights))

     if (!not.exp) {
             Cdiag[iv] <- C1diag[iv] <- beta.t[iv]
     } else {
             Cdiag[iv] <- DnotExp(beta[iv]); C1diag[iv] <- D2notExp(beta[iv])
       }

     X1 <- t(Cdiag*t(X)) 
     g.deriv <- 1/ mu.eta(eta)        # diag(G)
     w1 <- weights/(variance(mu)*g.deriv^2)    # diag(W1)
     alpha <- 1+(y-mu)*(dv$dvar(mu)/variance(mu)+dg$d2link(mu)/g.deriv) # alpha elements of W
     w <- w1*alpha          # diag(W)

     diag(E) <- drop((C1diag*t(X))%*%(w1*g.deriv*(y-mu))) # diagonal elements of E
     abs.w <- abs(w)      # absolute values of the diag(W)
     I.minus <- rep(0,nobs)  # define diagonal elements of the matrix I^{-}
     I.minus[w<0] <- 1
     wX1 <- sqrt(abs.w)[1:nobs]*X1 ## wX1 actually used later
     wX11 <- rbind(wX1,rS) # augmented model matrix 
     ## Faster version only does SVD when needed (and then only on QR factor)
        illcond <- FALSE
        Q <- qr(wX11,LAPACK=TRUE) 
        R <- qr.R(Q)
        rp <- 1:ncol(R)
        rp[Q$pivot] <- rp ## reverse pivot X=Q%*%R[,rp]
 
        R.out <- R[,rp]  ## unpivoted R, needed for summary function
        rank <- Rrank(R)

        if (rank==ncol(R)) { ## no need to truncate, can just use QR
           R.inv <- backsolve(R,diag(ncol(R)))[rp,] ## inverse of unpivoted R
           tR.inv <- t(R.inv)
        } else { ## need SVD step
           R <- R[,rp] ## unpivoted R
           svd.r <- svd(R)
           d.inv <- rep(0,q)  # initial vector of inverse singular values
           good1 <- svd.r$d >= max(svd.r$d)*.Machine$double.eps^.5
           d.inv[good1] <- 1/svd.r$d[good1]
           if (sum(!good1)>0) illcond <- TRUE
           R <- svd.r$d*t(svd.r$v)
           Q <- qr.qy(Q,rbind(svd.r$u,matrix(0,nobs,q))) ## this is inefficient don't need to extract Q really
           tR.inv <- d.inv*t(svd.r$v)    # inverse of transpose of R
           R.inv <- t(tR.inv)
         }
       QtQRER <- tR.inv%*%(diag(E)*R.inv)
         if (sum(I.minus)>0) {
              if (is.qr(Q)) { 
              QtQRER <- QtQRER + 2*crossprod(I.minus*qr.Q(Q)[1:nobs,])  
              } else {
              QtQRER <- QtQRER + 2*crossprod(I.minus*Q[1:nobs,])
              }
         }
      ei <- eigen(QtQRER,symmetric=TRUE)    # eigen-decomposition of QtQRER
      d <- ei$values        # vector of eigenvalues
      ok1 <- sum(d>1)>0
      if (ok1) { ## Fisher step in case of not positive semi-definiteness of penalized loglikelihood
                 ## set alpha =1 ...
         wX1<-sqrt(w1)[1:nobs]*X1
         wX11<-rbind(wX1,rS)
         Q <- qr(wX11,LAPACK=TRUE) 
         R <- qr.R(Q)
         rp <- 1:ncol(R)
         rp[Q$pivot] <- rp ## reverse pivot X=Q%*%R[,rp]

         R.out <- R[,rp]  ## unpivoted R, needed for summary function
         rank <- Rrank(R)
         if (rank==ncol(R)) { ## no need to truncate, can just use QR
             P <- backsolve(R,diag(ncol(R)))[rp,]
             K <- qr.Q(Q)[1:nobs,]
         } else { ## need SVD step
             R <- R[,rp] ## unpivoted R
             s1 <- svd(R)
             d.inv1 <- rep(0,q)
             good1 <- s1$d >= max(s1$d)*.Machine$double.eps^.5
             d.inv1[good1] <- 1/s1$d[good1]
             P <- t(d.inv1*t(s1$v))
             K <- qr.qy(Q,rbind(s1$u,matrix(0,nobs,q)))[1:nobs,]
          }
     } ## end of if (ok1)
     else  {   ## full Newton step
         Id.inv.r<-1/(1-d)^.5   # vector of inverse values of (1-sv)^1/2
         ii <- (1-d) < .Machine$double.eps
         Id.inv.r[ii] <- 0
         eidrop <- t(Id.inv.r*t(ei$vectors))
         P <- R.inv%*%eidrop  ## ei$vectors%*%diag(Id.inv.r) # define matrix P
         if (is.qr(Q)) {
           K <- qr.qy(Q,rbind(eidrop,matrix(0,nobs,q)))[1:nobs,]
         } else {
           K <- Q[1:nobs,]%*%eidrop ## (ei$vectors%*%diag(Id.inv.r))  # define matrix K 
         }
     }  ## end of if (!ok1)
     # end of calculation of the matrices at their converged values...
     

     Dp.g <- -t(X1)%*%(w1*g.deriv*(y-mu))+S.t%*%beta # the gradient vector of the penalized deviance
     Dp.gnorm<-max(abs(Dp.g)) 

     # calculating tr(A)...
     I.plus <- rep(1,nobs)   # define diagonal elements of the matrix I^{+}
     I.plus[w<0] <- -1
     L <- c(1/alpha)    # define diagonal elements of L=diag(1/alpha)
     ## NOTE PKt is O(np^2) and not needed --- can compute trA as side effect of gradiant
     KtILQ1R <- crossprod(L*I.plus*K,wX1) ## t(L*I.plus*K)%*%wX1 
     edf <- rowSums(P*t(KtILQ1R))
     trA <- sum(edf)

     ## below is needed for later calculations of derivative of trA
     wXC1 <- sqrt(abs.w)[1:nobs]*t(C1diag*t(X)) 
     KtILQ1R <-  if (!not.exp) KtILQ1R else crossprod(L*I.plus*K,wXC1)  ##t(L*I.plus*K)%*%wXC1
     KtIQ1R <-  if (!not.exp) crossprod(I.plus*K,wX1) 
                  else crossprod(I.plus*K,wXC1)
     C2diag <- rep(0,q)
     C2diag[iv] <- if (!not.exp)  C1diag[iv] else D3notExp(beta[iv])
     XC2 <- t(C2diag*t(X))
     XC1 <- t(C1diag*t(X))
     
     ############################
     ## some return values....
     dlink.mu  <- 1/mu.eta(eta); Var<- variance(mu)
     link <- family$linkfun(mu); d2link.mu <- dg$d2link(mu)
     dvar.mu <- dv$dvar(mu); d2var.mu <- dv$d2var(mu)
     d3link.mu <- dg$d3link(mu)
     z <- g.deriv*(y-mu)+X1%*%beta
     #############################
     #############################
     scale.est <- dev/(nobs-trA)  #  scale estimate
     ## re-transforming 'mu' back on original scale and hence residuals, in case of correlated errors...
     residuals <- rep.int(NA, nobs)
     residuals <- (y-mu)*g.deriv
   #  if (AR1.rho==0)  residuals <- (y-mu)*g.deriv
   #  else {
   #         eta <- as.numeric(G$X%*%beta.t + offset) # linear predictor on original scale
   #         mu <- linkinv(eta)    # fitted values on original scale 
   #         g.deriv <- 1/ mu.eta(eta)        # diag(G)
   #         residuals <- (G$y-mu)*g.deriv
   #  }
     ################################
     ## calculation of the derivatives of beta by the Implicit Function Theorem starts here
     dbeta.rho <- matrix(0,q,n.pen) # define matrix of the parameters derivatives
     if (n.pen>0) for (j in 1:n.pen) {
                     dbeta.rho[,j] <- - sp[j]*P%*%(t(P)%*%(S[[j]]%*%beta)) # derivative of beta wrt rho[j]
                  }
     # end of calculating the parameters derivatives
     aic.model <- aic(y, n, mu, weights, dev) +  2 * sum(edf)
     if (AR1.rho!=0) { ## correct aic for AR1 transform
        df <- 1 ## if (getARs) sum(b$model$"(AR.start)") else 1
        aic.model <- aic.model - 2*(n-df)*log(ld) 
     }
     assign("start",beta,envir=env)
     assign("dbeta.start",dbeta.rho,envir=env)
     assign("sp.last",sp,envir=env)
  } ### end if (!EMPTY) 

 list(L=L,C1diag=C1diag,E=E,iter=iter, old.beta=betaold, gcv=dev*nobs/(nobs-gamma *trA)^2, sp=sp,
      mu=mu,X=G$X,y=drop(G$y), X1=X1,beta=beta,beta.t=beta.t,iv=iv,S=S,S.t=S.t,rS=rS,
      P=P,K=K, C2diag=C2diag, KtILQ1R= KtILQ1R, KtIQ1R=KtIQ1R, ## XC1=XC1, XC2=XC2,
      dlink.mu=dlink.mu,Var=Var, abs.w=drop(abs.w),
      link=link,w=as.numeric(w),w1=drop(w1),d2link.mu=d2link.mu,wX1=wX1,I.plus=I.plus,
      dvar.mu=dvar.mu,d2var.mu=d2var.mu,deviance=dev,scale.est=scale.est,
      ok1=ok1,alpha=as.numeric(alpha),d3link.mu=d3link.mu,eta=eta,iter=iter,
      Dp.gnorm=Dp.gnorm, Dp.g=Dp.g,d=d, conv=conv, illcond=illcond,R=R.out, edf=edf,trA=trA,
      residuals=residuals,z=z,dbeta.rho=dbeta.rho, aic=aic.model,rank=rank) ##step=step, AR1.rho=AR1.rho,AR.start=mf$"(AR.start)")
} ## end of scam.fit




#######################################################################
## function to get null deviance and covariance matrices after fit   ##
#######################################################################


scam.fit.post<- function(G, object) ##,sig2,offset,intercept, weights,scale.known, not.exp)
{  ## Function to compute null deviance and covariance matrices after a scam fit.
   ## covariance matrix should use expected Hessian, so re-computation of factors 
   ## is required.  
   ## object - object from estimate.scam()
   y <- G$y; X <- G$X;
   sig2 <- G$sig2; offset <- G$offset; intercept <- G$intercept; 
   weights <- G$weights; scale.known <- G$scale.known; not.exp <- G$not.exp
   n <- nobs <- NROW(y) # number of observations
   q <- ncol(X)
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
        X <- rwMatrix(end,row,weight.r,X)
        y <- rwMatrix(end,row,weight.r,y)    
      }   
  
   linkinv <- object$family$linkinv
   dev.resids <- object$family$dev.resids
   dg <- fix.family.link(object$family)
   dv <- fix.family.var(object$family)

   eta <- as.numeric(X%*%object$beta.t + offset) # linear predictor
   mu <- linkinv(eta)          # fitted values
   dev <- sum(dev.resids(y,mu,weights)) # deviance of the final model
   
   wtdmu <- if (intercept) sum(weights * G$y)/sum(weights) 
              else linkinv(offset)

   null.dev <- sum(dev.resids(G$y, wtdmu, weights))
 
   n.ok <- nobs - sum(weights == 0)
   nulldf <- n.ok - as.integer(intercept)
   
   ## define matrices at their converged values from the full Newton method...
   Cdiag <- rep(1,q); C1diag <- rep(0,q)
   iv <- object$iv
   if (!not.exp) {
             Cdiag[iv] <- C1diag[iv] <- object$beta.t[iv]
     } else {
             Cdiag[iv] <- DnotExp(object$beta[iv]); C1diag[iv] <- D2notExp(object$beta[iv])
      }
    X1 <- t(Cdiag*t(X)) 
    g.deriv <- 1/object$family$mu.eta(eta)        # diag(G)
    w1 <- weights/(object$family$variance(mu)*g.deriv^2)    # diag(W1)
    alpha <- 1+(y-mu)*(dv$dvar(mu)/object$family$variance(mu)+dg$d2link(mu)/g.deriv) # alpha elements of W
    w <- w1*alpha         # diag(W)
    E <- matrix(0,q,q)
    diag(E) <- drop((C1diag*t(X))%*%(w1*g.deriv*(y-mu))) # diagonal elements of E
    abs.w <- abs(w)      # absolute values of the diag(W)
    I.minus <- rep(0,nobs)  # define diagonal elements of the matrix I^{-}
    I.minus[w<0] <- 1
    wX1 <- sqrt(abs.w)[1:nobs]*X1 
    wX11 <- rbind(wX1,object$rS) # augmented model matrix 
     ## Faster version only does SVD when needed (and then only on QR factor)
        illcond <- FALSE
        Q <- qr(wX11,LAPACK=TRUE) 
        R <- qr.R(Q)
        rp <- 1:ncol(R)
        rp[Q$pivot] <- rp ## reverse pivot X=Q%*%R[,rp]
 
        R.out <- R[,rp]  ## unpivoted R, needed for summary function
        rank <- Rrank(R)
        if (rank==ncol(R)) { ## no need to truncate, can just use QR
           R.inv <- backsolve(R,diag(ncol(R)))[rp,] ## inverse of unpivoted R
           tR.inv <- t(R.inv)
        } else { ## need SVD step
           R <- R[,rp] ## unpivoted R
           svd.r <- svd(R)
           d.inv <- rep(0,q)  # initial vector of inverse singular values
           good <- svd.r$d >= max(svd.r$d)*.Machine$double.eps^.5
           d.inv[good] <- 1/svd.r$d[good]
           if (sum(!good)>0) illcond <- TRUE
           R <- svd.r$d*t(svd.r$v)
           Q <- qr.qy(Q,rbind(svd.r$u,matrix(0,nobs,q))) ## this is inefficient don't need to extract Q really
           tR.inv <- d.inv*t(svd.r$v)    # inverse of transpose of R
           R.inv <- t(tR.inv)
         }
       QtQRER <- tR.inv%*%(diag(E)*R.inv)
         if (sum(I.minus)>0) {
              if (is.qr(Q)) { 
              QtQRER <- QtQRER + 2*crossprod(I.minus*qr.Q(Q)[1:nobs,])  
              } else {
              QtQRER <- QtQRER + 2*crossprod(I.minus*Q[1:nobs,])
              }
         }
      ei <- eigen(QtQRER,symmetric=TRUE)    # eigen-decomposition of QtQRER
      d <- ei$values        # vector of eigenvalues
      ok1 <- sum(d>1)>0
      if (ok1) { ## Fisher step in case of not positive semi-definiteness of penalized loglikelihood
                 ## set alpha =1 ...
         wX1<-sqrt(w1)[1:nobs]*X1
         wX11<-rbind(wX1,object$rS)
         Q <- qr(wX11,LAPACK=TRUE) 
         R <- qr.R(Q)
         rp <- 1:ncol(R)
         rp[Q$pivot] <- rp ## reverse pivot X=Q%*%R[,rp]

         R.out <- R[,rp]  ## unpivoted R, needed for summary function
         rank <- Rrank(R)
         if (rank==ncol(R)) { ## no need to truncate, can just use QR
             P <- backsolve(R,diag(ncol(R)))[rp,]
             K <- qr.Q(Q)[1:nobs,]
         } else { ## need SVD step
             R <- R[,rp] ## unpivoted R
             s1 <- svd(R)
             d.inv1 <- rep(0,q)
             good1 <- s1$d >= max(s1$d)*.Machine$double.eps^.5
             d.inv1[good1] <- 1/s1$d[good1]
             P <- t(d.inv1*t(s1$v))
             K <- qr.qy(Q,rbind(s1$u,matrix(0,nobs,q)))[1:nobs,]
          }
     } ## end of if (ok1)
     else  {   ## full Newton step
         Id.inv.r<-1/(1-d)^.5   # vector of inverse values of (1-sv)^1/2
         ii <- (1-d) < .Machine$double.eps
         Id.inv.r[ii] <- 0
         eidrop <- t(Id.inv.r*t(ei$vectors))
         P <- R.inv%*%eidrop  ## ei$vectors%*%diag(Id.inv.r) # define matrix P
         if (is.qr(Q)) {
           K <- qr.qy(Q,rbind(eidrop,matrix(0,nobs,q)))[1:nobs,]
         } else {
           K <- Q[1:nobs,]%*%eidrop ## (ei$vectors%*%diag(Id.inv.r))  # define matrix K 
         }
     }  ## end of if (!ok1)
     # end of calculation of the matrices at their converged values...

    ## calculating edf and trA...
    I.plus <- rep(1,nobs)   # define diagonal elements of the matrix I^{+}
    I.plus[w<0] <- -1
    L <- c(1/alpha)
    KtILQ1R <- crossprod(L*I.plus*K,wX1) ## t(object$L*object$I.plus*K)%*%object$wX1
    F <- P%*%(KtILQ1R)
    edf <- diag(F) ## effective degrees of freedom
    edf1 <- 2*edf - rowSums(t(F)*F) ## alternative
    trA <- sum(edf)
    
   ## calculating the approximate covariance matrices 
   ## (dealing with the expected Hessian of the log likelihood) ...
   ## get the inverse of the expected Hessian...
   if (!scale.known) sig2 <-  dev/(nobs-trA)  # scale estimate
   Vb <- tcrossprod(P) * sig2 
          ## P%*%t(P)*sig2 # Bayesian posterior covariance matrix for the parameters 
   Ve <- crossprod(K%*%t(P)) *sig2
        #PKt%*%t(PKt)*sig2 # covariance matrix of the parameter estimators 
   ## Delta method to get covariance matrix for the reparametrized parameters...
   df.p <- rep(1,q)
   df.p[object$iv] <- object$beta.t[object$iv]
   Vb.t <- t(df.p*t(df.p*Vb))
   Ve.t <- t(df.p*t(df.p*Ve))
   
   eta <- as.numeric(G$X%*%object$beta.t + offset)
   mu <- linkinv(eta)     # fitted values
   residuals <- rep.int(NA, nobs)
   g.deriv <- 1/object$family$mu.eta(eta) # diag(G)
   residuals <- (G$y-mu)*g.deriv # the working residuals for the fitted model
      
   aic.model <- object$family$aic(y, n, mu, weights, dev) +  2 * sum(edf)
   if (G$AR1.rho!=0) { ## correct aic for AR1 transform
        df <- 1 ## if (getARs) sum(G$AR.start) else 1
        aic.model <- aic.model - 2*(n-df)*log(1/sqrt(1-G$AR1.rho^2)) 
   }
    
   list (null.dev=null.dev, df.null=nulldf,Vb=Vb,Vb.t=Vb.t,Ve=Ve,Ve.t=Ve.t,rank=rank,
        sig2=sig2,edf=edf,edf1=edf1,trA=trA, deviance=dev,residuals=residuals, aic=aic.model, mu=mu, eta=eta)
} ## end of scam.fit.post



### the following three functions are for use in place of exp(beta)
### notExp() is similar to that in R package mgcv() of Simon N Wood
### in positivity ensuring beta parameters re-parameterization.... they have `better' 
### over/underflow characteristics, but is still continuous to second
### derivative. 
### DnotExp() calculates the first derivative
### D2notExp() gets the second derivative 

notExp <- function(x){
  f <- x
  ind <- x > 1
  f[ind] <- exp(1)*(x[ind]^2+1)/2
  ind <- (x <= 1)&(x > -1)
  f[ind] <- exp(x[ind])
  ind <- (x <= -1)
  f[ind] <-  exp(1)*(x[ind]^2+1)/2; f[ind]<-1/f[ind]
  f
}

DnotExp <- function(x) {
## first derivative of notExp()...
  f <- x
  ind <- x > 1
  f[ind] <- exp(1)*x[ind]
  ind <- (x <= 1)&(x > -1)
  f[ind] <- exp(x[ind])
  ind <- (x <= -1)
  f[ind] <-  -4*x[ind]/exp(1)/(x[ind]^2+1)^2
  f
}

D2notExp <- function(x) {
## second derivative of notExp()...
  f <- x
  ind <- x > 1
  f[ind] <- exp(1)
  ind <- (x <= 1)&(x > -1)
  f[ind] <- exp(x[ind])
  ind <- (x <= -1)
  f[ind] <-  (12*x[ind]^2-4)/exp(1)/(x[ind]^2+1)^3
  f
}


D3notExp <- function(x) {
## third derivative of notExp()...
  f <- x
  ind <- x > 1
  f[ind] <- 0
  ind <- (x <= 1)&(x > -1)
  f[ind] <- exp(x[ind])
  ind <- (x <= -1)
  f[ind] <-  48*x[ind]*(1-x[ind]^2)/exp(1)/(x[ind]^2+1)^4
  f
}




## checking derivatives...
#eps <- 1e-7
#x<-seq(-5,5,length.out=100)
#d1 <- (notExp(x+eps)-notExp(x))/eps
#range((DnotExp(x)-d1)/d1)
#d2 <- (DnotExp(x+eps)-DnotExp(x))/eps
#range((D2notExp(x)-d2)/d2)
#d3 <- (D2notExp(x+eps)-D2notExp(x))/eps
#range((D3notExp(x)-d3)/d3)


logLik.scam <- function (object,...)
{  # based on logLik.gam and logLik.glm 
    sc.p <- as.numeric(object$scale.estimated)
    p <- sum(object$edf) + sc.p
    val <- p - object$aic/2
    #if (fam %in% c("gaussian", "Gamma", "inverse.gaussian","Tweedie"))
    #    p <- p + 1
    np <- length(object$coefficients) + sc.p 
    if (p > np) p <- np 
    attr(val, "df") <- p
    class(val) <- "logLik"
    val
} ## logLik.scam


formula.scam <- function(x, ...)
# clone of formula.gam...
# formula.lm and formula.glm reconstruct the formula from x$terms, this is 
# problematic because of the way mgcv handles s() and te() terms 
{ x$formula
}


###############################################################################################
## below are functions from mgcv package, copied as they are not exported by 'namespace:mgcv' 
## Copyright (c) Simon N. Wood 2008-2019 simon.wood@r-project.org

gam.setup <- function(formula,pterms,
                     data=stop("No data supplied to gam.setup"),knots=NULL,sp=NULL,
                    min.sp=NULL,H=NULL,absorb.cons=TRUE,sparse.cons=0,select=FALSE,idLinksBases=TRUE,
                    scale.penalty=TRUE,paraPen=NULL,gamm.call=FALSE,drop.intercept=FALSE,
                    diagonal.penalty=FALSE,apply.by=TRUE,list.call=FALSE,modCon=0) 
## set up the model matrix, penalty matrices and auxilliary information about the smoothing bases
## needed for a gam fit.
## elements of returned object:
## * m - number of smooths
## * min.sp - minimum smoothing parameters
## * H supplied H matrix
## * pearson.extra, dev.extra, n.true --- entries to hold these quantities
## * pterms - terms object for parametric terms
## * intercept TRUE if intercept present
## * offset - the model offset
## * nsdf - number of strictly parameteric coefs
## * contrasts 
## * xlevels - records levels of factors
## * assign - indexes which parametric model matrix columns map to which term in pterms
## * smooth - list of smooths
## * S - penalties (non-zero block only)
## * off - first coef penalized by each element of S
## * cmX - col mean of X
## * P - maps parameters in fit constraint parameterization to those in prediction parameterization
## * X - model matrix
## * sp
## * rank
## * n.paraPen
## * L 
## * lsp0
## * y - response
## * C - constraint matrix - only if absorb.cons==FALSE
## * n - dim(y)
## * w - weights
## * term.names
## * nP
{ # split the formula if the object being passed is a formula, otherwise it's already split

  if (inherits(formula,"split.gam.formula")) split <- formula else
  if (inherits(formula,"formula")) split <- interpret.gam(formula) 
  else stop("First argument is no sort of formula!") 
  
  if (length(split$smooth.spec)==0) {
    if (split$pfok==0) stop("You've got no model....")
    m <- 0
  } else  m <- length(split$smooth.spec) # number of smooth terms
  
  G <- list(m=m,min.sp=min.sp,H=H,pearson.extra=0,
            dev.extra=0,n.true=-1,pterms=pterms) ## dev.extra gets added to deviance if REML/ML used in gam.fit3
  
  if (is.null(attr(data,"terms"))) # then data is not a model frame
  mf <- model.frame(split$pf,data,drop.unused.levels=FALSE) # must be false or can end up with wrong prediction matrix!
  else mf <- data # data is already a model frame

  G$intercept <-  attr(attr(mf,"terms"),"intercept")>0

  ## get any model offset. Complicated by possibility of offsets in multiple formulae...
  if (list.call) {
    offi <- attr(pterms,"offset")
    if (!is.null(offi)) {
      G$offset <- mf[[names(attr(pterms,"dataClasses"))[offi]]]
    }
  } else G$offset <- model.offset(mf)   # get any model offset including from offset argument
  
  if (!is.null(G$offset))  G$offset <- as.numeric(G$offset) 

  # construct strictly parametric model matrix.... 
  if (drop.intercept) attr(pterms,"intercept") <- 1 ## ensure there is an intercept to drop
  X <- model.matrix(pterms,mf)
  if (drop.intercept) { ## some extended families require intercept to be dropped 
    xat <- attributes(X);ind <- xat$assign>0 ## index of non intercept columns 
    X <- X[,ind,drop=FALSE] ## some extended families need to drop intercept
    xat$assign <- xat$assign[ind];xat$dimnames[[2]]<-xat$dimnames[[2]][ind];
    xat$dim[2] <- xat$dim[2]-1;attributes(X) <- xat
    G$intercept <- FALSE
  } 
  rownames(X) <- NULL ## save memory
  
  G$nsdf <- ncol(X)
  G$contrasts <- attr(X,"contrasts")
  G$xlevels <- .getXlevels(pterms,mf)
  G$assign <- attr(X,"assign") # used to tell which coeffs relate to which pterms

  ## now deal with any user supplied penalties on the parametric part of the model...
  PP <- parametricPenalty(pterms,G$assign,paraPen,sp)
  if (!is.null(PP)) { ## strip out supplied sps already used
    ind <- 1:length(PP$sp)
    if (!is.null(sp)) sp <- sp[-ind]
    if (!is.null(min.sp)) { 
      PP$min.sp <- min.sp[ind]
      min.sp <- min.sp[-ind]
    } 
  }
  
  # next work through smooth terms (if any) extending model matrix.....
  
  G$smooth <- list()
  G$S <- list()
 
  if (gamm.call) { ## flag that this is a call from gamm --- some smoothers need to know!
    if (m>0) for (i in 1:m) attr(split$smooth.spec[[i]],"gamm") <- TRUE
  }

  if (m>0 && idLinksBases) { ## search smooth.spec[[]] for terms linked by common id's
    id.list <- list() ## id information list
    for (i in 1:m) if (!is.null(split$smooth.spec[[i]]$id)) {
      id <- as.character(split$smooth.spec[[i]]$id)
      if (length(id.list)&&id%in%names(id.list)) { ## it's an existing id
        ni <- length(id.list[[id]]$sm.i) ## number of terms so far with this id
        id.list[[id]]$sm.i[ni+1] <- i    ## adding smooth.spec index to this id's list
        ## clone smooth.spec from base smooth spec....
        base.i <- id.list[[id]]$sm.i[1]
         
        split$smooth.spec[[i]] <- clone.smooth.spec(split$smooth.spec[[base.i]],
                                                      split$smooth.spec[[i]])
        
        ## add data for this term to the data list for basis setup...
        temp.term <- split$smooth.spec[[i]]$term
       
        ## note cbind deliberate in next line, as construction will handle matrix argument 
        ## correctly... 
        for (j in 1:length(temp.term)) id.list[[id]]$data[[j]] <- cbind(id.list[[id]]$data[[j]],
                                                          get.var(temp.term[j],data,vecMat=FALSE))
       
       } else { ## new id
        id.list[[id]] <- list(sm.i=i) ## start the array of indices of smooths with this id
        id.list[[id]]$data <- list()
        ## need to collect together all data for which this basis will be used,
        ## for basis setup...
        term <- split$smooth.spec[[i]]$term
        for (j in 1:length(term)) id.list[[id]]$data[[j]] <- get.var(term[j],data,vecMat=FALSE)
      } ## new id finished
    }
  } ## id.list complete

  G$off<-array(0,0)
  first.para<-G$nsdf+1
  sm <- list()
  newm <- 0
  if (m>0) for (i in 1:m) {
    # idea here is that terms are set up in accordance with information given in split$smooth.spec
    # appropriate basis constructor is called depending on the class of the smooth
    # constructor returns penalty matrices model matrix and basis specific information
  
    id <- split$smooth.spec[[i]]$id
    if (is.null(id)||!idLinksBases) { ## regular evaluation
      sml <- smoothCon(split$smooth.spec[[i]],data,knots,absorb.cons,scale.penalty=scale.penalty,
                       null.space.penalty=select,sparse.cons=sparse.cons,
                       diagonal.penalty=diagonal.penalty,apply.by=apply.by,modCon=modCon) 
    } else { ## it's a smooth with an id, so basis setup data differs from model matrix data
      names(id.list[[id]]$data) <- split$smooth.spec[[i]]$term ## give basis data suitable names
      sml <- smoothCon(split$smooth.spec[[i]],id.list[[id]]$data,knots,
                       absorb.cons,n=nrow(data),dataX=data,scale.penalty=scale.penalty,
                       null.space.penalty=select,sparse.cons=sparse.cons,
                       diagonal.penalty=diagonal.penalty,apply.by=apply.by,modCon=modCon)
    }
    for (j in 1:length(sml)) {
      newm <- newm + 1
      sm[[newm]] <- sml[[j]]
    }
  }
  
  G$m <- m <- newm ## number of actual smooths

  ## at this stage, it is neccessary to impose any side conditions required
  ## for identifiability
  if (m>0) { 
    sm <- gam.side(sm,X,tol=.Machine$double.eps^.5)
    if (!apply.by) for (i in 1:length(sm)) { ## restore any by-free model matrices
      if (!is.null(sm[[i]]$X0)) { ## there is a by-free matrix to restore 
        ind <- attr(sm[[i]],"del.index") ## columns, if any to delete
        sm[[i]]$X <- if (is.null(ind)) sm[[i]]$X0 else sm[[i]]$X0[,-ind,drop=FALSE] 
      }
    }
  }

  ## The matrix, L, mapping the underlying log smoothing parameters to the
  ## log of the smoothing parameter multiplying the S[[i]] must be
  ## worked out...
  idx <- list() ## idx[[id]]$c contains index of first col in L relating to id
  L <- matrix(0,0,0)
  lsp.names <- sp.names <- rep("",0) ## need a list of names to identify sps in global sp array
  if (m>0) for (i in 1:m) {
    id <- sm[[i]]$id
    ## get the L matrix for this smooth...
    length.S <- length(sm[[i]]$S)
    if (is.null(sm[[i]]$L)) Li <- diag(length.S) else Li <- sm[[i]]$L 
     
    if (length.S > 0) { ## there are smoothing parameters to name
       if (length.S == 1) lspn <- sm[[i]]$label else {
          Sname <- names(sm[[i]]$S)
          lspn <- if (is.null(Sname)) paste(sm[[i]]$label,1:length.S,sep="") else
                  paste(sm[[i]]$label,Sname,sep="") ## names for all sp's
       }
       spn <- lspn[1:ncol(Li)] ## names for actual working sps
    }

    ## extend the global L matrix...
    if (is.null(id)||is.null(idx[[id]])) { ## new `id'     
      if (!is.null(id)) { ## create record in `idx'
        idx[[id]]$c <- ncol(L)+1   ## starting column in L for this `id'
        idx[[id]]$nc <- ncol(Li)   ## number of columns relating to this `id'
      }
      L <- rbind(cbind(L,matrix(0,nrow(L),ncol(Li))),
                 cbind(matrix(0,nrow(Li),ncol(L)),Li))
      if (length.S > 0) { ## there are smoothing parameters to name
        sp.names <- c(sp.names,spn) ## extend the sp name vector
        lsp.names <- c(lsp.names,lspn) ## extend full.sp name vector
      }
    } else { ## it's a repeat id => shares existing sp's
      L0 <- matrix(0,nrow(Li),ncol(L))
      if (ncol(Li)>idx[[id]]$nc) {
        stop("Later terms sharing an `id' can not have more smoothing parameters than the first such term")
      }
      L0[,idx[[id]]$c:(idx[[id]]$c+ncol(Li)-1)] <- Li
      L <- rbind(L,L0)
      if (length.S > 0) { ## there are smoothing parameters to name
        lsp.names <- c(lsp.names,lspn) ## extend full.sp name vector
      }
    }
  }

  ## create the model matrix...

  Xp <- NULL ## model matrix under prediction constraints, if given
  if (m>0) for (i in 1:m) {
    n.para<-ncol(sm[[i]]$X)
    # define which elements in the parameter vector this smooth relates to....
    sm[[i]]$first.para<-first.para     
    first.para<-first.para+n.para
    sm[[i]]$last.para<-first.para-1
    ## termwise offset handling ...
    Xoff <- attr(sm[[i]]$X,"offset")
    if (!is.null(Xoff)) { 
      if (is.null(G$offset)) G$offset <- Xoff
      else G$offset <- G$offset + Xoff
    }
    ## model matrix accumulation ...
    
    ## alternative version under alternative constraint first (prediction only)
    if (is.null(sm[[i]]$Xp)) {
      if (!is.null(Xp)) Xp <- cbind2(Xp,sm[[i]]$X)
    } else { 
      if (is.null(Xp)) Xp <- X
      Xp <- cbind2(Xp,sm[[i]]$Xp);sm[[i]]$Xp <- NULL
    }
    ## now version to use for fitting ...
    X <- cbind2(X,sm[[i]]$X);sm[[i]]$X<-NULL
   
    G$smooth[[i]] <- sm[[i]]   
  }

  if (is.null(Xp)) {
    G$cmX <- colMeans(X) ## useful for componentwise CI construction 
  } else {
    G$cmX <- colMeans(Xp)
    ## transform from fit params to prediction params...
    ## G$P <- qr.coef(qr(Xp),X) ## old code assumes always full rank!!
    
    qrx <- qr(Xp,LAPACK=TRUE)
    R <- qr.R(qrx)
    p <- ncol(R)
    rank <- Rrank(R) ## rank of Xp/R    
    QtX <- qr.qty(qrx,X)[1:rank,]
    if (rank<p) { ## rank deficient  
      R <- R[1:rank,]
      qrr <- qr(t(R),tol=0)
      R <- qr.R(qrr)
      G$P <- forwardsolve(t(R),QtX)
    } else {
      G$P <- backsolve(R,QtX)
    }
    if (rank<p) {
      G$P <- qr.qy(qrr,rbind(G$P,matrix(0,p-rank,p)))
    }
    G$P[qrx$pivot,] <- G$P
  }
  ## cmX relates to computation of CIs incorportating uncertainty about the mean
  ## It may make more sense to incorporate all uncertainty about the mean,
  ## rather than just the uncertainty in the fixed effects mean. This means
  ## incorporating the mean of random effects and unconstrained smooths. Hence
  ## comment out the following.
  #if (G$nsdf>0) G$cmX[-(1:G$nsdf)] <- 0 ## zero the smooth parts here 
  #else G$cmX <- G$cmX * 0
  G$X <- X;rm(X)
  n.p <- ncol(G$X) 
  # deal with penalties


  ## min.sp must be length nrow(L) to make sense
  ## sp must be length ncol(L) --- need to partition
  ## L into columns relating to free log smoothing parameters,
  ## and columns, L0, corresponding to values supplied in sp.
  ## lsp0 = L0%*%log(sp[sp>=0]) [need to fudge sp==0 case by
  ## setting log(0) to log(effective zero) computed case-by-case]

  ## following deals with supplied and estimated smoothing parameters...
  ## first process the `sp' array supplied to `gam'...
  
  if (!is.null(sp)) { # then user has supplied fixed smoothing parameters
   ok <- TRUE 
   if (length(sp) < ncol(L)) { 
      warning("Supplied smoothing parameter vector is too short - ignored.")
      ok <- FALSE
    }
    if (sum(is.na(sp))) { 
      warning("NA's in supplied smoothing parameter vector - ignoring.")
      ok <- FALSE
    }
  } else ok <- FALSE
  G$sp <- if (ok) sp[1:ncol(L)] else rep(-1,ncol(L))
  
  names(G$sp) <- sp.names

  ## now work through the smooths searching for any `sp' elements
  ## supplied in `s' or `te' terms.... This relies on `idx' created 
  ## above...
  
  k <- 1 ## current location in `sp' array
  if (m>0) for (i in 1:m) {
    id <- sm[[i]]$id
    if (is.null(sm[[i]]$L)) Li <- diag(length(sm[[i]]$S)) else Li <- sm[[i]]$L 
    if (is.null(id)) { ## it's a smooth without an id
      spi <- sm[[i]]$sp
      if (!is.null(spi)) { ## sp supplied in `s' or `te'
        if (length(spi)!=ncol(Li)) stop("incorrect number of smoothing parameters supplied for a smooth term")
        G$sp[k:(k+ncol(Li)-1)] <- spi
      }       
      k <- k + ncol(Li) 
    } else { ## smooth has an id
      spi <- sm[[i]]$sp
      if (is.null(idx[[id]]$sp.done)) { ## not already dealt with these sp's
        if (!is.null(spi)) { ## sp supplied in `s' or `te'
          if (length(spi)!=ncol(Li)) stop("incorrect number of smoothing parameters supplied for a smooth term")
          G$sp[idx[[id]]$c:(idx[[id]]$c+idx[[id]]$nc-1)] <- spi
        }
        idx[[id]]$sp.done <- TRUE ## only makes sense to use supplied `sp' from defining term
        k <- k + idx[[id]]$nc 
      }
    }
  } ## finished processing `sp' vectors supplied in `s' or `te' terms

  ## copy initial sp's back into smooth objects, so there is a record of
  ## fixed and free...
  k <- 1 
  if (length(idx)) for (i in 1:length(idx)) idx[[i]]$sp.done <- FALSE
  if (m>0) for (i in 1:m) { ## work through all smooths
    id <- sm[[i]]$id 
    if (!is.null(id)) { ## smooth with id
      if (idx[[id]]$nc>0) { ## only copy if there are sp's
        G$smooth[[i]]$sp <- G$sp[idx[[id]]$c:(idx[[id]]$c+idx[[id]]$nc-1)]
      }   
      if (!idx[[id]]$sp.done) { ## only update k on first encounter with this smooth
        idx[[id]]$sp.done <- TRUE
        k <- k + idx[[id]]$nc
      }
    
    } else { ## no id, just work through sp 
      if (is.null(sm[[i]]$L)) nc <- length(sm[[i]]$S) else nc <- ncol(sm[[i]]$L)
      if (nc>0) G$smooth[[i]]$sp <- G$sp[k:(k+nc-1)]
      k <- k + nc
    }
  } ## now all elements of G$smooth have a record of initial sp. 


  if (!is.null(min.sp)) { # then minimum s.p.'s supplied
    if (length(min.sp)<nrow(L)) stop("length of min.sp is wrong.")
    min.sp <- min.sp[1:nrow(L)]
    if (sum(is.na(min.sp))) stop("NA's in min.sp.")
    if (sum(min.sp<0)) stop("elements of min.sp must be non negative.")
  }

  k.sp <- 0 # count through sp and S
  G$rank <- array(0,0)
  if (m>0) for (i in 1:m) {
    sm<-G$smooth[[i]]
    if (length(sm$S)>0)
    for (j in 1:length(sm$S)) {  # work through penalty matrices
      k.sp <- k.sp+1
      G$off[k.sp] <- sm$first.para 
      G$S[[k.sp]] <- sm$S[[j]]
      G$rank[k.sp]<-sm$rank[j]
      if (!is.null(min.sp)) {
        if (is.null(H)) H<-matrix(0,n.p,n.p)
        H[sm$first.para:sm$last.para,sm$first.para:sm$last.para] <-
        H[sm$first.para:sm$last.para,sm$first.para:sm$last.para]+min.sp[k.sp]*sm$S[[j]] 
      }           
    } 
  }
 
  ## need to modify L, lsp.names, G$S, G$sp, G$rank and G$off to include any penalties
  ## on parametric stuff, at this point....
  if (!is.null(PP)) { ## deal with penalties on parametric terms
    L <- rbind(cbind(L,matrix(0,nrow(L),ncol(PP$L))),
                 cbind(matrix(0,nrow(PP$L),ncol(L)),PP$L))
    G$off <- c(PP$off,G$off)
    G$S <- c(PP$S,G$S)
    G$rank <- c(PP$rank,G$rank)
    G$sp <- c(PP$sp,G$sp)
    lsp.names <- c(PP$full.sp.names,lsp.names)
    G$n.paraPen <- length(PP$off)
    if (!is.null(PP$min.sp)) { ## deal with minimum sps
      if (is.null(H)) H <- matrix(0,n.p,n.p)
      for (i in 1:length(PP$S)) {
        ind <- PP$off[i]:(PP$off[i]+ncol(PP$S[[i]])-1)
        H[ind,ind] <- H[ind,ind] + PP$min.sp[i] * PP$S[[i]]
      }
    } ## min.sp stuff finished
  } else G$n.paraPen <- 0


  ## Now remove columns of L and rows of sp relating to fixed 
  ## smoothing parameters, and use removed elements to create lsp0

  fix.ind <- G$sp>=0

  if (sum(fix.ind)) {
    lsp0 <- G$sp[fix.ind]
    ind <- lsp0==0 ## find the zero s.p.s
    ef0 <- indi <- (1:length(ind))[ind]
    if (length(indi)>0) for (i in 1:length(indi)) {
      ## find "effective zero" to replace each zero s.p. with
      ii <- G$off[i]:(G$off[i]+ncol(G$S[[i]])-1) 
      ef0[i] <- norm(G$X[,ii],type="F")^2/norm(G$S[[i]],type="F")*.Machine$double.eps*.1
    }
    lsp0[!ind] <- log(lsp0[!ind])
    lsp0[ind] <- log(ef0) ##log(.Machine$double.xmin)*1000 ## zero fudge
    lsp0 <- as.numeric(L[,fix.ind,drop=FALSE]%*%lsp0)

    L <- L[,!fix.ind,drop=FALSE]  
    G$sp <- G$sp[!fix.ind]
  } else {lsp0 <- rep(0,nrow(L))}

  G$H <- H

  if (ncol(L)==nrow(L)&&!sum(L!=diag(ncol(L)))) L <- NULL ## it's just the identity

  G$L <- L;G$lsp0 <- lsp0
  names(G$lsp0) <- lsp.names ## names of all smoothing parameters (not just underlying)

  if (absorb.cons==FALSE) {  ## need to accumulate constraints 
    G$C <- matrix(0,0,n.p)
    if (m>0) {
      for (i in 1:m) {
        if (is.null(G$smooth[[i]]$C)) n.con<-0 
        else n.con<- nrow(G$smooth[[i]]$C)
        C <- matrix(0,n.con,n.p)
        C[,G$smooth[[i]]$first.para:G$smooth[[i]]$last.para]<-G$smooth[[i]]$C
        G$C <- rbind(G$C,C)
        G$smooth[[i]]$C <- NULL
      }
      rm(C)
    }
  } ## absorb.cons == FALSE
 
  G$y <- data[[split$response]]
         
  ##data[[deparse(split$full.formula[[2]],backtick=TRUE)]]
  
  G$n <- nrow(data)

  if (is.null(data$"(weights)")) G$w <- rep(1,G$n)
  else G$w <- data$"(weights)"  

  ## Create names for model coefficients... 

  if (G$nsdf > 0) term.names <- colnames(G$X)[1:G$nsdf] else term.names<-array("",0)
  n.smooth <- length(G$smooth)
  if (n.smooth)
  ## create coef names, if smooth has any coefs, and create a global indicator of non-linear parameters
  ## g.index, if needed
  for (i in 1:n.smooth) {
    k <- 1
    jj <- G$smooth[[i]]$first.para:G$smooth[[i]]$last.para
    if (G$smooth[[i]]$df > 0) for (j in jj) {
      term.names[j] <- paste(G$smooth[[i]]$label,".",as.character(k),sep="")
      k <- k+1
    }
    if (!is.null(G$smooth[[i]]$g.index)) {
      if (is.null(G$g.index)) G$g.index <- rep(FALSE,n.p)
      G$g.index[jj] <- G$smooth[[i]]$g.index
    } 
  }
  G$term.names <- term.names

  ## Deal with non-linear parameterizations...
 

  G$pP <- PP ## return paraPen object, if present

  G
} ## gam.setup


clone.smooth.spec <- function(specb,spec) {
## produces a version of base smooth.spec, `specb', but with 
## the variables relating to `spec'. Used by `gam.setup' in handling 
## of linked smooths.
 ## check dimensions same...
 if (specb$dim!=spec$dim) stop("`id' linked smooths must have same number of arguments") 
 ## Now start cloning...
 if (inherits(specb,c("tensor.smooth.spec","t2.smooth.spec"))) { ##`te' or `t2' generated base smooth.spec
    specb$term <- spec$term
    specb$label <- spec$label 
    specb$by <- spec$by
    k <- 1
    for (i in 1:length(specb$margin)) {
      if (is.null(spec$margin)) { ## sloppy user -- have to construct margin info...
         for (j in 1:length(specb$margin[[i]]$term)) {
           specb$margin[[i]]$term[j] <- spec$term[k]
           k <- k + 1
         }
         specb$margin[[i]]$label <- ""
 
      } else { ## second term was at least `te'/`t2', so margin cloning is easy
        specb$margin[[i]]$term <- spec$margin[[i]]$term
        specb$margin[[i]]$label <- spec$margin[[i]]$label
        specb$margin[[i]]$xt <- spec$margin[[i]]$xt
      }
    }

  } else { ## `s' generated case
    specb$term <- spec$term
    specb$label <- spec$label 
    specb$by <- spec$by
    specb$xt <- spec$xt ## don't generally know what's in here => don't clone
  }
  specb ## return clone
} ## clone.smooth.spec


all.vars1 <- function(form) {
## version of all.vars that doesn't split up terms like x$y into x and y
  vars <- all.vars(form)
  vn <- all.names(form)
  vn <- vn[vn%in%c(vars,"$","[[")] ## actual variable related names
  if ("[["%in%vn) stop("can't handle [[ in formula")
  ii <- which(vn%in%"$") ## index of '$'
  if (length(ii)) { ## assemble variable names
    vn1 <- if (ii[1]>1) vn[1:(ii[1]-1)]
    go <- TRUE
    k <- 1
    while (go) {
      n <- 2; 
      while(k<length(ii) && ii[k]==ii[k+1]-1) { k <- k + 1;n <- n + 1 }
      vn1 <- c(vn1,paste(vn[ii[k]+1:n],collapse="$"))
      if (k==length(ii)) {
        go <- FALSE
	ind <- if (ii[k]+n<length(vn)) (ii[k]+n+1):length(vn) else rep(0,0) 
      } else {
        k <- k +  1
	ind <- if (ii[k-1]+n<ii[k]-1) (ii[k-1]+n+1):(ii[k]-1) else rep(0,0)
      }
      vn1 <- c(vn1,vn[ind])
    }
  } else vn1 <- vn
  vn1
} ## all.vars1


variable.summary <- function(pf,dl,n) {
## routine to summarize all the variables in dl, which is a list
## containing raw input variables to a model (i.e. no functions applied)
## pf is a formula containing the strictly parametric part of the
## model for the variables in dl. A list is returned, with names given by 
## the variables. For variables in the parametric part, then the list elements
## may be:
## * a 1 column matrix with elements set to the column medians, if variable 
##   is a matrix.
## * a 3 figure summary (min,median,max) for a numeric variable.
## * a factor variable, with the most commonly occuring factor (all levels)
## --- classes are as original data type, but anything not numeric, factor or matrix
## is coerced to numeric. 
## For non-parametric variables, any matrices are coerced to numeric, otherwise as 
## parametric.      
## medians in the above are always observed values (to deal with variables coerced to 
## factors in the model formulae in a nice way).
## variables with less than `n' entries are discarded
   v.n <- length(dl)
   ## if (v.n) for (i in 1:v.n) if (length(dl[[i]])<n) dl[[i]] <- NULL 
   
   v.name <- v.name1 <- names(dl)
   if (v.n) ## need to strip out names of any variables that are too short.
   { k <- 0 ## counter for retained variables
     for (i in 1:v.n) if (length(dl[[i]])>=n) { 
       k <- k+1
       v.name[k] <- v.name1[i] ## save names of variables of correct length
     }
     if (k>0) v.name <- v.name[1:k] else v.name <- rep("",k)
   }

   ## v.name <- names(dl)    ## the variable names
   p.name <- all.vars(pf[-2]) ## variables in parametric part (not response)
   vs <- list()
   v.n <- length(v.name)
   if (v.n>0) for (i in 1:v.n) {
     if (v.name[i]%in%p.name) para <- TRUE else para <- FALSE ## is variable in the parametric part?

     if (para&&is.matrix(dl[[v.name[i]]])&&ncol(dl[[v.name[i]]])>1) { ## parametric matrix --- a special case
       x <- matrix(apply(dl[[v.name[i]]],2,quantile,probs=0.5,type=3,na.rm=TRUE),1,ncol(dl[[v.name[i]]])) ## nearest to median entries
     } else { ## anything else
       x <- dl[[v.name[i]]]
       if (is.character(x)) x <- as.factor(x)
       if (is.factor(x)) {
         x <- x[!is.na(x)]
         lx <- levels(x)
         freq <- tabulate(x)
         ii <- min((1:length(lx))[freq==max(freq)])
         x <- factor(lx[ii],levels=lx) 
       } else {
         x <- as.numeric(x)
         x <- c(min(x,na.rm=TRUE),as.numeric(quantile(x,probs=.5,type=3,na.rm=TRUE)) ,max(x,na.rm=TRUE)) ## 3 figure summary
       }
     }
     vs[[v.name[i]]] <- x
   }
   vs
} ## end variable.summary



parametricPenalty <- function(pterms,assign,paraPen,sp0) {
## routine to process any penalties on the parametric part of the model.
## paraPen is a list whose items have names corresponding to the 
## term.labels in pterms. Each list item may have named elements 
## L, rank and sp. All other elements should be penalty coefficient matrices.
  S <- list()     ## penalty matrix list
  off <- rep(0,0) ## offset array
  rank <- rep(0,0) ## rank array
  sp <- rep(0,0)    ## smoothing param array
  full.sp.names <- rep("",0) ## names for sp's multiplying penalties (not underlying)
  L <- matrix(0,0,0) 
  k <- 0
  tind <- unique(assign) ## unique term indices
  n.t <- length(tind)
  if (n.t>0) for (j in 1:n.t) if (tind[j]>0) {
    term.label <- attr(pterms[tind[j]],"term.label")
    P <- paraPen[[term.label]] ## get any penalty information for this term
    if (!is.null(P)) { ## then there is information
      ind <- (1:length(assign))[assign==tind[j]] ## index of coefs involved here
      Li <- P$L;P$L <- NULL
      spi <- P$sp;P$sp <- NULL
      ranki <- P$rank;P$rank <- NULL
      ## remaining terms should be penalty matrices...
      np <- length(P)

      if (!is.null(ranki)&&length(ranki)!=np) stop("`rank' has wrong length in `paraPen'") 
      if (np) for (i in 1:np) { ## unpack penalty matrices, offsets and ranks
        k <- k + 1
        S[[k]] <- P[[i]]
        off[k] <- min(ind) ## index of first coef penalized by this term
        if ( ncol(P[[i]])!=nrow(P[[i]])||nrow(P[[i]])!=length(ind)) stop(" a parametric penalty has wrong dimension")
        if (is.null(ranki)) {
          ev <- eigen(S[[k]],symmetric=TRUE,only.values=TRUE)$values
          rank[k] <- sum(ev>max(ev)*.Machine$double.eps*10) ## estimate rank
        } else rank[k] <- ranki[i]
      }
      ## now deal with L matrices
      if (np) { ## only do this stuff if there are any penalties!
        if (is.null(Li)) Li <- diag(np)
        if (nrow(Li)!=np) stop("L has wrong dimension in `paraPen'")
        L <- rbind(cbind(L,matrix(0,nrow(L),ncol(Li))),
                   cbind(matrix(0,nrow(Li),ncol(L)),Li))
        ind <- (length(sp)+1):(length(sp)+ncol(Li))
        ind2 <- (length(sp)+1):(length(sp)+nrow(Li)) ## used to produce names for full sp array
        if (is.null(spi)) {
          sp[ind] <- -1 ## auto-initialize
        } else {
          if (length(spi)!=ncol(Li)) stop("`sp' dimension wrong in `paraPen'")
          sp[ind] <- spi
        }
        ## add smoothing parameter names....
        if (length(ind)>1) names(sp)[ind] <- paste(term.label,ind-ind[1]+1,sep="") 
        else names(sp)[ind] <- term.label
        
        if (length(ind2)>1) full.sp.names[ind2] <- paste(term.label,ind2-ind2[1]+1,sep="") 
        else full.sp.names[ind2] <- term.label
      }
    } ## end !is.null(P)  
  } ## looped through all terms
  if (k==0) return(NULL)
  if (!is.null(sp0)) {
    if (length(sp0)<length(sp)) stop("`sp' too short")
    sp0 <- sp0[1:length(sp)]
    sp[sp<0] <- sp0[sp<0]
  }
  ## S is list of penalty matrices, off[i] is index of first coefficient penalized by each S[[i]]
  ## sp is array of underlying smoothing parameter (-ve to estimate), L is matrix mapping log
  ## underlying smoothing parameters to log smoothing parameters, rank[i] is the rank of S[[i]].
  list(S=S,off=off,sp=sp,L=L,rank=rank,full.sp.names=full.sp.names)
} ## parametricPenalty


## copied from bam(mgcv)..
## (c) Simon N. Wood 2009-2019

rwMatrix <- function(stop,row,weight,X,trans=FALSE) {
## Routine to recombine the rows of a matrix X according to info in 
## stop, row and weight. Consider the ith row of the output matrix 
## ind <- 1:stop[i] if i==1 and ind <- (stop[i-1]+1):stop[i]
## otherwise. The ith output row is then X[row[ind],]*weight[ind]
  if (is.matrix(X)) { n <- nrow(X);p<-ncol(X);ok <- TRUE} else { n<- length(X);p<-1;ok<-FALSE}
  stop <- stop - 1;row <- row - 1 ## R indices -> C indices
  oo <-.C(C_rwMatrix,as.integer(stop),as.integer(row),as.double(weight),X=as.double(X),
          as.integer(n),as.integer(p),trans=as.integer(trans),work=as.double(rep(0,n*p)))
  if (ok) return(matrix(oo$X,n,p)) else
  return(oo$X) 
} ## rwMatrix

AR.resid <- function(rsd,rho=0,AR.start=NULL) {
## standardised residuals for AR1 model
  if (rho==0) return(rsd)
  ld <- 1/sqrt(1-rho^2) ## leading diagonal of root inverse correlation
  sd <- -rho*ld         ## sub diagonal
  N <- length(rsd)    
  ## see rwMatrix() for how following are used...
  ar.row <- c(1,rep(1:N,rep(2,N))[-c(1,2*N)]) ## index of rows to reweight
  ar.weight <- c(1,rep(c(sd,ld),N-1))     ## row weights
  ar.stop <- c(1,1:(N-1)*2+1)    ## (stop[i-1]+1):stop[i] are the rows to reweight to get ith row
  if (!is.null(AR.start)) { ## need to correct the start of new AR sections...
    ii <- which(AR.start==TRUE)
    if (length(ii)>0) {
          if (ii[1]==1) ii <- ii[-1] ## first observation does not need any correction
          ar.weight[ii*2-2] <- 0 ## zero sub diagonal
          ar.weight[ii*2-1] <- 1 ## set leading diagonal to 1
    }
  }
  rwMatrix(ar.stop,ar.row,ar.weight,rsd)
} ## AR.resid

###############################################################
## loading functions, copied from mgcv() package of Simon Wood
#################################################################

print.scam.version <- function()
{ library(help=scam)$info[[1]] -> version
  version <- version[pmatch("Version",version)]
  um <- strsplit(version," ")[[1]]
  version <- um[nchar(um)>0][2]
  hello <- paste("This is scam ",version,".",sep="")
  packageStartupMessage(hello)
}


.onAttach <- function(...) { 
  print.scam.version()
}

##.onUnload <- function(libpath) library.dynam.unload("scam", libpath)








