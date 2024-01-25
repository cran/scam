## (c) Natalya Pya (2012-2024). Provided under GPL 2.
## routines for checking scam()...
## based on gam.check routines (c) Simon N Wood
## (2024) qq.scam() added as a short version of qq.gam (c) Simon N Wood


scam.check <- function(b,type=c("deviance","pearson","response"),old.style=FALSE, pch=".",
                       ## arguments passed to qq.scam():
		       rep=0, level=.9, rl.col=3, rep.col="gray80", ...)
## takes a fitted scam object and produces some standard diagnostic plots
{   ## old.par<-par(mfrow=c(2,2))
    if (is.null(.Platform$GUI) || .Platform$GUI != "RStudio") old.par <- par(mfrow=c(2,2))
    sc.name<-b$method
    type <- match.arg(type)
    ylab <- paste(type,"residuals") 
    resid <- residuals(b,type=type)

    if (old.style){
       qqnorm(resid,pch=pch,ylab=ylab,...)
       qqline(resid,col=rl.col,...)
    } else
       qq.scam(b, rep=rep, level=level, type=type, rl.col=rl.col, rep.col=rep.col, ...)
    
    plot(b$linear.predictors,resid,main="Resids vs. linear pred.",
         xlab="linear predictor",ylab="residuals",...);
    hist(resid,xlab="Residuals",main="Histogram of residuals",...);
    plot(fitted(b),b$y,xlab="Fitted Values",ylab="Response",main="Response vs. Fitted Values",...)
    
    ## now summarize convergence information 
    cat("\nMethod:", b$method, "  Optimizer:", b$optimizer)
    if (b$optimizer[1] == "optim"){ 
        cat("\nOptim Method:", b$optim.method[1])
        if (is.na(b$optim.method[2]))
              cat("\n Finite-difference approximation of the GCV/UBRE gradient was used.")
    }
    if (!is.null(b$bfgs.info)) { ## summarize BFGS convergence information
        boi <- b$bfgs.info
        cat("\nNumber of iterations of smoothing parameter selection performed was",boi$iter,".")
        cat("\n",boi$conv,".",sep="")
        cat("\nGradient range: [",min(boi$grad),",",max(boi$grad),"]",sep="")
        cat("\n(score ",formatC(b$gcv.ubre, digits = 5)," & scale ",formatC(b$sig2, digits = 5),")",sep="")
     }
     else if (!is.null(b$optim.info)) { ## summarize optim() convergence information
        boi <- b$optim.info
        cat("\nNumber of iterations of smoothing parameter selection performed was",boi$iter[1],".")
        cat("\n",boi$conv,".",sep="")
        cat("\n(score ",formatC(b$gcv.ubre, digits = 5)," & scale ",formatC(b$sig2, digits = 5),")",sep="")
     }
     else if (!is.null(b$nlm.info)) { ## summarize nlm() convergence information
        boi <- b$nlm.info
        cat("\nNumber of iterations of smoothing parameter selection performed was",boi$iter,".")
        cat("\n",boi$conv,".",sep="")
        cat("\nGradient range: [",min(boi$grad),",",max(boi$grad),"]",sep="")
        cat("\n(score ",formatC(b$gcv.ubre, digits = 5)," & scale ",formatC(b$sig2, digits = 5),")",sep="")
     }
     else if (!is.null(b$efs.info)) { ## summarize 'efs' convergence information
        boi <- b$efs.info
        cat("\nNumber of iterations of smoothing parameter selection performed was",boi$iter[1],".")
        cat("\n",boi$conv,".",sep="")
        cat("\n(score ",formatC(b$gcv.ubre, digits = 5)," & scale ",formatC(b$sig2, digits = 5),")",sep="")
     }
     else {
       if (length(b$sp)==0) ## no sp's estimated  
        cat("\nModel required no smoothing parameter selection")
     }
    ## print the estimated smoothing parameters...
    if (length(b$sp)!=0)
        cat("\nThe optimal smoothing parameter(s):",round(b$sp,5),".")
    cat("\n")
    par(old.par)
}


qq.scam <- function(object, rep=0, level=.9,s.rep=10,
                   type=c("deviance","pearson","response"),
                   pch=".", rl.col=3, rep.col="gray80",...) 
{
## get deviance residual quantiles under good fit
  type <- match.arg(type)
  ylab <- paste(type,"residuals")

  if (inherits(object,c("glm","scam"))) {
    if (is.null(object$sig2)) object$sig2 <- summary(object)$dispersion
  } else stop("object is not a glm or scam")

  ## in case of NA & na.action="na.exclude", we need the "short" residuals:
  object$na.action <- NULL
  D <- residuals(object,type=type)

 
  lim <- Dq <- NULL
  if (rep==0) { 
    fam <- fix.family.qf(object$family)
    if (is.null(fam$qf))
      rep <- 50 ## try simulation if quantile function not available
    level <- 0
  } 
  n <- length(D)
  if (rep > 0) { ## simulate quantiles
    fam <- fix.family.rd(object$family)
    if (!is.null(fam$rd)) {
      ##d <- rep(0,0)
      ## simulate deviates... 
      dm <- matrix(0,n,rep)
      for (i in 1:rep) { 
        yr <- fam$rd(object$fitted.values, object$prior.weights, object$sig2)
        #di <- fam$dev.resids(yr,object$fitted.values,object$prior.weights)^.5*
        #       sign(yr-object$fitted.values)
        object$y <- yr
        dm[,i] <- sort(residuals(object,type=type))
        #d <- c(d,sort(di))
      }
      # n <- length(D)
      Dq <- quantile(as.numeric(dm),(1:n - .5)/n) 
    
      ## now get simulation limits on QQ plot
      #dm <- matrix(d,length(Dq),rep)
      alpha <- (1-level)/2
      if (alpha>.5||alpha<0) alpha <- .05
      if (level>0&&level<1) lim <- apply(dm,1,FUN=quantile,p=c(alpha,1-alpha)) else
      if (level >= 1) lim <- level 
    }
  } else {
    U <- (1:n-.5)/n
    if (!is.null(fam$qf)) { 
      dm <- matrix(0,n,s.rep)
      for (i in 1:s.rep) { 
        U <- sample(U,n) ## randomize uniform quantiles w.r.t. obs
        q0 <- fam$qf(U,object$fitted.values,object$prior.weights,object$sig2)
        object$y <- q0
        dm[,i] <- sort(residuals(object,type=type)) ## original proposal
      }
      Dq <- sort(rowMeans(dm))
    }
  }
 
  if (!is.null(Dq))  
  { qqplot(Dq,D,ylab=ylab,xlab="theoretical quantiles",ylim=range(c(lim,D)),
           pch=pch,...)
    abline(0,1,col=rl.col)
    if (!is.null(lim)) {
      if (level>=1) for (i in 1:rep) lines(Dq,dm[,i],col=rep.col) else {
        n <- length(Dq)
        polygon(c(Dq,Dq[n:1],Dq[1]),c(lim[1,],lim[2,n:1],lim[1,1]),col=rep.col,border=NA)
      }
      abline(0,1,col=rl.col) 
    }
    points(Dq,sort(D),pch=pch,...)
    return(invisible(Dq))
  } else qqnorm(D,ylab=ylab,pch=pch,...)
} ## qq.scam


