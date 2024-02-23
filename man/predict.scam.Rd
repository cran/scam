\name{predict.scam}
\alias{predict.scam}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Prediction from fitted SCAM model}
\description{ 
This function is a clone of the \code{mgcv} library code \code{\link{predict.gam}} with some modifications 
to adopt shape preserving smooth terms.
It takes a fitted \code{scam} object produced by \code{scam()} 
and produces predictions given a new set of values for the model covariates 
or the original values used for the model fit. Predictions can be accompanied
by standard errors, based on the posterior distribution of the model
coefficients. 

It now alows prediction outside the range of knots, and use linear extrapolation in this case.

%Also for extrapolation of the response variable values in case of a single univariate shape constrained %%term see \code{\link{extrapolate.uni.scam}}. 
}

\usage{
\method{predict}{scam}(object,newdata,type="link",se.fit=FALSE,terms=NULL,exclude=NULL,
    block.size=NULL,newdata.guaranteed=FALSE,na.action=na.pass,...)

}
%- maybe also `usage' for other objects documented here.
\arguments{ 
  
 \item{object}{ a fitted \code{scam} object as produced by \code{scam()}.
                }
 \item{newdata}{ A data frame or list containing the values of the model covariates at which predictions
  are required. If this is not provided then predictions corresponding to the
                original data are returned. If \code{newdata} is provided then
                it should contain all the variables needed for prediction: a
                warning is generated if not. } 
 
\item{type}{ When this has the value \code{"link"} (default) the linear predictor (possibly with
associated standard errors) is returned. When \code{type="terms"} each component of the 
linear predictor is returned seperately (possibly with standard errors): this includes 
parametric model components, followed by each smooth component, but excludes
any offset and any intercept. \code{type="iterms"} is the same, except that any standard errors 
returned for unconstrained smooth components will include the uncertainty about the intercept/overall mean.  When 
\code{type="response"} predictions 
on the scale of the response are returned (possibly with approximate
standard errors). When \code{type="lpmatrix"} then a matrix is returned
which yields the values of the linear predictor (minus any offset) when
postmultiplied by the
parameter vector (in this case \code{se.fit} is ignored). The latter
option is most useful for getting variance estimates for quantities derived from
the model: for example integrated quantities, or derivatives of smooths. A
linear predictor matrix can also be used to implement approximate prediction
outside \code{R} (see example code, below). }

\item{se.fit}{ when this is TRUE (not default) standard error estimates are returned for each prediction.}

\item{terms}{if \code{type=="terms"} then only results for the terms given in this array
will be returned.}

\item{exclude}{if \code{type=="terms"} or \code{type="iterms"} then terms (smooth or parametric) named
 in this array will not be returned. Otherwise any smooth terms named in this
 array will be set to zero. If \code{NULL} then no terms are excluded.}

\item{block.size}{maximum number of predictions to process per call to underlying
code: larger is quicker, but more memory intensive. Set to < 1 to use total number
of predictions as this.}

\item{newdata.guaranteed}{Set to \code{TRUE} to turn off all checking of
\code{newdata} except for sanity of factor levels: this can speed things up
for large prediction tasks, but \code{newdata} must be complete, with no
\code{NA} values for predictors required in the model. }

\item{na.action}{what to do about \code{NA} values in \code{newdata}. With the
default \code{na.pass}, any row of \code{newdata} containing \code{NA} values
for required predictors, gives rise to \code{NA} predictions (even if the term concerned has no
\code{NA} predictors). \code{na.exclude} or \code{na.omit} result in the
dropping of \code{newdata} rows, if they contain any \code{NA} values for
required predictors. If \code{newdata} is missing then \code{NA} handling is 
determined from \code{object$na.action}.}

\item{...}{ other arguments.}

}

\value{ If \code{type=="lpmatrix"} then a matrix is returned which will
  give a vector of linear predictor values (minus any offest) at the supplied covariate
  values, when applied to the model coefficient vector. 
Otherwise, if \code{se.fit} is \code{TRUE} then a 2 item list is returned with items (both arrays) \code{fit}
and \code{se.fit} containing predictions and associated standard error estimates, otherwise an 
array of predictions is returned. The dimensions of the returned arrays depends on whether 
\code{type} is \code{"terms"} or not: if it is then the array is 2 dimensional with each 
term in the linear predictor separate, otherwise the array is 1 dimensional and contains the 
linear predictor/predicted values (or corresponding s.e.s). The linear predictor returned termwise will 
not include the offset or the intercept.

\code{newdata} can be a data frame, list or model.frame: if it's a model frame
then all variables must be supplied.

}

\details{See \code{\link{predict.gam}} for details.
}

\references{

Chambers and Hastie (1993) Statistical Models in S. Chapman & Hall.

Wood S.N. (2006) Generalized Additive Models: An Introduction with R. Chapman and Hall/CRC Press.

Pya, N. and Wood, S.N. (2015) Shape constrained additive models. Statistics and Computing, 25(3), 543-559

Pya, N. (2010) Additive models with shape constraints. PhD thesis. University of Bath. Department of Mathematical Sciences

}
\author{ 
Natalya Pya <nat.pya@gmail.com> based partly on \code{mgcv} by Simon Wood

}

\seealso{  \code{\link{scam}}, \code{\link{plot.scam}}}

\examples{
\dontrun{
library(scam)
set.seed(2)
n <- 200
x1 <- runif(n)*6-3
f1 <- 3*exp(-x1^2) # unconstrained term
x2 <- runif(n)*4-1;
f2 <- exp(4*x2)/(1+exp(4*x2)) # monotone increasing smooth
f <- f1+f2
y <- f+rnorm(n)*0.2
dat <- data.frame(x1=x1,x2=x2,y=y)
b <- scam(y~s(x1,k=15,bs="cr")+s(x2,k=30,bs="mpi"),data=dat)

newd <- data.frame(x1=seq(-3,3,length.out=20),x2=seq(-1,3,length.out=20))
pred <- predict(b,newd)
pred
predict(b,newd,type="terms",se=TRUE)

## prediction on the original data...
pr <- predict(b,type="terms")
x<-sort(x1,index=TRUE)
old.par <- par(mfrow=c(2,2))
plot(x$x,(pr[,1])[x$ix],type="l",col=3,xlab="x1")
z<-sort(x2,index=TRUE)
plot(z$x,(pr[,2])[z$ix],type="l",col=3,xlab="x2")
plot(b,select=1,scale=0,se=FALSE)
plot(b,select=2,scale=0,se=FALSE)
par(old.par)

## linear extrapolation with predict.scam()...
set.seed(3)
n <- 100
x <- sort(runif(n)*3-1)
f <- exp(-1.3*x)
y <- rpois(n,exp(f))
dat <- data.frame(x=x,y=y)
b <- scam(y~s(x,k=15,bs="mpd"),family=poisson(link="log"),data=dat)
newd <- data.frame(x=c(2.3,2.7,3.2))
fe <- predict(b,newd,type="link",se=TRUE)
ylim<- c(min(y,exp(fe$fit)),max(y,exp(fe$fit)))
plot(c(x,newd[[1]]),c(y,NA,NA,NA),ylim=ylim,ylab="y",xlab="x")
lines(c(x,newd[[1]]),c(b$fitted,exp(fe$fit)),col=3)

## prediction on the original data...
pr <- predict(b)
plot(x,y)
lines(x,exp(pr),col=3)


## Gaussian model ....
## simulating data...
 set.seed(2)
 n <- 200
 x <- sort(runif(n)*4-1)
 f <- exp(4*x)/(1+exp(4*x)) # monotone increasing smooth
 y <- f+rnorm(n)*0.1
 dat <- data.frame(x=x,y=y)
 b <- scam(y~ s(x,k=25,bs="mpi"),data=dat)
 newd <- data.frame(x=c(3.2,3.3,3.6))
 fe <- predict(b,newd)
 plot(c(x,newd[[1]]),c(y,NA,NA,NA),ylab="y",xlab="x")
 lines(c(x,newd[[1]]),c(b$fitted,fe),col=3)

### passing observed data + new data...
 newd <- data.frame(x=c(x,3.2,3.3,3.6))
 fe <- predict(b,newd,se=TRUE)
 plot(newd[[1]],c(y,NA,NA,NA),ylab="y",xlab="x")
 lines(newd[[1]],fe$fit,col=2)
 lines(newd[[1]],fe$fit+2*fe$se.fit,col=3)
 lines(newd[[1]],fe$fit-2*fe$se.fit,col=4)

## prediction with CI...
 newd <- data.frame(x=seq(-1.2,3.5,length.out=100))
 fe <- predict(b,newd,se=TRUE)
 ylim<- c(min(y,fe$se.fit),max(y,fe$se.fit))
 plot(newd[[1]],fe$fit,type="l",ylim=ylim,ylab="y",xlab="x")
 lines(newd[[1]],fe$fit+2*fe$se.fit,lty=2)
 lines(newd[[1]],fe$fit-2*fe$se.fit,lty=2)
 
## prediction on the original data...
pr <- predict(b)
plot(x,y)
lines(x,pr,col=3)

## bivariate example...
   set.seed(2)
   n <- 30
   x1 <- sort(runif(n));  x2 <- sort(runif(n)*4-1)
   f <- matrix(0,n,n)
   for (i in 1:n) for (j in 1:n) 
        f[i,j] <- 2*sin(pi*x1[i]) +exp(4*x2[j])/(1+exp(4*x2[j]))
   f <- as.vector(t(f));  
   y <- f+rnorm(length(f))*0.1
   x11 <-  matrix(0,n,n); x11[,1:n] <- x1;  x11 <- as.vector(t(x11))
   x22 <- rep(x2,n)
   dat <- list(x1=x11,x2=x22,y=y)
   b <- scam(y~s(x1,x2,k=c(10,10),bs="tesmi2"),data=dat,optimizer="efs")
   par(mfrow=c(2,2),mar=c(4,4,2,2))
   plot(b,se=TRUE);   plot(b,pers=TRUE,theta = 80, phi = 40)

   n.out <- 20
   xp <- seq(0,1.4,length.out=n.out) 
   zp <- seq(-1,3.4,length.out=n.out)
   xp1 <-  matrix(0,n.out,n.out);   xp1[,1:n.out] <- xp
   xp1 <- as.vector(t(xp1));   xp2 <- rep(zp,n.out)
   newd <- data.frame(x1=xp1,x2=xp2)
   fe <- predict(b,newd)
   fc <- t(matrix(fe,n.out,n.out))
   persp(xp,zp,fc,expand= 0.85,ticktype = "simple",xlab="x1",
     ylab="x2",zlab="f^",main="", theta = 80, phi = 40)


## obtaining a 'prediction matrix'...
newd <- data.frame(x1=c(-2,-1),x2=c(0,1))
Xp <- predict(b,newdata=newd,type="lpmatrix")
fv <- Xp\%*\% b$beta.t
fv
 }
}

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..

