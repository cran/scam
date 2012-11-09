\name{extrapolate.uni.scam}
\alias{extrapolate.uni.scam}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Extrapolation from the fitted SCAM model}
\description{ 
This function does extrapolations of the response variable values to the right-hand side of the data in case of a single univariate shape constrained term.
It takes a fitted \code{scam} object produced by \code{scam()} 
and produces extrapolations given a new set of values for the model covariates to the ``right'' of the data set by adding knots and giving the same increment to parameter vector which results in linear extrapolation. It also produces predictions given a new set of values within the range of the model covariate or the original values used for the model fit.
Extrapolations are accompanied by standard errors, based on the posterior distribution of the model coefficients.

For scam predictions and linear extrapolations see \code{\link{predict.scam}}. 

} 

\usage{
extrapolate.uni.scam(object, newdata)
}
\arguments{ 
  
 \item{object}{ a fitted \code{scam} object as produced by \code{scam()}.
                }
 \item{newdata}{ A data frame or list containing the values of the model covariates at which extrapolations
  are required.} 
}

\value{ 
A list is returned with the following items:

  \item{fit}{Extrapolations of the response variable.}
  \item{se}{standard error estimates of extrapolated values.}

}

\note{
  Extrapolation for "micv" and "mdcx" will flatten the smooth, the very short forecasting is recommended (parameters extensions differ from those for other smooths due to the Sigma matrix construction)
}

\references{

Pya, N. (2010) \emph{Additive models with shape constraints}. PhD thesis. University of Bath. Department of Mathematical Sciences

}
\author{ 
Natalya Pya <nat.pya@gmail.com> 
}

\seealso{  \code{\link{scam}}, \code{\link{predict.scam}}}

\examples{
\dontrun{
 library(scam)
#########################################################
## Example on extrapolating a single univariate monotone increasing smooth...
##########################################################

set.seed(12)
   n <- 100
   x <- sort(runif(n)*4-1)
   f <- 4*exp(4*x)/(1+exp(4*x))
   y <- rpois(n,exp(f))
   dat <- data.frame(x=x,y=y)

   b <- scam(y~s(x,k=15,bs="mpi",m=2),family=poisson(link="log"),
       data=dat,sp=NULL)

   plot(x,y,xlab="x",ylab="y")
   lines(x,exp(f))      ## the true function
   lines(x,b$fitted.values,col=2)  ## monotone fit

   newd <- data.frame(x=c(3.2,3.3,3.6))
   fe <- extrapolate.uni.scam(b,newd)

   plot(c(x,newd[[1]]),c(y,NA,NA,NA))
   lines(c(x,newd[[1]]),c(b$fitted,exp(fe$fit)),col=2)

### passing observed data + new data...
   newd <- data.frame(x=c(x,3.2,3.3,3.6))
   fe <- extrapolate.uni.scam(b,newd)

   plot(newd[[1]],c(y,NA,NA,NA))
   lines(newd[[1]],exp(fe$fit),col=2)


## Gaussian model ....
   ## simulating data...
 set.seed(2)
 n <- 200
 x <- sort(runif(n)*4-1)
 f <- exp(4*x)/(1+exp(4*x)) # monotone increasing smooth
 y <- f+rnorm(n)*0.1
 dat <- data.frame(x=x,y=y)
 b <- scam(y~ s(x,k=25,bs="mpi",m=2),family=gaussian(link="identity"),data=dat)

 newd <- data.frame(x=c(3.2,3.3,3.6))
 fe <- extrapolate.uni.scam(b,newd)

 plot(c(x,newd[[1]]),c(y,NA,NA,NA))
 lines(c(x,newd[[1]]),c(b$fitted,fe$fit),col=2)

### passing observed data + new data...
   newd <- data.frame(x=c(x,3.2,3.3,3.6))
   fe <- extrapolate.uni.scam(b,newd)

   plot(newd[[1]],c(y,NA,NA,NA))
   lines(newd[[1]],fe$fit,col=2)
   lines(newd[[1]],fe$fit+2*fe$se,col=3)
   lines(newd[[1]],fe$fit-2*fe$se,col=4)

###########################################################
## Example on extrapolating a single univariate monotone increasing and concave smooth...
##########################################################
  
set.seed(3)
   n <- 100
   x <- sort(runif(n)*99+1)
   f1 <- log(x)
   f <- (f1-min(f1))/(max(f1)-min(f1))
   y <- f+rnorm(n)*0.10
   dat <- data.frame(x=x,y=y)
 ## fit model ...
   b <- scam(y~s(x,k=15,bs="micv",m=2),family=gaussian(link="identity"),data=dat)

   newd <- data.frame(x=c(100,101,113))
   fe <- extrapolate.uni.scam(b,newd)

   plot(c(x,newd[[1]]),c(y,NA,NA,NA),ylim=c(min(y),max(c(y,fe$fit))))
   lines(c(x,newd[[1]]),c(b$fitted,fe$fit),col=2)
 
###########################################################
## Example on extrapolating a single univariate monotone decreasing smooth...
##########################################################

set.seed(3)
n <- 100
x <- sort(runif(n)*3-1)
f <- exp(-1.3*x)
y <- rpois(n,exp(f))
dat <- data.frame(x=x,y=y)
b <- scam(y~s(x,k=15,bs="mpd",m=2),family=poisson(link="log"),
       data=dat,sp=NULL)


newd <- data.frame(x=c(2.3,2.7,3.2))
fe <- extrapolate.uni.scam(b,newd)

ylim<- c(min(y,exp(fe$fit)),max(y,exp(fe$fit)))
plot(c(x,newd[[1]]),c(y,NA,NA,NA),ylim=ylim)
lines(c(x,newd[[1]]),c(b$fitted,exp(fe$fit)),col=2)

### passing observed data + new data...
newd <- data.frame(x=c(x,2.3,2.7,3.2))
fe <- extrapolate.uni.scam(b,newd)

plot(newd[[1]],c(y,NA,NA,NA),ylim=ylim)
lines(newd[[1]],exp(fe$fit),col=2)

## Gaussian model ....
set.seed(3)
n <- 100
x <- sort(runif(n)*3-1)
f <- exp(-1.3*x)
y <- f+rnorm(n)*0.3
dat <- data.frame(x=x,y=y)
b <- scam(y~s(x,k=15,bs="mpd",m=2),family=gaussian,
       data=dat,sp=NULL)


newd <- data.frame(x=c(2.3,2.7,3.2))
fe <- extrapolate.uni.scam(b,newd)

ylim<- c(min(y,fe$fit),max(y,fe$fit))
plot(c(x,newd[[1]]),c(y,NA,NA,NA),ylim=ylim)
lines(c(x,newd[[1]]),c(b$fitted,fe$fit),col=2)


###########################################################
## Example on extrapolating a single univariate convex smooth...
##########################################################
  set.seed(1)
   n <- 100
   x <- sort(2*runif(n)-1)
   f <- 4*x^2
   y <- f + rnorm(n)*0.4
   dat <- data.frame(x=x,y=y)
   b <- scam(y~s(x,k=15,bs="cx",m=2),family=gaussian,data=dat)

newd <- data.frame(x=c(1.1,1.2,1.3))
fe <- extrapolate.uni.scam(b,newd)

ylim<- c(min(y,fe$fit),max(y,fe$fit))
plot(c(x,newd[[1]]),c(y,NA,NA,NA),ylim=ylim)
lines(c(x,newd[[1]]),c(b$fitted,fe$fit),col=2)

### passing observed data + new data...
newd <- data.frame(x=c(x,1.1,1.2,1.3))
fe <- extrapolate.uni.scam(b,newd)

plot(newd[[1]],c(y,NA,NA,NA),ylim=ylim)
lines(newd[[1]],fe$fit,col=2)

 }
}

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..

