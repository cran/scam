\name{linear.functional.terms}
\alias{linear.functional.terms}
\alias{function.predictors}
\alias{signal.regression}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Linear functionals of a smooth in GAMs}

\description{Since \code{scam} uses the model setup  of \code{gam} of the \code{mgcv} package, in the same way as in \code{gam} \code{\link{scam}} allows the response variable to depend on linear 
functionals of smooth terms in the s with additional shape constraints.

See \code{linear.functional.terms(mgcv)}.

}

%- maybe also `usage' for other objects documented here.

%\author{ Simon N. Wood \email{simon.wood@r-project.org}}

\examples{

###########################################
## similar to a "signal" regression 
## example from mgcv() ...
###########################################
library(scam)
## decreasing smooth...
set.seed(4)
rf <- function(x=seq(-1,3,length=100)) { ## from mgcv package
## generates random functions...
  m <- ceiling(runif(1)*5) ## number of components
  f <- x*0;
  mu <- runif(m,min(x),max(x)); sig <- (runif(m)+.5)*(max(x)-min(x))/10
  for (i in 1:m) f <- f+ dnorm(x,mu[i],sig[i])
  f
}

## simulate 200 functions and store in rows of L...
L <- matrix(NA,200,100) 
for (i in 1:200) L[i,] <- rf()  ## simulate the functional predictors

x <- seq(-1,3,length=100) ## evaluation points
f2 <- function(x) { ## the coefficient function
     -4*exp(4*x)/(1+exp(4*x))  
}
f <- f2(x) 
plot(x,f ,type="l")
y <- L\%*\%f + rnorm(200)*20 ## simulated response data
X <- matrix(x,200,100,byrow=TRUE) 

b <- scam(y~s(X,by=L,k=20,bs="mpdBy")) 
par(mfrow=c(1,2))
plot(b,shade=TRUE);lines(x,f,col=2); 
## compare with gam() of mgcv package...
require(mgcv)
g <- gam(y~s(X,by=L,k=20)) 
plot(g,shade=TRUE);lines(x,f,col=2)


## increasing smooth....
L <- matrix(NA,200,100) 
for (i in 1:200) L[i,] <- rf()  ## simulate the functional predictors
x <- seq(-1,3,length=100) ## evaluation points
f2 <- function(x) { ## the coefficient function
     4*exp(4*x)/(1+exp(4*x))  
}
f <- f2(x) 
plot(x,f ,type="l")
y <- L\%*\%f + rnorm(200)*20 ## simulated response data
X <- matrix(x,200,100,byrow=TRUE) 
b <- scam(y~s(X,by=L,k=20,bs="mpiBy")) 
par(mfrow=c(1,2))
plot(b,shade=TRUE);lines(x,f,col=2); 
## compare with unconstrained fit...
g <- scam(y~s(X,by=L,k=20)) 
plot(g,shade=TRUE);lines(x,f,col=2)


## convex smooth...
 ## simulate 200 functions and store in rows of L...
set.seed(4)
L <- matrix(NA,200,100) 
for (i in 1:200) L[i,] <- rf(x=sort(2*runif(100)-1))  ## simulate the functional predictors

x <- sort(runif(100,-1,1)) ## evaluation points
f2 <- function(x){4*x^2 } ## the coefficient function
f <- f2(x) 
plot(x,f ,type="l")
y <- L\%*\%f + rnorm(200)*30 ## simulated response data
X <- matrix(x,200,100,byrow=TRUE) 

b <- scam(y~s(X,by=L,k=20,bs="cxBy")) 
par(mfrow=c(1,2))
plot(b,shade=TRUE);lines(x,f,col=2); 

g <- scam(y~s(X,by=L,k=20)) 
plot(g,shade=TRUE);lines(x,f,col=2)

 
}
\keyword{models} \keyword{regression}%-- one or more ..


