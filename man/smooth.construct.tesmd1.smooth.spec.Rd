\name{smooth.construct.tesmd1.smooth.spec}
%\Rdversion{1.0}
\alias{smooth.construct.tesmd1.smooth.spec}
\title{Tensor product smoothing constructor for a bivariate function monotone decreasing
       in the first covariate
}
\description{This is a special method function
  for creating tensor product bivariate smooths monotone decreasing in the first covariate which is built by 
  the \code{mgcv} constructor function for smooth terms, \code{smooth.construct}. 
  It is constructed from a pair of single penalty marginal smooths. 
  This tensor product is specified by model terms such as \code{s(x1,x2,k=c(q1,q2),bs="tesmd1",m=2)}. 
  The default basis for the second marginal smooth is P-spline. Cyclic cubic regression spline (\code{"cc"}) is implemented in addition to the P-spline. See an example below on how to call for it.   
}
\usage{
\method{smooth.construct}{tesmd1.smooth.spec}(object, data, knots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A smooth specification object, generated by an \code{s} term in a GAM formula.} 

  \item{data}{A data frame or list containing the values of the elements of \code{object$term},
     with names given by \code{object$term}.}
 
  \item{knots}{An optional list containing the knots corresponding to \code{object$term}. 
          If it is \code{NULL} then the knot locations are generated automatically.}
 
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
An object of class \code{"tesmd1.smooth"}. In addition to the usual 
elements of a smooth class documented under \code{smooth.construct} of the \code{mgcv} library, 
this object contains:
\item{p.ident}{A vector of 0's and 1's for model parameter identification: 
     1's indicate parameters which will be exponentiated, 0's - otherwise.}
\item{Zc}{A matrix of identifiability constraints.}
\item{margin.bs}{A two letter character string indicating the (penalized) smoothing basis to use for the second unconstrained marginal smooth. (eg \code{"ps"} for P-splines).}

}
\references{
Pya, N. and Wood, S.N. (2015) Shape constrained additive models. Statistics and Computing, 25(3), 543-559

Pya, N. (2010) Additive models with shape constraints. PhD thesis. University of Bath. Department of Mathematical Sciences

}
\author{ 
   Natalya Pya <nat.pya@gmail.com>

}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{smooth.construct.tesmd2.smooth.spec}}
}
\examples{
 \dontrun{
## tensor product `tesmd1' example 
  ## simulating data...
   require(scam)
   set.seed(2)
   n <- 30
   x1 <- sort(runif(n)*4-1); x2 <- sort(runif(n))
   f <- matrix(0,n,n)
   for (i in 1:n) for (j in 1:n) 
        f[i,j] <- -exp(4*x1[i])/(1+exp(4*x1[i]))+2*sin(pi*x2[j])
   f <- as.vector(t(f))
   y <- f+rnorm(length(f))*.2
   x11 <-  matrix(0,n,n)
   x11[,1:n] <- x1
   x11 <- as.vector(t(x11))
   x22 <- rep(x2,n)
   dat <- list(x1=x11,x2=x22,y=y)
## fit model ...
   b <- scam(y~s(x1,x2,bs="tesmd1",k=10),data=dat)
## plot results ...
   old.par <- par(mfrow=c(2,2),mar=c(4,4,2,2))
   plot(b,se=TRUE)
   plot(b,pers=TRUE,theta = 30, phi = 40)
   plot(y,b$fitted.values,xlab="Simulated data",ylab="Fitted data")
   par(old.par)
   vis.scam(b,theta=40,phi=20)

## example with cyclic cubic regression spline along the second covariate...
   set.seed(2)
   n <- 30
   x1 <- sort(runif(n)*4-1); x2 <- sort(runif(n))
   f <- matrix(0,n,n)
   for (i in 1:n) for (j in 1:n) 
        f[i,j] <- -exp(4*x1[i])/(1+exp(4*x1[i]))+sin(2*pi*x2[j])
   f <- as.vector(t(f))
   y <- f+rnorm(length(f))*.2
   x11 <-  matrix(0,n,n)
   x11[,1:n] <- x1
   x11 <- as.vector(t(x11))
   x22 <- rep(x2,n)
   dat <- list(x1=x11,x2=x22,y=y)
## fit model ...
   b1 <- scam(y~s(x1,x2,bs="tesmd1",xt=list("cc"),k=10), data=dat)
## plot results ...
   old.par <- par(mfrow=c(2,2))
   plot(b1,se=TRUE)
   plot(b1,pers=TRUE,theta = 30, phi = 40)
   plot(y,b1$fitted.values,xlab="Simulated data",ylab="Fitted data")
   par(old.par)
   vis.scam(b1,theta=40,phi=20)
  }
} 
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models} \keyword{regression}%-- one or more ..

