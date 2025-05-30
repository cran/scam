\name{smooth.construct.tesmi2.smooth.spec}
%\Rdversion{1.0}
\alias{smooth.construct.tesmi2.smooth.spec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Tensor product smoothing constructor for a bivariate function monotone increasing
        in the second covariate
}
\description{This is a special method function
  for creating tensor product bivariate smooths monotone increasing in the second covariate which is built by 
  the \code{mgcv} constructor function for smooth terms, \code{smooth.construct}. 
  It is constructed from a pair of single penalty
  marginal smooths. This tensor product is specified by model terms such as \code{s(x1,x2,k=c(q1,q2),bs="tesmi2",m=c(2,2))}. The default basis for the first marginal smooth is P-spline. Cyclic cubic regression spline (\code{"cc"}) is implemented in addition to the P-spline. See an example below on how to call for it. 
}
\usage{
\method{smooth.construct}{tesmi2.smooth.spec}(object, data, knots)
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
\value{An object of class \code{"tesmi2.smooth"}. In addition to the usual 
elements of a smooth class documented under \code{smooth.construct} of the \code{mgcv} library, 
this object contains:
\item{p.ident}{A vector of 0's and 1's for model parameter identification: 
     1's indicate parameters which will be exponentiated, 0's - otherwise.}
\item{Zc}{A matrix of identifiability constraints.}
\item{margin.bs}{A two letter character string indicating the (penalized) smoothing basis to use for the first unconstrained marginal smooth. (eg \code{"cc"} for cyclic cubic regression spline).}


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
 \code{\link{smooth.construct.tesmi1.smooth.spec}}
}
\examples{
\donttest{
## tensor product `tesmi2' example 
require(scam)
simu <- function(x,z) {2*sin(pi*x) +exp(4*z)/(1+exp(4*z)) }
xs <- seq(0,1,length=30); zs <- seq(-1,3,length=30)
pr <- data.frame(x=rep(xs,30),z=rep(zs,rep(30,30)))
truth <- matrix(simu(pr$x,pr$z),30,30)
   
set.seed(24)
n <- 500
x <- runif(n)
z <- runif(n)*4-1
f <- simu(x,z)
y <- f + rnorm(n)*.3
## fit model ...
b <- scam(y~s(x,z,bs="tesmi2",k=c(10,10)))
old.par <- par(mfrow=c(2,2))
persp(xs,zs,truth,theta=50,phi=20);title("truth")
vis.scam(b,theta=50,phi=20);title("tesmi2")
plot(b,se=TRUE)
plot(y,b$fitted.values,xlab="Simulated data",ylab="Fitted data")
par(old.par)

## example with cyclic cubic regression spline along the first covariate...
simu2 <- function(x,z) {sin(2*pi*x)+ exp(4*z)/(1+exp(4*z))}
xs <- seq(0,1,length=30); zs <- seq(-1,3,length=30)
pr <- data.frame(x=rep(xs,30),z=rep(zs,rep(30,30)))
truth2 <- matrix(simu2(pr$x,pr$z),30,30)
    
set.seed(2)
n <- 500
x <- runif(n)
z <- runif(n)*4-1
f <- simu2(x,z)
y <- f + rnorm(n)*.3
## fit model ...
b1 <- scam(y~s(x,z,bs="tesmi2",xt=list("cc"),k=10))
old.par <- par(mfrow=c(2,2))
plot(b1,se=TRUE)
vis.scam(b1,theta=50,phi=20);title("tesmi2, cyclic")
plot(y,b1$fitted.values,xlab="Simulated data",ylab="Fitted data")
persp(xs,zs,truth2,theta = 50, phi = 20);title("truth")
par(old.par)    
} 
}

\keyword{models} \keyword{regression}%-- one or more ..
