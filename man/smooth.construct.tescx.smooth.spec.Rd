\name{smooth.construct.tescx.smooth.spec}
\alias{smooth.construct.tescx.smooth.spec}

\title{Tensor product smoothing constructor for a bivariate function convex
        in the second covariate
}
\description{This is a special method function
  for creating tensor product bivariate smooths convex in the second covariate which is built by the \code{mgcv} constructor function for smooth terms, \code{smooth.construct}. 
  It is constructed from a pair of single penalty
  marginal smooths. This tensor product is specified by model terms such as \code{s(x1,x2,k=c(q1,q2),bs="tescx",m=c(2,2))}, 
  where the basis for the first marginal smooth is specified in the second element of \code{bs}.

}
\usage{
\method{smooth.construct}{tescx.smooth.spec}(object, data, knots)
}
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
\value{An object of class \code{"tescx.smooth"}. In addition to the usual 
elements of a smooth class documented under \code{smooth.construct} of the \code{mgcv} library, this object contains:
\item{p.ident}{A vector of 0's and 1's for model parameter identification: 
     1's indicate parameters which will be exponentiated, 0's - otherwise.}
\item{Zc}{A matrix of identifiability constraints.}
%\item{margin.bs}{A two letter character string indicating the (penalized) smoothing basis to use for the first unconstrained marginal smooth. (eg \code{"ps"} for P-splines).}

}
\references{
Pya, N. and Wood, S.N. (2015) Shape constrained additive models. Statistics and Computing, 25(3), 543-559

}
\author{
Natalya Pya <nat.pya@gmail.com>

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{smooth.construct.temicv.smooth.spec}}
  \code{\link{smooth.construct.temicx.smooth.spec}}

  \code{\link{smooth.construct.tedecv.smooth.spec}}
  \code{\link{smooth.construct.tedecx.smooth.spec}}

  \code{\link{smooth.construct.tescv.smooth.spec}}
}
\examples{
 
## tensor product `tescx' example 
require(scam)
simu <- function(x,z) { sin(x) + 2*z^2  } 
xs <-seq(0,1,length=30); zs <- seq(-1,1,length=30)
pr <- data.frame(x=rep(xs,30),z=rep(zs,rep(30,30)))
truth <- matrix(simu(pr$x,pr$z),30,30)
   
set.seed(5)
n <- 500
x <- runif(n)
z <- 2*runif(n)-1
f <- simu(x,z)
y <- f + rnorm(n)*.2
## fit model ...
b <- scam(y~s(x,z,bs="tescx"))
summary(b)
old.par <- par(mfrow=c(2,2),mar=c(4,4,2,2))
plot(b,se=TRUE)
plot(b,pers=TRUE,theta = 50, phi = 20);title("tescx")
plot(y,b$fitted.values,xlab="Simulated data",ylab="Fitted data",pch=".",cex=3)
persp(xs,zs,truth,theta = 50, phi = 20);title("truth")
par(old.par)  
vis.scam(b,theta = 50, phi = 20) 

}

\keyword{models} \keyword{regression}%-- one or more ..
