\name{smooth.construct.tecvcv.smooth.spec}
\alias{smooth.construct.tecvcv.smooth.spec}

\title{Tensor product smoothing constructor for bivariate function subject to double concavity constraint
}
\description{This is a special method function
  for creating tensor product bivariate smooths subject to double concavity constraint, i.e. concavity constraint wrt both the first and the second covariates. This is built by 
  the \code{mgcv} constructor function for smooth terms, \code{smooth.construct}. 
  It is constructed from a pair of single penalty marginal smooths which are represented using the B-spline basis functions.  
  This tensor product is specified by model terms such as \code{s(x1,x2,k=c(q1,q2),bs="tecvcv",m=c(2,2))}, 
  where \code{q1} and \code{q2} denote the basis dimensions for the marginal smooths.
}
\usage{
\method{smooth.construct}{tecvcv.smooth.spec}(object, data, knots)
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
\value{
An object of class \code{"tecvcv.smooth"}. In addition to the usual 
elements of a smooth class documented under \code{smooth.construct} of the \code{mgcv} library, 
this object contains:
\item{p.ident}{A vector of 0's and 1's for model parameter identification: 
     1's indicate parameters which will be exponentiated, 0's - otherwise.}
\item{Zc}{A matrix of identifiability constraints.}

}
\references{
Pya, N. and Wood, S.N. (2015) Shape constrained additive models. Statistics and Computing, 25(3), 543-559


}
\author{ 
   Natalya Pya <nat.pya@gmail.com>

}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{smooth.construct.tedmd.smooth.spec}}
\code{\link{smooth.construct.temicx.smooth.spec}}

\code{\link{smooth.construct.tedecx.smooth.spec}}
\code{\link{smooth.construct.tecxcx.smooth.spec}}

\code{\link{smooth.construct.tecxcv.smooth.spec}}
}
\examples{
\donttest{
## tensor product `tecvcv' example 
require(scam)
simu <- function(x,z) {-4*(x^2+z^2)   } 
xs <- zs <- seq(-1,1,length=30)
pr <- data.frame(x=rep(xs,30),z=rep(zs,rep(30,30)))
truth <- matrix(simu(pr$x,pr$z),30,30)
   
set.seed(3)
n <- 500
x <- 2*runif(n)-1
z <- 2*runif(n)-1
f <- simu(x,z)
y <- f + rnorm(n)*.1
## fit model ...
b <- scam(y~s(x,z,k=c(10,10),bs="tecvcv"))
summary(b)
old.par <- par(mfrow=c(2,2),mar=c(4,4,2,2))
plot(b,se=TRUE)
plot(b,pers=TRUE,theta = 30, phi = 40);title("tecvcv")
plot(y,b$fitted.values,xlab="Simulated data",ylab="Fitted data",pch=".",cex=3)
persp(xs,zs,truth,theta = 30, phi = 40);title("truth")
par(old.par)  
vis.scam(b,theta = 30, phi = 40) 
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models} \keyword{regression}%-- one or more ..


