\name{smooth.construct.cv.smooth.spec}
%\Rdversion{1.1-2}
\alias{smooth.construct.cv.smooth.spec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Constructor for concave P-splines in SCAMs
}
\description{This is a special method function
  for creating smooths subject to concavity constraint which is built by 
  the \code{mgcv} constructor function for smooth terms, \code{smooth.construct}. 
  It is constructed using concave P-splines. This smooth is specified via model terms such as 
  \code{s(x,k,bs="cv",m=2)}, 
  where \code{k} denotes the basis dimension and \code{m+1} is the order of the B-spline basis.

}
\usage{
smooth.construct.cv.smooth.spec(object, data, knots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A smooth specification object, generated by an \code{s} term in a GAM formula.} 

  \item{data}{A data frame or list containing the data required by this term,
     with names given by \code{object$term}. The \code{by} variable is the last element.}
 
  \item{knots}{An optional list containing the knots supplied for basis setup.  
          If it is \code{NULL} then the knot locations are generated automatically.}

}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{An object of class \code{"cv.smooth"}. 
}
\references{
Pya, N. (2010) \emph{Additive models with shape constraints}. PhD thesis. University of Bath. Department of Mathematical Sciences

}
\author{ 
   Natalya Pya <nat.pya@gmail.com>

}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{smooth.construct.cx.smooth.spec}}, 
\code{\link{smooth.construct.mpi.smooth.spec}}, \code{\link{smooth.construct.mdcv.smooth.spec}}, 
\code{\link{smooth.construct.mdcx.smooth.spec}}, \code{\link{smooth.construct.micx.smooth.spec}}, 
\code{\link{smooth.construct.mpd.smooth.spec}}

}
\examples{
 \dontrun{
## Concave P-splines example 
  ## simulating data...
   set.seed(1)
   n <- 100
   x <- sort(2*runif(n)-1)
   f <- -4*x^2
   y <- f + rnorm(n)*0.45
   dat <- data.frame(x=x,y=y)
   b <- scam(y~s(x,k=15,bs="cv",m=2),family=gaussian,data=dat)
   # UNCONSTRAINED FIT *****************
   b1 <- scam(y~s(x,k=15,bs="cr",m=2),family=gaussian, data=dat)

## plot results ...
   plot(x,y,xlab="x",ylab="y")
   lines(x,f)      ## the true function
   lines(x,b$fitted,col=2)  ## constrained fit 
   lines(x,b1$fitted,col=3) ## unconstrained fit 

## Poisson version...
   y <- rpois(n,15*exp(f))
   dat <- data.frame(x=x,y=y)
 ## fit model ...
   b <- scam(y~s(x,k=15,bs="cv",m=2),family=poisson(link="log"),data=dat)

# UNCONSTRAINED FIT *****************
   b1 <- scam(y~s(x,k=15,bs="cr",m=2),family=poisson(link="log"), data=dat)

## plot results ...
   plot(x,y,xlab="x",ylab="y")
   lines(x,15*exp(f))      ## the true function
   lines(x,b$fitted,col=2)  ## constrained fit 
   lines(x,b1$fitted,col=3) ## unconstrained fit 

## plotting on log scale...
   plot(x,log(15*exp(f)),type="l")      ## the true function
   lines(x,log(b$fitted),col=2)  ## constrained fit 
   lines(x,log(b1$fitted),col=3) ## unconstrained fit 
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models} \keyword{regression}%-- one or more ..







