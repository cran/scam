\name{derivative.scam}
\alias{derivative.scam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Derivative of the univariate smooth model terms
}
\description{Function to get derivatives of the smooth model terms (currently only of the univariate smooths).
 Analytical derivatives for SCOP-splines, finite difference approximation is used for all others                                   
}
\usage{
derivative.scam(object,smooth.number=1,deriv=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{fitted scam object} 

  \item{smooth.number}{ordered number of the smooth model term (1,2,...), ordered as in the formula, which derivative is needed to be calculated.}
 
  \item{deriv}{either 1 if the 1st derivative is required, or 2 if the 2nd}

}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\item{d}{values of the derivative of the smooth term.} 

\item{se.d}{standard errors of the derivative.}
}
\references{

Pya, N. and Wood, S.N. (2014) Shape constrained additive models. Statistics and Computing, DOI \url{10.1007/s11222-013-9448-7}

}
\author{ 
   Natalya Pya <nat.pya@gmail.com>
}
%% ~Make other sections like Warning with \section{Warning}{....} ~

\seealso{\code{\link{scam}}
}
\examples{

\dontrun{
 set.seed(2)
n <- 200
x1 <- runif(n)*4-1;
f1 <- exp(4*x1)/(1+exp(4*x1)) # monotone increasing smooth
x2 <- sort(runif(n)*3-1)      # decreasing smooth
f2 <- exp(-1.3*x2)
f <- f1+ f2 
y <- f+ rnorm(n)*0.2
## fit model, results, and plot...
b <- scam(y~ s(x1,k=20,bs="mpi")+s(x2,k=15,bs="mpd"),family=gaussian)

d1 <- derivative.scam(b,smooth.number=1,deriv=1)

par(mfrow=c(1,2))

xx <- sort(x1,index=TRUE)
plot(xx$x,d1$d[xx$ix],type="l",xlab=expression(x[1]),
     ylab=expression(df[1]/dx[1]))

d2 <- derivative.scam(b,smooth.number=2,deriv=1)

xx <- sort(x2,index=TRUE)
plot(xx$x,d2$d[xx$ix],type="l",xlab=expression(x[2]),
     ylab=expression(df[2]/dx[2]))
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models} \keyword{regression}%-- one or more ..






