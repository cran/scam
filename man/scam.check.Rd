\name{scam.check}
%\Rdversion{1.0}
\alias{scam.check}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Some diagnostics for a fitted scam object
}
\description{Takes a fitted \code{scam} object produced by \code{scam()} and produces some diagnostic information
about the fitting procedure and results. This function is almost the same as \code{\link[mgcv]{gam.check}} of the \code{mgcv} library. The default is to produce four residual plots and some information about the
convergence of the smoothness selection optimization.
}
\usage{
scam.check(b,type=c("deviance","pearson","response"),old.style=FALSE, pch=".",
                       rep=0, level=.9, rl.col=3, rep.col="gray80",\dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{b}{ a fitted \code{scam} object as produced by \code{scam()}.}
 \item{old.style}{produces qq-norm plots as it was in scam versions < 1.2-15 when set to \code{TRUE}.}
  \item{type}{type of residuals, see \code{\link{residuals.scam}}, used in
    all plots.}
  \item{rep, level, rep.col}{
    arguments passed to \code{\link{qq.scam}()} when \code{old.style} is
    \code{FALSE} (default).} 
  \item{rl.col}{color for the reference line on the quantile-quantile plot.}
  \item{pch}{plot character to use for the quantile-quantile plot.}
  \item{\dots}{extra graphics parameters to pass to plotting functions.}

}
\details{As for \code{mgcv(gam)} plots 4 standard diagnostic plots, and some other
convergence diagnostics. The
printed information relates to the optimization process used to select smoothing
parameters.
}

\references{
Wood S.N. (2006) Generalized Additive Models: An Introduction with R. Chapman and Hall/CRC Press.

}
\author{

Natalya Pya \email{nat.pya@gmail.com} based partly on \code{mgcv} by Simon N Wood
}
\seealso{  \code{\link{scam}}
}

\examples{
library(scam)
set.seed(2)
n <- 200
x1 <- runif(n)*4-1;
f1 <- exp(4*x1)/(1+exp(4*x1)) # monotone increasing smooth
x2 <- runif(n)*3-1;
f2 <- exp(-3*x2)/15  # monotone decreasing and convex smooth
f <- f1+f2
y <- f+ rnorm(n)*0.2
dat <- data.frame(x1=x1,x2=x2,y=y)
b <- scam(y~ s(x1,k=25,bs="mpi")+s(x2,k=25,bs="mdcx"),data=dat)
plot(b,pages=1)
scam.check(b)
 
}

\keyword{models} 
\keyword{smooth} 
\keyword{regression} %-- one or more ..

