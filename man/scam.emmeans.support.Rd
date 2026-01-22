\name{scam.emmeans.support}
\alias{scam.emmeans.support}
\alias{recover_data.scam}
\title{Support Functions for \pkg{emmeans}}
\usage{
recover_data.scam(object, ...)
}
\arguments{
\item{object}{A fitted model object of the class \code{scam} as produced by \code{scam()}.}

\item{...}{Additional parameters that may be supported by the method.}

}
\description{
Function required for compatibility of \pkg{scam} with the \pkg{emmeans} package.
The function is not normally called directly, but rather service routines for \code{emmeans} function
of the \pkg{emmeans} package.
}
\references{
Searle, S. R., Speed, F. M., and Milliken, G. A. (1980). Population Marginal Means in the Linear Model: An Alternative to Least Squares Means. The American Statistician, 34(4), 216-221. https://doi.org/10.1080/00031305.1980.10483031
}

\examples{ 
\donttest{
require(scam)
require(emmeans)
set.seed(2)
n <- 100
x <- sort(runif(n,-.5,.5))
f <- -20*x^2
fac <- as.factor(sample(1:3,n,replace=TRUE))
y <- ifelse(fac == 1, 4, ifelse(fac == 2, 1, -2)) + 
      f + rnorm(n, 0, 1.5)
dat <- data.frame(y=y,x=x,fac=fac)
m <- scam(y ~ fac+s(x,bs="cv"), data=dat)
em <- emmeans(m, ~ fac)
em
pairs(em)
 }
}

\keyword{models} 
\keyword{factor} 
\keyword{regression}
\keyword{internal}
