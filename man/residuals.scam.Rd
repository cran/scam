\name{residuals.scam}
\alias{residuals.scam}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{SCAM residuals}
\description{This function is a copy of the \code{mgcv} library code \code{residuals.gam}, but with an additional type of residuals, (randomized) quantile residuals, implemented. It returns residuals for a fitted \code{scam} model
  object. Pearson, deviance, working, response, and (randomized) quantile residuals are available. 
}

\usage{
\method{residuals}{scam}(object, type =  c("deviance", "pearson","scaled.pearson", "working", 
                        "response", "rquantile"), setseed=NULL, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{object}{ a \code{scam} fitted model object. }
  \item{type}{the type of residuals.  Usually one of \code{"deviance"}, \code{"pearson"}, \code{"scaled.pearson"}, 
\code{"working"}, \code{"response"}, or \code{"rquantile"}.}
  \item{setseed}{seed to be used in getting randomized quantile residuals.}
  \item{...}{other arguments.}
}
\details{
\code{"deviance"}: deviance residuals defined by the model family.

\code{"scaled.pearson"}: scaled Pearson residuals are raw residuals (data minus fitted values) divided by the standard deviation of the data according to the model mean variance relationship and estimated scale parameter.

\code{"pearson"}:  Pearson residuals are the same as scaled Pearson, but multiplied by the square root of the scale parameter (so they are independent of the scale parameter): \eqn{(y-\mu)/\sqrt{V(\mu)}}{(y-\mu)/\sqrt{V(\mu)}}, where
\eqn{y}{y} is data, \eqn{\mu}{\mu} is model fitted value, and \eqn{V}{V} is model mean-variance relationship.

\code{"working"}: working residuals are the residuals returned from model fitting at convergence.

\code{"response"}: response residuals are the raw residuals (data minus fitted values).

 \code{"rquantile"}: (randomized) quantile residuals (RQRs) as suggested by Dunn and Smyth (1996). Based on the idea of inverting the estimated distribution function for each observation to get approximately normally distributed residuals. RQRs are the theoretical standard normal quantiles evaluated at the estimated CDF at each observation. For models with a continuous response distribution they are simply quantile residuals; for discrete distributions, one realization of randomized quantile residuals is produced. In case of normal responses, quantile residuals are equivalent to deviance and Pearson residuals. In non-normal situations, RQRs might be a better choice for performing model checking. 
}
\value{ An array of residuals.
}

\references{
 Dunn, P. K. and Smyth, G. K. (1996) Randomised quantile residuals,
    \emph{J. Comput. Graph. Statist.}, \bold{5}, 236--244
}    
\author{ 
Natalya Pya <nat.pya@gmail.com> based partly on \code{residual.gam} by Simon N Wood
}
\seealso{   \code{\link{scam}}, \code{\link{scam.check}}
}
\examples{
\donttest{
require(scam)
## Poisson model ...
set.seed(2)
n <- 200
x1 <- runif(n)*6-3
f1 <- 3*exp(-x1^2) # unconstrained term
x2 <- runif(n)*4-1;
f2 <- exp(4*x2)/(1+exp(4*x2)) # monotone increasing smooth
f <- f1+f2
y <- rpois(n,exp(f))
b <- scam(y~s(x1,bs="cr")+s(x2,bs="mpi"),
      family=poisson(link="log"),optimizer=c("efs","bfgs"))
scam.check(b)
## residuals plots using randomized quantile residuals with seed set...
scam.check(b, type="rquantile", setseed=2)
}
}
\keyword{models} \keyword{smooth} \keyword{regression}%












