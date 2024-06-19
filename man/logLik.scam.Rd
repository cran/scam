\name{logLik.scam}
\alias{logLik.scam}
\alias{AIC.scam}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Log likelihood for a fitted SCAM, for AIC
}
\description{Function to extract the log-likelihood for a fitted \code{scam}
model (fitted by penalized likelihood maximization). 
Used by \code{AIC}.

 The function is a clone of \code{\link[mgcv]{logLik.gam}} of the \code{mgcv} package. 

The documentation below is similar to that of object \code{\link[mgcv]{logLik.gam}}.
}
\usage{
\method{logLik}{scam}(object,...)
}
%- maybe also `usage' for other objects documented here.
\arguments{ 
\item{object}{ fitted model objects of class \code{scam} as produced by \code{scam()}.}
\item{...}{unused in this case} 

}

\details{ see \code{\link[mgcv]{logLik.gam}} for details.
}

\value{Standard \code{logLik} object: see \code{\link{logLik}}.
}
\references{
Hastie and Tibshirani, 1990, Generalized Additive Models.

Wood, S.N. (2008) Fast stable direct fitting and smoothness selection for
generalized additive models. J.R.Statist. Soc. B 70(3):495-518

}

%\author{ Simon N. Wood \email{simon.wood@r-project.org} }

\seealso{  \code{\link{AIC}}
}

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ...






