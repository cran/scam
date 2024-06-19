\name{residuals.scam}
\alias{residuals.scam}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{SCAM residuals}
\description{This function is a clone of the \code{mgcv} library code \code{\link[mgcv]{residuals.gam}}.
It returns residuals for a fitted \code{scam} model
  object. Pearson, deviance, working and response residuals are
  available. 
}

\usage{
\method{residuals}{scam}(object, type = c("deviance", "pearson","scaled.pearson", 
                        "working", "response"),...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{object}{ a \code{scam} fitted model object. }
  \item{type}{the type of residuals wanted. }
   \item{...}{other arguments.}
}
\details{See  \code{\link[mgcv]{residuals.gam}} for details.
}
\value{ An array of residuals.
}

\author{ 
Natalya Pya <nat.pya@gmail.com> 
}

\seealso{   \code{\link{scam}}}

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..












