\name{formula.scam}
\alias{formula.scam}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{SCAM formula}
\description{ Description of \code{\link{scam}} formula (see \code{gam} of the \code{mgcv} package for Details), and how to extract it from a fitted \code{scam} object.

The function is a clone of \code{\link[mgcv]{formula.gam}} of the \code{mgcv} package. 
}
\usage{
\method{formula}{scam}(x,...)
}
%- maybe also `usage' for other objects documented here.
\arguments{ 
\item{x}{ fitted model objects of class \code{scam} as produced by \code{scam()}.}
\item{...}{un-used in this case} 
}

\value{ Returns the model formula, \code{x$formula}. Provided so that \code{anova} methods
print an appropriate description of the model.
}

\details{ see \code{\link[mgcv]{formula.gam}} for details.
}


%\author{ Simon N. Wood \email{simon.wood@r-project.org}}


\seealso{  \code{\link{scam}}}


\keyword{models} \keyword{smooth} \keyword{regression} \keyword{internal} %-- one or more ...






