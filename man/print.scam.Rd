\name{print.scam}
%\Rdversion{1.0}
\alias{print.scam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Print a SCAM object
}
\description{The default print method for a \code{scam} object. The code is a clone of \code{print.gam}
of the \code{mgcv} package with a slight simplification since only two methods of smoothing parameter 
selection (by GCV or UBRE) was implemented for \code{scam}.
}
\usage{
\method{print}{scam}(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{ fitted model objects of class \code{scam} as produced by \code{scam()}.}

\item{...}{ other arguments.}

}
\details{As for \code{mgcv(gam)} prints out the family, model formula, effective degrees of freedom for each smooth term, 
   and optimized value of the smoothness selection criterion used.
}



\references{
Wood S.N. (2006) Generalized Additive Models: An Introduction with R. Chapman and Hall/CRC Press.

}
\author{

Natalya Pya <nat.pya@gmail.com> 
}
\seealso{   \code{\link{scam}}, \code{\link{summary.scam}}}

\keyword{models} \keyword{smooth} \keyword{regression} %-- one or more ..

