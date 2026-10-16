\name{model.matrix.scam}
\alias{model.matrix.scam}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Extract model matrix from SCAM fit}
\description{Obtains the model matrix from a fitted \code{scam} object.
}
\usage{
\method{model.matrix}{scam}(object, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{ 
\item{object}{ fitted model object of class \code{scam} as produced by \code{scam()}.}

\item{...}{ other arguments, passed to \code{\link{predict.scam}}.}
}
\details{Calls \code{\link{predict.scam}} with no \code{newdata} argument and
\code{type="lpmatrix"} in order to obtain the model matrix of \code{object}.
}

\value{ A model matrix.
}

\seealso{ \code{\link{scam}}}

\examples{ 
require(scam)
set.seed(4)
n <- 200
x <- runif(n)*4-1
y <- exp(4*x)/(1+exp(4*x)) +rnorm(n)*.5
m <- scam(y ~ s(x,bs="mpi"))
X <- model.matrix(m)
}

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ...






