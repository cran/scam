\name{vcov.scam}
\alias{vcov.scam}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Parameter estimator covariance matrix from SCAM fit}
\description{ Extracts, from a fitted \code{scam} object, either the Bayesian posterior covariance matrix of the
transformed parameters (a mix of exponentiated and un-exponentiated parameters) or the frequentist covariance matrix of the parameter estimators.
}
\usage{
\method{vcov}{scam}(object, freq = FALSE, untransformed=FALSE, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{ 
\item{object}{fitted model object of class \code{scam} as produced by \code{scam()}.}
\item{freq}{ \code{TRUE} to return the frequentist covariance matrix of the (transformed or untransformed)
parameter estimators, \code{FALSE} to return the Bayesian posterior covariance
matrix of the parameters. } 
\item{untransformed}{ if \code{TRUE} then the covariance matrix of the untransformed parameters is returned. }
\item{...}{ other arguments, currently ignored.}
}
\details{ Extracts, from a fitted \code{scam} object, the Bayesian posterior covariance matrix of the transformed parameters (default; \code{object$Vp.t}), the frequentist covariance matrix of the transformed parameters (\code{object$Ve.t}), or the covariance matrix of the untransformed parameters (\code{object$Vp} for Bayesian and \code{object$Ve} for frequentist inference).
}

\value{ A matrix corresponding to the estimated posterior covariance matrix of the (un)transformed model parameter estimators/coefficients, or the estimated frequentist covariance matrix of the parameters, depending on the arguments \code{freq} and \code{unstransformed}.
}

\author{ Natalya Pya <nat.pya@gmail.com>
}

\references{
Wood, S.N. (2017) Generalized Additive Models: An Introductio with R (2nd ed) CRC Press

Pya, N. and Wood, S.N. (2015) Shape constrained additive models. Statistics and Computing, 25(3), 543-559
}

\seealso{ \code{\link{scam}}}

\examples{ 
require(scam)
set.seed(2)
n <- 200
x <- runif(n)*4-1;
y <- exp(4*x)/(1+exp(4*x)) +rnorm(n)*.2
m <- scam(y~s(x,bs="mpi"))
diag(vcov(m))

}

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ...






