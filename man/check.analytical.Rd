\name{check.analytical}
%\Rdversion{1.0}
\alias{check.analytical}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Checking the analytical gradient of the GCV/UBRE score
}
\description{This function calculates the finite-difference approximation of the GCV/UBRE gradient for the fitted 
         model and compares it with the analytical gradient. 
}
\usage{
check.analytical(object, data, del=1e-6,control)
}
\arguments{
 
  \item{object}{A fitted \code{scam} object.}

  \item{data}{An original data frame or list containing the model response variable and covariates.}

  \item{del}{A positive scalar (default is 1e-6) giving an increment for finite difference approximation.} 

   \item{control}{Control option list as returned by \code{\link{scam.control}}.}
}


\value{A list is returned with the following items:

  \item{dgcv.ubre.fd}{The finite-difference approximation of the gradient.}

  \item{check.grad}{The relative difference in percentage between the analytical
           and finite differenced derivatives.}

}
\author{

Natalya Pya <nat.pya@gmail.com>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{scam}}
}


\keyword{models} 
\keyword{smooth} 
\keyword{regression}

