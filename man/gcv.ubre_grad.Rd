\name{gcv.ubre_grad}
%\Rdversion{1.0}
\alias{gcv.ubre_grad}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{The GCV/UBRE score value and its gradient
}
\description{
For the estimation of the SCAM smoothing parameters the GCV/UBRE score is optimized outer to the Newton-Raphson 
procedure of the model fitting. 
This function returns the value of the GCV/UBRE score  and calculates its first derivative 
with respect to the log smoothing parameter using the method of Wood (2009).

The function is not normally called directly, but rather service routines for \code{\link{bfgs_gcv.ubre}}.
}
\usage{
gcv.ubre_grad(rho, G, env, control)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 
  \item{rho}{log of the initial values of the smoothing parameters.}

  \item{G}{a list of items needed to fit a SCAM.}

%  \item{gamma}{A constant multiplier to inflate the model degrees of freedom in the GCV or UBRE/AIC score.}

  \item{env}{Get the enviroment for the model coefficients, their derivatives and the smoothing parameter.}
  
  \item{control}{A list of fit control parameters as returned by \code{scam.control}.}

}


\value{A list is returned with the following items:

  \item{dgcv.ubre}{The value of GCV/UBRE gradient.}

  \item{gcv.ubre}{The GCV/UBRE score value.}

  \item{scale.est}{The value of the scale estimate.}

  \item{object}{The elements of the fitting procedure \code{monogam.fit} for a given value of the smoothing parameter.}

  \item{dgcv.ubre.check}{If \code{check.analytical=TRUE} this returns the finite-difference approximation of the gradient.}

  \item{check.grad}{If \code{check.analytical=TRUE} this returns the relative difference (in %) between the analytical
           and finite differenced derivatives.}

}
\references{

Pya, N. and Wood, S.N. (2015) Shape constrained additive models. Statistics and Computing, 25(3), 543-559

Pya, N. (2010) Additive models with shape constraints. PhD thesis. University of Bath. Department of Mathematical Sciences

Wood S.N. (2006) Generalized Additive Models: An Introduction with R. Chapman and Hall/CRC Press.

Wood, S.N. (2011) Fast stable restricted maximum likelihood and marginal likelihood estimation of semiparametric generalized 
            linear models. Journal of the Royal Statistical Society: Series B. 73(1): 1--34
}
\author{

Natalya Pya <nat.pya@gmail.com>

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{scam}}, \code{\link{scam.fit}}, \code{\link{bfgs_gcv.ubre}}
}


\keyword{models} 
\keyword{smooth} 
\keyword{regression}
