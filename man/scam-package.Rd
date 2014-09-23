\name{scam-package}
%\Rdversion{1.1-8}
\alias{scam-package}
\docType{package}
\title{Shape constrained additive models 
}
\description{
\code{scam} provides functions for generalized additive modelling under shape constraints on the component functions of the linear predictor of the GAM. Models can contain multiple shape constrained and unconstrained terms as well as bivariate smooths with double or single monotonicity. Univariate smooths under  eight possible shape constraints such as monotonically increasing/decreasing, convex/concave, increasing/decreasing and convex, increasing/decreasing and concave, are available as model terms.

The model set up is the same as in \code{mgcv(gam)} with the added shape constrained smooths, so the unconstrained smooths can be of more than one variable, and other 
user defined smooths can be included.  Penalized log likelihood maximization is used to fit the model
together with the automatic smoothness selection.
}
\details{
\tabular{ll}{
Package: \tab scam\cr
Type: \tab Package\cr
Version: \tab 1.1-8\cr
Date: \tab 2014-09-23\cr
License: \tab GPL (version 2 or later) See file LICENSE\cr
LazyLoad: \tab yes\cr
}
The package provides generalized additive modelling under shape constraints on the component functions of the linear predictor. 
\code{\link{scam}} and \code{\link{plot.scam}} functions are based on the functions of the 
unconstrained GAM \code{mgcv(gam)} and \code{mgcv(plot.gam)} and similar in use.  
\code{\link{summary.scam}} allows to extract the results of the model fitting in the same way 
as in \code{summary.gam}. A Bayesian approach is used to obtain a covariance matrix of the model
coefficients and credible intervals for each smooth. 
}
\author{

Natalya Pya <nat.pya@gmail.com> based partly on \code{mgcv} by Simon Wood

Maintainer: Natalya Pya <nat.pya@gmail.com> 

}
\references{

Pya, N. and Wood, S.N. (2014) Shape constrained additive models. Statistics and Computing, DOI \url{http://dx.doi.org/10.1007/s11222-013-9448-7}

Pya, N. (2010) Additive models with shape constraints. PhD thesis. University of Bath. Department of Mathematical Sciences

Wood S.N. (2006) Generalized Additive Models: An Introduction with R. Chapman and Hall/CRC Press.

Wood, S.N. (2008) Fast stable direct fitting and smoothness selection for generalized
additive models. Journal of the Royal Statistical Society (B) 70(3):495-518 

Wood, S.N. (2011) Fast stable restricted maximum likelihood 
and marginal likelihood estimation of semiparametric generalized linear 
models. Journal of the Royal Statistical Society (B) 73(1):3-36

The package was part funded by EPSRC grant EP/I000917/1.

}
\keyword{package}
\keyword{models} 
\keyword{smooth} 
\keyword{regression}
\keyword{monotonicity}
\keyword{convexity}


\examples{
## see examples for scam 
}
