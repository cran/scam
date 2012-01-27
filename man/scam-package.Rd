\name{scam-package}
%\Rdversion{1.1}
\alias{scam-package}
\docType{package}
\title{Shape constrained additive models 
}
\description{
\code{scam} provides functions for generalized additive modelling under shape constraints on the component functions of the linear predictor of the GAM. Models can contain multiple shape constrained and unconstrained terms as well as bivariate smooths with double or single monotonicity. The model set up is the same as in \code{mgcv(gam)} with the added 
shape constrained smooths, so the unconstrained smooths can be of more than one variable, and other 
user defined smooths can be included.  Penalized log likelihood maximization is used to fit the model
together with the automatic smoothness selection.
}
\details{
\tabular{ll}{
Package: \tab scam\cr
Type: \tab Package\cr
Version: \tab 1.1-1\cr
Date: \tab 2012-27-01\cr
License: \tab GPL (version 2 or later) See file LICENSE\cr
LazyLoad: \tab yes\cr
}
The package provides generalized additive modelling under shape constraints on the component functions of the linear predictor. 
\code{\link{scam}} and \code{\link{plot.scam}}. These functions are based on the functions of the 
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
Pya, N. (2010) \emph{Additive models with shape constraints}. PhD thesis. University of Bath. Department of Mathematical Sciences

Wood S.N. (2006) \emph{Generalized Additive Models: An Introduction with R}. Chapman and Hall/CRC Press.

Wood, S.N. (2008) Fast stable direct fitting and smoothness selection for generalized
additive models. J.R.Statist.Soc.B 70(3):495-518 

Wood, S.N. (2011) Fast stable restricted maximum likelihood and marginal likelihood estimation of semiparametric generalized 
            linear models. \emph{Journal of the Royal Statistical Society: Series B}. 73(1): 1--34

}
\keyword{ package }
\keyword{models} 
\keyword{smooth} 
\keyword{regression}
\keyword{monotonicity}


\examples{
## see examples for scam 
}
