\name{scam-package}
\alias{scam-package}
\docType{package}
\title{
\packageTitle{scam}
}
\description{
\code{scam} provides functions for generalized additive modelling under shape constraints on the component functions of the linear predictor of the GAM. Models can contain multiple univariate and bivariate shape constrained terms, unconstrained terms and parametric terms. A wide variety of shape constrained smooths covered in \code{\link{shape.constrained.smooth.terms}} are provided.

The model set-up is similar to that of \code{gam()} of the package \code{mgcv}, so unconstrained smooths of one or more variables of the \code{mgcv} can be included in SCAMs. User-defined smooths can be added as well. SCAM is estimated by penalized log likelihood maximization and provides automatic smoothness selection by minimizing generalized cross validation or similar. A Bayesian approach is used to obtain a covariance matrix of the model coefficients and credible intervals for each smooth. Linear functionals of smooth functions with shape constraints, parametric model terms, simple linear random effects terms, bivariate interaction smooths with increasing/decreasing constraints (smooth ANOVA), and identity link Gaussian models with AR1 residuals are supported.
}

\details{
  %\packageDESCRIPTION{scam}
  %\packageIndices{scam}
\code{scam} provides generalized additive modelling under shape constraints functions \code{\link{scam}}, \code{\link{summary.scam}}, \code{\link{plot.scam}}, \code{\link{scam.check}}, \code{\link{predict.scam}}, \code{\link{anova.scam}}, and \code{\link{vis.scam}}. These are based on the functions of the unconstrained GAM of the package \code{mgcv} and are similar in use.  

The use of \code{scam()} is much like the use of \code{gam()}, except that within a \code{scam} model formula, shape constrained smooths of one or two predictors can be specified using \code{s} terms with a type of shape constraints used specified as a letter character string of the argument \code{bs}, e.g. \code{s(x, bs="mpi")} for smooth subject to increasing constraint. See \code{\link{shape.constrained.smooth.terms}} for a complete overview of what is available. scam model estimation is performed by penalized likelihood maximization, with smoothness selection by GCV, UBRE/AIC criteria. See \code{\link{scam}}, \code{\link{linear.functional.terms}} for a short discussion of model specification and some examples. See \code{\link{scam}} arguments \code{optimizer} and \code{optim.method}, and \code{\link{scam.control}} for detailed control of scam model fitting. For checking and visualization, see \code{\link{scam.check}}, \code{\link{plot.scam}}, \code{\link{qq.scam}} and \code{\link{vis.scam}}. For extracting fitting results, see \code{\link{summary.scam}} and \code{\link{anova.scam}}.

A Bayesian approach to smooth modelling is used to obtain covariance matrix of the model coefficients and credible intervals for each smooth. \code{Vp} element of the fitted object of class \code{scam} returns the Bayesian covariance matrix, \code{Ve} returns the frequentist estimated covariance matrix for the parameter estimators. The 
frequentist estimated covariance matrix for the reparametrized parameter estimators (obtained using the delta method) is returned in \code{Ve.t}, which is particularly useful for testing individual smooth
terms for equality to the zero function (not so useful for CI's as smooths are usually biased).
\code{Vp.t} returns the Bayesian covariance matrix for the reparametrized parameters. Frequentist approximations can be used for hypothesis testing based on model comparison; see \code{\link{anova.scam}} and \code{\link{summary.scam}} for info on hypothesis testing.

For a complete list of functions type \code{library(help=scam)}.
  
}

\author{
\packageAuthor{scam} based partly on \code{mgcv} by Simon Wood

Maintainer: \packageMaintainer{scam}
}

\references{
Pya, N. and Wood, S.N. (2015) Shape constrained additive models. Statistics and Computing, 25(3), 543-559

Pya, N. (2010) Additive models with shape constraints. PhD thesis. University of Bath. Department of Mathematical Sciences

Wood S.N. (2017) Generalized Additive Models: An Introduction with R (2nd edition). Chapman and Hall/CRC Press

Wood, S.N. (2008) Fast stable direct fitting and smoothness selection for generalized
additive models. Journal of the Royal Statistical Society (B) 70(3):495-518 

Wood, S.N. (2011) Fast stable restricted maximum likelihood 
and marginal likelihood estimation of semiparametric generalized linear 
models. Journal of the Royal Statistical Society (B) 73(1):3-36

The package was part supported by EPSRC grants EP/I000917/1, EP/K005251/1 and the Science Committee of the Ministry of Science and Education of the Republic of Kazakhstan grant #2532/GF3.
}
\keyword{ package }
\keyword{models} 
\keyword{smooth} 
\keyword{regression}
\keyword{monotonicity}
\keyword{convexity}
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
\examples{
## see examples for scam 
}
