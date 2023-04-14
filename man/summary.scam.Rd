\name{summary.scam}
%\Rdversion{1.1-6}
\alias{summary.scam}
\alias{print.summary.scam}

\title{Summary for a SCAM fit
}
\description{Takes a fitted \code{scam} object produced by \code{scam()} and produces various useful
summaries from it. The same code as in \code{summary.gam} of the \code{mgcv} package is used here with slight modifications
to accept the exponentiated parameters of the shape constrained smooth terms and the corresponding covariance matrix. 
}
\usage{
\method{summary}{scam}(object, dispersion=NULL, freq=FALSE, ...)

\method{print}{summary.scam}(x,digits = max(3, getOption("digits") - 3), 
                  signif.stars = getOption("show.signif.stars"),...)
}

\arguments{
\item{object}{a fitted \code{scam} object as produced by \code{scam()}.}

\item{x}{a \code{summary.scam} object produced by \code{summary.scam()}.} 

\item{dispersion}{A known dispersion parameter. \code{NULL} to use estimate or
                  default (e.g. 1 for Poisson).}

\item{freq}{By default p-values for individual terms are calculated using the Bayesian estimated 
 covariance matrix of the parameter estimators. If this is set to TRUE then
the frequentist covariance matrix of the parameters is used instead. }

\item{digits}{controls number of digits printed in output.}

\item{signif.stars}{Should significance stars be printed alongside output.}

\item{...}{ other arguments.}

}
% \details{
%%  ~~ If necessary, more details than the description above ~~
% }

\value{\code{summary.scam} produces the same list of summary information for a fitted \code{scam} object as in the 
  unconstrained case \code{summary.gam} except for the last element \code{BFGS termination condition}. 

\item{p.coeff}{is an array of estimates of the strictly parametric model coefficients.}

\item{p.t}{is an array of the \code{p.coeff}'s divided by their standard errors.}

\item{p.pv}{is an array of p-values for the null hypothesis that the corresponding parameter is zero. 
Calculated with reference to the t distribution with the estimated residual
degrees of freedom for the model fit if the dispersion parameter has been
estimated, and the standard normal if not.}

\item{m}{The number of smooth terms in the model.}

\item{chi.sq}{An array of test statistics for assessing the significance of
model smooth terms. See details.}

\item{s.pv}{An array of approximate p-values for the null hypotheses that each
smooth term is zero. Be warned, these are only approximate.}

\item{se}{array of standard error estimates for all parameter estimates.}

\item{r.sq}{The adjusted r-squared for the model. Defined as the proportion of variance explained, where original variance and 
residual variance are both estimated using unbiased estimators. This quantity can be negative if your model is worse than a one 
parameter constant model, and can be higher for the smaller of two nested models! Note that proportion null deviance 
explained is probably more appropriate for non-normal errors.}

\item{dev.expl}{The proportion of the null deviance explained by the model.}

\item{edf}{array of estimated degrees of freedom for the model terms.}

\item{residual.df}{estimated residual degrees of freedom.}

\item{n}{number of data.}

\item{gcv}{minimized GCV score for the model, if GCV used.}

\item{ubre}{minimized UBRE score for the model, if UBRE used.}

\item{scale}{estimated (or given) scale parameter.}

\item{family}{the family used.}

\item{formula}{the original scam formula.}

\item{dispersion}{the scale parameter.}

\item{pTerms.df}{the degrees of freedom associated with each parameteric term
(excluding the constant).}

\item{pTerms.chi.sq}{a Wald statistic for testing the null hypothesis that the
each parametric term is zero.}

\item{pTerms.pv}{p-values associated with the tests that each term is
zero. For penalized fits these are approximate. The reference distribution 
is an appropriate chi-squared when the
scale parameter is known, and is based on an F when it is not.}

\item{cov.unscaled}{The estimated covariance matrix of the parameters (or
estimators if \code{freq=TRUE}), divided
by scale parameter.}

\item{cov.scaled}{The estimated covariance matrix of the parameters
(estimators if \code{freq=TRUE}).}

\item{p.table}{significance table for parameters}

\item{s.table}{significance table for smooths}

\item{p.Terms}{significance table for parametric model terms}
\item{BFGS termination condition}{the value of the maximum component of the scaled GCV/UBRE gradient used as stopping 
     condition. This value is printed if 
       the termination code of the BFGS optimization process is not `1' (not full convergence)
   (see \code{\link{bfgs_gcv.ubre}}for details)}
}

\references{
Wood S.N. (2006) \emph{Generalized Additive Models: An Introduction with R}. Chapman and Hall/CRC Press.

Pya, N. and Wood, S.N. (2015) Shape constrained additive models. Statistics and Computing, 25(3), 543-559

Pya, N. (2010) Additive models with shape constraints. PhD thesis. University of Bath. Department of Mathematical Sciences

}
\author{
Natalya Pya <nat.pya@gmail.com> based on \code{mgcv} by Simon Wood

}

\section{WARNING }{ The p-values are approximate.
} 

\seealso{
 \code{\link{scam}} 
}
\examples{
 \dontrun{
## simulating data...
require(scam)
n <- 200
set.seed(1)
x1 <- runif(n)*6-3
f1 <- 3*exp(-x1^2) # unconstrained smooth term
x2 <- runif(n)*4-1;
f2 <- exp(4*x2)/(1+exp(4*x2)) # monotone increasing smooth
x3 <- runif(n)*5;
f3 <- -log(x3)/5  # monotone decreasing smooth
f <- f1+f2+f3
y <- f + rnorm(n)*.3
dat <- data.frame(x1=x1,x2=x2,x3=x3,y=y)
## fit model ...
b <- scam(y~s(x1,k=15,bs="cr",m=2)+s(x2,k=30,bs="mpi",m=2)+s(x3,k=30,bs="mpd",m=2),
     data=dat)

summary(b)
plot(b,pages=1,shade=TRUE)  
 }  
}
\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ...

