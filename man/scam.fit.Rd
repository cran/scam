\name{scam.fit}
%\Rdversion{1.0}
\alias{scam.fit}

\title{Newton-Raphson method to fit SCAM}

\description{
This routine estimates SCAM coefficients given log smoothing parameters using the Newton-Raphson method.
The estimation of the smoothing parameters by the GCV/UBRE score optimization is outer to the model fitting. Routine 
\code{\link{gcv.ubre_grad}} evaluates the first derivatives of the smoothness selection scores with respect to the 
log smoothing parameters. Routine \code{\link{bfgs_gcv.ubre}} estimates the smoothing parameters using the BFGS method. 

The function is not normally called directly, but rather service routines for \code{\link{scam}}. 
}
\usage{
scam.fit(G,sp, SVD=TRUE,ee,eb,esp, maxit=200,epsilon=1e-8,start=NULL,    
          etastart=NULL, mustart=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{G}{A list of items needed to fit a SCAM.}

  \item{sp}{The vector of smoothing parameters.}

  \item{SVD}{Logical, if TRUE (default) then svd is applied to the augmented working model matrix, 
  otherwise the qr decomposition will be used (not recommended).}

  \item{ee}{Get the enviroment for the model coefficients.}

  \item{eb}{Get the enviroment for the model coefficients derivatives.}

  \item{esp}{Get the enviroment for the smoothing parameter.}

   \item{maxit}{Maximum iterations in the Newton-Raphson procedure.}

  \item{epsilon}{A positive scalar giving the tolerance at which the scaled distance between
        two successive penalized deviances is considered close enough to zero to terminate the algorithm.}
  
   \item{start}{Initial values for the model coefficients}

  \item{etastart}{Initial values for the linear predictor}

  \item{mustart}{Initial values for the expected values}
  
}
\details{
The routine applies step halving to any step that increases the
penalized deviance substantially.
}


\references{
Pya, N. (2010) \emph{Additive models with shape constraints}. PhD thesis. University of Bath. Department of Mathematical Sciences

Wood, S.N. (2008) Fast stable direct fitting and smoothness selection for
generalized additive models. J.R.Statist. Soc. B 70(3):495-518

Wood, S.N. (2009) Fast stable REML estimation of semiparametric GLMs. (to be published)
}
\author{
Natalya Pya <nat.pya@gmail.com>

}

\seealso{\code{\link{scam}}}

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ...

