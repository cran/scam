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
scam.fit(G, sp, maxit=200, devtol=1e-8, steptol=1e-8, gamma=1, start=NULL,    
          etastart=NULL, mustart=NULL, env = env)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{G}{A list of items needed to fit a SCAM.}

  \item{sp}{The vector of smoothing parameters.}

  \item{env}{Get the enviroment for the model coefficients, their derivatives and the smoothing parameter.}

  \item{maxit}{Maximum iterations in the Newton-Raphson procedure.}

  \item{devtol}{A positive scalar giving the tolerance at which the scaled distance between
        two successive penalized deviances is considered close enough to zero to terminate the algorithm.}
  
   \item{steptol}{A positive scalar giving the tolerance at which the scaled distance between
        two successive iterates is considered close enough to zero to terminate the algorithm.}

   \item{gamma}{This constant allows to inflate the model degrees of
          freedom in the GCV or UBRE/AIC score.}  

   \item{start}{Initial values for the model coefficients}

  \item{etastart}{Initial values for the linear predictor}

  \item{mustart}{Initial values for the expected values}

    
}
\details{
The routine applies step halving to any step that increases the
penalized deviance substantially.
}


\references{
Pya, N. and Wood, S.N. (2015) Shape constrained additive models. Statistics and Computing, 25(3), 543-559

Pya, N. (2010) Additive models with shape constraints. PhD thesis. University of Bath. Department of Mathematical Sciences

Wood, S.N. (2008) Fast stable direct fitting and smoothness selection for
generalized additive models. Journal of the Royal Statistical Society (B) 70(3):495-518

Wood, S.N. (2011) Fast stable restricted maximum likelihood 
and marginal likelihood estimation of semiparametric generalized linear 
models. Journal of the Royal Statistical Society (B) 73(1):3-36

}
\author{
Natalya Pya <nat.pya@gmail.com>

}

\seealso{\code{\link{scam}}}

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ...

