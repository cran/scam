\name{scam.control}
\alias{scam.control}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Setting SCAM fitting defaults}
\description{ This is an internal function of package \code{scam} which allows 
control of the numerical options for fitting a SCAM. 
}
\usage{
scam.control(maxit = 200, maxHalf=30, devtol.fit=1e-7, steptol.fit=1e-7,
            keepData=FALSE,efs.lspmax=15,efs.tol=.1, nlm=list(),optim=list(),
            bfgs=list(), trace =FALSE, print.warn=FALSE) 
}
\arguments{ 

\item{maxit}{Maximum number of IRLS iterations to perform used in \code{\link{scam.fit}}.}

\item{maxHalf}{If a step of  the BFGS optimization method leads 
to a worse penalized deviance, then the step length of the model coefficients is halved. This is
the number of halvings to try before giving up used in \code{\link{bfgs_gcv.ubre}}.}

\item{devtol.fit}{A positive scalar giving the convergence control for the model fitting algorithm in \code{\link{scam.fit}}.}

\item{steptol.fit}{A positive scalar giving the tolerance at which the scaled distance between
        two successive iterates is considered close enough to zero to terminate the model fitting algorithm in \code{\link{scam.fit}}.}

\item{keepData}{Should a copy of the original \code{data} argument be kept in the \code{scam} 
object?}

\item{efs.lspmax}{maximum log smoothing parameters to allow under extended Fellner Schall
smoothing parameter optimization.}

\item{efs.tol}{change in GCV to count as negligible when testing for EFS convergence. If the
step is small and the last 3 steps led to a GCV change smaller than this, then stop.}

\item{nlm}{list of control parameters to pass to \code{\link{nlm}} if this is
used for outer estimation of smoothing parameters (not default). }

\item{optim}{list of control parameters to pass to \code{\link{optim}} if this
is used for outer estimation of smoothing parameters (not default).}

\item{bfgs}{list of control parameters to pass to default BFGS optimizer
used for outer estimation of log smoothing parameters.}

\item{trace}{turns on or off some de-bugging information.}

\item{print.warn}{ when set to \code{FALSE} turns off printing warning messages for step halving under non-finite exponentiated coefficients,  non-finite deviance and/or if \code{mu} or \code{eta} are out of bounds.}
}

\details{ Outer iteration is used to estimate smoothing parameters of SCAM by GCV/UBRE score optimization. The default procedure is the built-in BFGS method which is controlled by the list \code{bfgs} with the following elements: \code{steptol.bfgs} (default 1e-7) is the relative convergence tolerance; 
\code{gradtol.bfgs} (default 6.0554*1e-6) is a tolerance at which the gradient is considered to be close enougth to 0 to terminate the BFGS algorithm; 
\code{maxNstep} is a positive scalar which gives the maximum allowable step length (default 5);
\code{maxHalf} gives the maximum number of step halving in "backtracking" to permit before giving up(default 30);
\code{check.analytical} is logical whether the analytical gradient of GCV/UBRE should be checked numerically (default \code{FALSE});
\code{del} is an increment for finite differences when checking analytical gradients (default 1e-4). 

If outer iteration using \code{nlm} is used for fitting, then the control list \code{nlm} stores control arguments
for calls to routine \code{nlm}. As in \code{gam.control} the list has the following named elements: \code{ndigit} is the number
of significant digits in the GCV/UBRE score; \code{gradtol} is the tolerance used to judge convergence of the gradient of the GCV/UBRE score to zero 
(default 1e-6); \code{stepmax} is the maximum allowable log smoothing parameter step (default 2);
 \code{steptol} is the minimum allowable step length (default 1e-4); 
\code{iterlim} is the maximum number of optimization steps allowed (default 200);
\code{check.analyticals} indicates whether the built in exact derivative calculations should be checked numerically (default \code{FALSE}).
Any of these which are not supplied and named in the list are set to their default values.

Outer iteration using \code{optim} is controlled using list \code{optim}, which currently has one element: \code{factr}
which takes default value 1e7.


}


\references{
Pya, N. and Wood, S.N. (2015) Shape constrained additive models. Statistics and Computing, 25(3), 543-559

Pya, N. (2010) Additive models with shape constraints. PhD thesis. University of Bath. Department of Mathematical Sciences

Wood, S.N. (2011) Fast stable restricted maximum likelihood 
and marginal likelihood estimation of semiparametric generalized linear 
models. Journal of the Royal Statistical Society (B) 73(1):3-36

}
\author{ Natalya Pya Arnqvist \email{nat.pya@gmail.com}  based partly on \code{\link{gam.control}} by Simon Wood}


\seealso{\code{\link{scam}}, \code{\link{scam.fit}}, \code{\link{gam.control}} }

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..




