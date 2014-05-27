\name{bfgs_gcv.ubre}
\alias{bfgs_gcv.ubre}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Multiple Smoothing Parameter Estimation by GCV/UBRE}
\description{
Function to efficiently estimate smoothing parameters of SCAM by GCV/UBRE score optimization.
The procedure is outer to the model fitting by the Newton-Raphson method. 
The function uses the BFGS method where the Hessian matrix is updated iteratively at each step.
Backtracking is included to satisfy the sufficient decrease condition.

The function is not normally called directly, but rather service routines for \code{\link{scam}}.
}
\usage{
bfgs_gcv.ubre(fn=gcv.ubre_grad, rho, ini.fd=TRUE, G, gamma=1, ee, eb, 
         esp, n.pen=length(rho), typx=rep(1,n.pen), typf=1, 
        steptol= 1e-7, gradtol = 6.0554*1e-06, maxNstep = 5, 
        maxHalf = 30, check.analytical, del) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fn}{GCV/UBRE Function which returs the GCV/UBRE value and its derivative wrt log smoothing parameter.}

  \item{rho}{log of the initial values of the smoothing parameters.}

  \item{ini.fd}{If TRUE, a finite difference to the Hessian is used to find the initial 
                inverse Hessian, otherwise the initial 
                inverse Hessian is a diagonal matrix `100*I'.}

  \item{G}{A list of items needed to fit a SCAM.}

  \item{gamma}{An ad hoc parametrer of the GCV/UBRE score.}

  \item{ee}{Get the enviroment for the model coefficients.}

  \item{eb}{Get the enviroment for the model coefficients derivatives.}

  \item{esp}{Get the enviroment for the smoothing parameter.}

  \item{n.pen}{Smoothing parameter dimension.}

  \item{typx}{A vector whose component is a positive scalar specifying the typical magnitude of sp.}

  \item{typf}{A positive scalar estimating the magnitude of the gcv near the minimum.}

  \item{steptol}{A positive scalar giving the tolerance at which the scaled distance between
        two successive iterates is considered close enough to zero to terminate the algorithm.}
 
  \item{gradtol}{A scalar giving a tolerance at which the gradient is considered
           to be close enougth to 0 to terminate the algorithm.}

  \item{maxNstep}{A positive scalar which gives the maximum allowable step length.}


  \item{maxHalf}{A positive scalar which gives the maximum number of step halving 
                in "backtracking".}

  \item{check.analytical}{If this is \code{TRUE} then finite difference derivatives of GCV/UBRE score will be calculated.} 

  \item{del}{A positive scalar (default is 1e-4) giving an increment for finite difference approximation when 
             \code{check.analytical=TRUE}.}
 
}
% \details{
%%  ~~ If necessary, more details than the description above ~~
% }
\value{A list is returned with the following items:

  \item{gcv.ubre}{The optimal value of GCV/UBRE.}
  \item{rho}{The best value of the log smoothing parameter.}
  \item{dgcv.ubre}{The gradient of the GCV/UBRE.}
  \item{iterations}{The number of iterations taken until convergence.}
  \item{conv.bfgs}{Convergence information indicating why the BFGS terminated (given below).}
  \item{termcode}{ An integer code indicating why the optimization process
          terminated.

          1: relative gradient is close to zero, current iterate 
              probably is a solution.

          2: scaled distance between last two steps less than `steptol', current iterate 
              probably is a local minimizer, but it's possible that the algorithm is making very slow progress, 
              or `steptol' is too large.

          3: last global step failed to locate a point lower than
              estimate.  Either estimate is an approximate
              local minimum of the function or \code{steptol} is too
              small.

          4: iteration limit exceeded.

          5: five consecutive steps of length \code{maxNstep} have been taken, it's possible that `maxstep' is too small.
              }
  \item{object}{A list of elements returned by the fitting procedure \code{scam.fit} for an optimal value of the smoothing parameter.}
  \item{dgcv.ubre.check}{If \code{check.analytical=TRUE} this is the finite-difference approximation of the gradient calculated by
                            \code{\link{gcv.ubre_grad}}, otherwise \code{NULL}.} 
  \item{check.grad}{If \code{check.analytical=TRUE} this is the relative difference (in %) between the analytical
           and finite differenced derivatives calculated by \code{\link{gcv.ubre_grad}}, otherwise \code{NULL}.} 
}
\references{

Pya, N. and Wood, S.N. (2014) Shape constrained additive models. Statistics and Computing, DOI \url{10.1007/s11222-013-9448-7}

Pya, N. (2010) Additive models with shape constraints. PhD thesis. University of Bath. Department of Mathematical Sciences

Wood, S.N. (2011) Fast stable restricted maximum likelihood and marginal likelihood estimation of semiparametric generalized 
            linear models. Journal of the Royal Statistical Society: Series B. 73(1): 1--34
}
\author{

Natalya Pya <nat.pya@gmail.com>
}


\seealso{
\code{\link{scam}}
}

\keyword{models} 
\keyword{smooth} 
\keyword{regression}



