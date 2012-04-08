\name{shape.constrained.smooth.terms}
\alias{shape.constrained.smooth.terms}
\title{Shape preserving smooth terms in SCAM}
\description{
As in \code{mgcv(gam)}, shape preserving smooth terms are specified in a \code{\link{scam}} formula using \code{\link{s}} 
terms. All the shape constrained smooth terms are constructed using the B-splines basis proposed by Eilers and Marx (1996) with 
a discrete penalty on the basis coefficients.

The univariate single penalty built in shape constrained smooth classes are summarized as follows
\itemize{
\item{Monotone increasing P-splines}{ \code{bs="mpi"}. To achieve monotone increasing smooths these reparametrize the coefficients 
so that they form an increasing sequence. 
For details see 
\code{\link{smooth.construct.mpi.smooth.spec}}.} 

\item{Monotone decreasing P-splines}{ \code{bs="mpd"}. To achieve monotone decreasing smooths these reparametrize the coefficients 
so that they form a decreasing sequence. A first order difference penalty applied to the basis coefficients starting with the 
second is used for the monotone increasing and decreasing cases.}

\item{Convex P-splines}{ \code{bs="cx"}. These reparametrize the coefficients so that the second order differences of the basis coefficients are greater than zero. For details see 
\code{\link{smooth.construct.cx.smooth.spec}}.} 

\item{Concave P-splines}{ \code{bs="cv"}. These reparametrize the coefficients so that the second order differences of the basis coefficients are less than zero. For details see 
\code{\link{smooth.construct.cv.smooth.spec}}.} 

\item{Monotone increasing and convex P-splines}{ \code{bs="micx"}. These reparametrize the coefficients 
so that the first and the second order differences of the basis coefficients are greater than zero. 
For details see 
\code{\link{smooth.construct.micx.smooth.spec}}.} 

\item{Monotone increasing and concave P-splines}{ \code{bs="micv"}. These reparametrize the coefficients 
so that the first order differences of the basis coefficients are greater than zero while the second order difference
are less than zero.} 

\item{Monotone decreasing and convex P-splines}{ \code{bs="mdcx"}. These reparametrize the coefficients 
so that the first order differences of the basis coefficients are less than zero while the second order difference
are greater. 
For details see \code{\link{smooth.construct.mdcx.smooth.spec}}.} 

\item{Monotone decreasing and concave P-splines}{ \code{bs="mdcv"}. These reparametrize the coefficients 
so that the first and the second order differences of the basis coefficients are less than zero.}

For all four types of the mixed constrained smoothing a first order difference 
penalty applied to the basis coefficients starting with the third one is used. 
}

Using the concept of the tensor product spline bases bivariate smooths under monotonicity constraint where monotonicity may be 
assumed on only one of the covariates (single monotonicity) or both of them (double monotonicity) are added as the smooth terms 
of the SCAM. Bivariate B-spline is constructed by expressing the coefficients of one of the marginal univariate 
B-spline bases as the B-spline of the other covariate. Double or single monotonicity is achieved by the corresponding 
re-parametrization of the bivariate basis coefficients to satisfy the sufficient conditions formulated in terms of the first 
order differences of the coefficients. The following explains the built in bivariate monotonic smooth classes.
\itemize{
\item{Double monotone increasing P-splines}{ \code{bs="tedmi"}. 
   See \code{\link{smooth.construct.tedmi.smooth.spec}} for details.} 

\item{Double monotone decreasing P-splines}{ \code{bs="tedmd"}. 
    }

\item{Single monotone increasing P-splines along the first covariate direction}{ \code{bs="tesmi1"}. 
     }

\item{Single monotone increasing P-splines along the second covariate direction}{ \code{bs="tesmi2"}. 
     }

\item{Single monotone decreasing P-splines along the first covariate direction}{ \code{bs="tesmd1"}. 
     }

\item{Single monotone decreasing P-splines along the second covariate direction}{ \code{bs="tesmd2"}. 
     }

Double penalties for the monotonic tensor product smooths are obtained from the penalties of the marginal smooths.
}

}

\seealso{\code{\link{s}}, \code{\link{smooth.construct.mpi.smooth.spec}}, 
\code{\link{smooth.construct.mpd.smooth.spec}},
\code{\link{smooth.construct.cx.smooth.spec}}, 
\code{\link{smooth.construct.cv.smooth.spec}},
\code{\link{smooth.construct.micx.smooth.spec}},
\code{\link{smooth.construct.micv.smooth.spec}},
\code{\link{smooth.construct.mdcx.smooth.spec}}, 
\code{\link{smooth.construct.mdcv.smooth.spec}}, 
\code{\link{smooth.construct.tedmi.smooth.spec}},
\code{\link{smooth.construct.tedmd.smooth.spec}},
\code{\link{smooth.construct.tesmi1.smooth.spec}}, 
\code{\link{smooth.construct.tesmi2.smooth.spec}}, 
\code{\link{smooth.construct.tesmd1.smooth.spec}}, 
\code{\link{smooth.construct.tesmd2.smooth.spec}}} 

\author{
Natalya Pya <nat.pya@gmail.com>
}
\references{
Pya, N. (2010) \emph{Additive models with shape constraints}. PhD thesis. University of Bath. Department of Mathematical Sciences

Eilers, P.H.C. and B.D. Marx (1996) Flexible Smoothing with B-splines and Penalties. 
Statistical Science, 11(2):89-121

Wood S.N. (2006a) \emph{Generalized Additive Models: An Introduction with R}. Chapman and Hall/CRC Press.

Wood, S.N. (2006b) Low rank scale invariant tensor product smooths for
generalized additive mixed models. Biometrics 62(4):1025-1036

}

\examples{
## see examples for scam 
}
\keyword{regression}
