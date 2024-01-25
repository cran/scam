\name{shape.constrained.smooth.terms}
\alias{shape.constrained.smooth.terms}
\title{Shape preserving smooth terms in SCAM}
\description{
As in \code{mgcv(gam)}, shape preserving smooth terms are specified in a \code{\link{scam}} formula using \code{\link{s}} 
terms. All the shape constrained smooth terms (SCOP-splines) are constructed using the B-splines basis proposed by Eilers and Marx (1996) with a discrete penalty on the basis coefficients.

The univariate single penalty built-in shape constrained smooth classes are summarized as follows.
\itemize{
\item Monotone increasing SCOP-splines: \code{bs="mpi"}. To achieve monotone increasing smooths this reparameterizes the coefficients so that they form an increasing sequence. 

For details see \code{\link{smooth.construct.mpi.smooth.spec}}. 

\item Monotone decreasing SCOP-splines: \code{bs="mpd"}. To achieve monotone decreasing smooths this reparameterizes the coefficients 
so that they form a decreasing sequence. A first order difference penalty applied to the basis coefficients starting with the 
second is used for the monotone increasing and decreasing cases.

\item Convex SCOP-splines: \code{bs="cx"}. This reparameterizes the coefficients so that the second order differences of the basis coefficients are greater than zero. 

For details see \code{\link{smooth.construct.cx.smooth.spec}}. 

\item Concave SCOP-splines: \code{bs="cv"}. This reparameterizes the coefficients so that the second order differences of the basis coefficients are less than zero.

 For details see \code{\link{smooth.construct.cv.smooth.spec}}. 

\item Increasing and convex SCOP-splines: \code{bs="micx"}. This reparameterizes the coefficients 
so that the first and the second order differences of the basis coefficients are greater than zero. 
For details see 
\code{\link{smooth.construct.micx.smooth.spec}}.

\item Increasing and concave SCOP-splines: \code{bs="micv"}. This reparameterizes the coefficients 
so that the first order differences of the basis coefficients are greater than zero while the second order difference
are less than zero.

\item Decreasing and convex SCOP-splines: \code{bs="mdcx"}. This reparameterizes the coefficients 
so that the first order differences of the basis coefficients are less than zero while the second order difference
are greater. 
For details see \code{\link{smooth.construct.mdcx.smooth.spec}}.

\item Decreasing and concave SCOP-splines: \code{bs="mdcv"}. This reparameterizes the coefficients 
so that the first and the second order differences of the basis coefficients are less than zero.

\item Increasing with an additional 'finish-at-zero' constraint SCOP-splines: \code{bs="mifo"}. This sets the last \code{(m+1)} spline coefficients to zero. According to the B-spline basis functions properties, the value of the spline, \code{f(x)}, is determined by \code{m+2} non-zero basis functions, and only \code{m+1} B-splines are non-zero at knots. Only \code{m+2} B-splines are non-zero on any \code{(k_i, k_{i+1})}, and the sum of these \code{m+2} basis functions is 1. 

For details see \code{\link{smooth.construct.mifo.smooth.spec}}. 

\item Increasing with an additional 'start-at-zero' constraint SCOP-spline: \code{bs="miso"}. This sets the first \code{(m+1)} spline coefficients to zero. According to the B-spline basis functions properties, the value of the spline, \code{f(x)}, is determined by \code{m+2} non-zero basis functions, and only \code{m+1} B-splines are non-zero at knots. Only \code{m+2} B-splines are non-zero on any \code{(k_i, k_{i+1})}, and the sum of these \code{m+2} basis functions is 1. 

For details see \code{\link{smooth.construct.miso.smooth.spec}}. 

\item SCOP-spline with positivity constraint: \code{bs="po"}. This reparameterizes the coefficients so that there are positive. For details see \code{\link{smooth.construct.po.smooth.spec}}. 

\item Decreasing/increasing SCOP-splines used with numeric 'by' variable: \code{bs="mpdBy"}, \code{bs="mpiBy"}. These work similar to \code{mpd.smooth.spec}, \code{mpi.smooth.spec},  but without applying an identifiability constraint ('zero intercept' constraint). Use when the smooth term has a numeric \code{by} variable  that takes more than one value.
 For details see \code{\link{smooth.construct.mpd.smooth.spec}}, \code{\link{smooth.construct.mpi.smooth.spec}}. 

\item Convex/concave SCOP-splines used with numeric 'by' variable: \code{bs="cxdBy"}, \code{bs="cvBy"}. These work similar to \code{cx.smooth.spec}, \code{cv.smooth.spec},  but without applying an identifiability constraint ('zero intercept' constraint). Use when the smooth term has a numeric \code{by} variable  that takes more than one value.

 For details see \code{\link{smooth.construct.cx.smooth.spec}}, 

\code{\link{smooth.construct.cv.smooth.spec}}. 

\item Decreasing/increasing and convex/concave SCOP-splines used with numeric 'by' variable: \code{bs="mdcxBy"}, \code{bs="mdcvBy"}, \code{bs="micxBy"}, \code{bs="micvBy"}. 

These work similar to \code{mdcx.smooth.spec}, \code{mdcv.smooth.spec}, \code{micx.smooth.spec}, \code{micv.smooth.spec}, but without applying an identifiability constraint ('zero intercept' constraint). Use when the smooth term has a numeric \code{by} variable  that takes more than one value.
 For details see \code{\link{smooth.construct.mdcx.smooth.spec}}, 

\code{\link{smooth.construct.mdcv.smooth.spec}},

\code{\link{smooth.construct.micx.smooth.spec}}, 

\code{\link{smooth.construct.micv.smooth.spec}}. 

}

For all types of the mixed constrained smoothing a first order difference penalty applied to the basis coefficients starting with the third one is used. 


Using the concept of the tensor product spline bases bivariate smooths under monotonicity constraint where monotonicity may be 
assumed on only one of the covariates (single monotonicity) or both of them (double monotonicity) are added as the smooth terms 
of the SCAM. Bivariate B-spline is constructed by expressing the coefficients of one of the marginal univariate 
B-spline bases as the B-spline of the other covariate. Double or single monotonicity is achieved by the corresponding 
re-parametrization of the bivariate basis coefficients to satisfy the sufficient conditions formulated in terms of the first 
order differences of the coefficients. The following explains the built in bivariate shape constrained smooth classes.
\itemize{
\item Double monotone increasing SCOP-splines: \code{bs="tedmi"}. 

   See \code{\link{smooth.construct.tedmi.smooth.spec}} for details.

\item Double monotone decreasing SCOP-splines: \code{bs="tedmd"}. 

\item Single monotone increasing SCOP-splines along the first covariate direction: \code{bs="tesmi1"}. 

\item Single monotone increasing SCOP-splines along the second covariate direction: \code{bs="tesmi2"}. 
    
\item Single monotone decreasing SCOP-splines along the first covariate direction: \code{bs="tesmd1"}. 
     
\item Single monotone decreasing SCOP-splines along the second covariate direction: \code{bs="tesmd2"}. 

\item SCOP-splines with double concavity constraint: \code{bs="tecvcv"}.

 See \code{\link{smooth.construct.tecvcv.smooth.spec}} for details.

\item SCOP-splines with double convexity constraint: \code{bs="tecxcx"}. 

See \code{\link{smooth.construct.tecxcx.smooth.spec}} for details.

\item SCOP-splines with convexity wrt the first covariate and
concavity wrt the second covariate: \code{bs="tecxcv"}. See \code{\link{smooth.construct.tecxcv.smooth.spec}} for details.

\item Decreasing along the first covariate and
concave along the second covariate SCOP-splines: \code{bs="tedecv"}. 
See \code{\link{smooth.construct.tedecv.smooth.spec}} for details.

\item Decreasing along the first covariate and convex along the second covariate SCOP-splines: \code{bs="tedecx"}. 
See \code{\link{smooth.construct.tedecx.smooth.spec}} for details.

\item Increasing along the first covariate and concave along the second covariate SCOP-splines: \code{bs="temicv"}. 
See \code{\link{smooth.construct.temicv.smooth.spec}} for details.

\item Increasing along the first covariate and convex along the second covariate SCOP-splines: \code{bs="temicx"}. 
See \code{\link{smooth.construct.temicx.smooth.spec}} for details.

\item Convex along the second covariate SCOP-splines: \code{bs="tescx"}. 

See \code{\link{smooth.construct.tescx.smooth.spec}} for details.      

\item Concave along the second covariate SCOP-splines: \code{bs="tescv"}. 

See \code{\link{smooth.construct.tescv.smooth.spec}} for details.

\item Tensor product interaction with increasing constraint along the first
 covariate and unconstrained along the second covariate: \code{bs="tismi"}. 

See \code{\link{smooth.construct.tismi.smooth.spec}} for details.

\item Tensor product interaction with decreasing constraint along the first
 covariate and unconstrained along the second covariate: \code{bs="tismd"}. 

See \code{\link{smooth.construct.tismd.smooth.spec}} for details.
 }
Double penalties for the shape constrained tensor product smooths are obtained from the penalties of the marginal smooths.

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
\code{\link{smooth.construct.tesmd2.smooth.spec}},

\code{\link{smooth.construct.tismi.smooth.spec}},
\code{\link{smooth.construct.tismd.smooth.spec}}

} 

\author{
Natalya Pya \email{nat.pya@gmail.com}
}
\references{
Pya, N. and Wood, S.N. (2015) Shape constrained additive models. Statistics and Computing, 25(3), 543-559

Pya, N. (2010) Additive models with shape constraints. PhD thesis. University of Bath. Department of Mathematical Sciences

Eilers, P.H.C. and B.D. Marx (1996) Flexible Smoothing with B-splines and Penalties. 
Statistical Science, 11(2):89-121

Wood S.N. (2017) Generalized Additive Models: An Introduction with R (2nd edition). Chapman
and Hall/CRC Press

Wood, S.N. (2006) Low rank scale invariant tensor product smooths for
generalized additive mixed models. Biometrics 62(4):1025-1036

}

\examples{
## see examples for scam 
}
\keyword{regression}
