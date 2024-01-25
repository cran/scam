\name{qq.scam}
\alias{qq.scam}
\title{QQ plots for scam model residuals}
\description{
  Takes a fitted \code{scam} object produced by \code{scam()} and produces
  QQ plots of its residuals (conditional on the fitted model
  coefficients and scale parameter). This is an adapted short version of \code{qq.gam()} of \code{mgcv} package of Simon N Wood.
}
\usage{
qq.scam(object, rep=0, level=.9,s.rep=10,
       type=c("deviance","pearson","response"),
       pch=".", rl.col=3, rep.col="gray80", \dots)
}
\arguments{ 
\item{object}{ a fitted \code{scam} object as produced by \code{scam()} (or a \code{glm} object).}
\item{rep}{How many replicate datasets to generate to simulate quantiles
  of the residual distribution.  \code{0} results in an efficient
  simulation free method for direct calculation, if this is possible for
  the object family.}
\item{level}{If simulation is used for the quantiles, then reference intervals can be provided for the QQ-plot, this specifies the level. 
             0 or less for no intervals, 1 or more to simply plot the QQ plot for each replicate generated.}
\item{s.rep}{how many times to randomize uniform quantiles to data under direct computation.}
\item{type}{what sort of residuals should be plotted?  See
  \code{\link{residuals.scam}}.}
\item{pch}{plot character to use. 19 is good.}
\item{rl.col}{color for the reference line on the plot.}
\item{rep.col}{color for reference bands or replicate reference plots.}
\item{\dots}{extra graphics parameters to pass to plotting functions.}
}
\details{QQ-plots of the the model residuals can be produced in one of two ways. The cheapest method generates reference quantiles by 
associating a quantile of the uniform distribution with each datum, and feeding these uniform quantiles into the quantile function associated with each datum. The resulting quantiles are then used in place of each datum to generate approximate quantiles of residuals.
The residual quantiles are averaged over \code{s.rep} randomizations of the uniform quantiles to data. 


The second method is to use direct simulatation. For each replicate, data are simulated from the fitted model, and the corresponding residuals computed. This is repeated \code{rep} times.
Quantiles are readily obtained from the empirical distribution of residuals so obtained. From this method reference bands are also computable.  

Even if \code{rep} is set to zero, the routine will attempt to simulate quantiles if no quantile function is available for the family. If no random deviate generating function family is available (e.g. for the quasi families), then a normal QQ-plot is produced. The routine conditions on the fitted model coefficents and the scale parameter estimate. 

The plots are very similar to those proposed in Ben and Yohai (2004), but are substantially cheaper to produce (the interpretation of residuals for binary data in Ben and Yohai is not recommended).  

}

\references{

N.H. Augustin, E-A Sauleaub, S.N. Wood (2012) On quantile quantile plots for generalized linear models
Computational Statistics & Data Analysis. 56(8), 2404-2409.

M.G. Ben and V.J. Yohai (2004) JCGS 13(1), 36-47.

}
\author{ Simon N. Wood \email{simon.wood@r-project.org}

Natalya Pya \email{nat.pya@gmail.com} adapted for usage with \code{scam}
}
\seealso{ \code{\link{scam}}}

%\examples{
%}
\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ...




