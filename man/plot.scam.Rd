\name{plot.scam}
%\Rdversion{1.0}
\alias{plot.scam}
\title{SCAM plotting
}
\description{The function is a clone of the \code{plot.gam} of the \code{mgcv} package with the differences
   in the construction of the Bayesian confidence intervals of the shape constrained smooth terms. The function
   takes a fitted \code{scam} object produced by \code{scam()} and plots the 
   component smooth functions that make it up, on the scale of the linear
   predictor. Optionally produces term plots for parametric model components
   as well.

   Note: The fitted shape constrained smooth functions are centred when plotted, which is done in order to be in line with plots of unconstrained smooths (as in gam()). Although 'zeroed intercept' constraints are applied to deal with identifiability of the scop-splines.
}

\usage{
\method{plot}{scam}(x,residuals=FALSE,rug=TRUE,se=TRUE,pages=0,select=NULL,scale=-1,
         n=100,n2=40,pers=FALSE,theta=30,phi=30,jit=FALSE,xlab=NULL,
         ylab=NULL,main=NULL,ylim=NULL,xlim=NULL,too.far=0.1,
         all.terms=FALSE,shade=FALSE,shade.col="gray80",
         shift=0,trans=I,seWithMean=FALSE,unconditional = FALSE, 
         by.resids = FALSE,scheme=0,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{The list of the arguments is the same as in \code{plot.gam} of the \code{mgcv} package.
  \item{x}{ a fitted \code{gam} object as produced by \code{gam()}.}

\item{residuals}{If \code{TRUE} then partial residuals are added to plots of 1-D smooths. If \code{FALSE} 
then no residuals are added. If this is an array of the correct length then it is used as the array of 
residuals to be used for producing partial residuals. If \code{TRUE} then the
residuals are the working residuals from the IRLS iteration weighted by the
IRLS weights. Partial residuals for a smooth term are the
residuals that would be obtained by dropping the term concerned from the model, while leaving all other 
estimates fixed (i.e. the estimates for the term plus the residuals).}

\item{rug}{ when TRUE (default) then the covariate to which the plot applies is displayed as a rug plot
 at the foot of each plot of a 1-d smooth, and the locations of the
 covariates are plotted as points on the contour plot representing a 2-d
 smooth.} 
  
\item{se}{ when TRUE (default) upper and lower lines are added to the
  1-d plots at 2 standard errors
        above and below the estimate of the smooth being plotted while for
	2-d plots, surfaces at +1 and -1 standard errors are contoured
	and overlayed on the contour plot for the estimate. If a
	positive number is supplied then this number is multiplied by
	the standard errors when calculating standard error curves or
  surfaces. See also \code{shade}, below. }

\item{pages}{ (default 0) the number of pages over which to spread the output. For example, 
if \code{pages=1} then all terms will be plotted on one page with the layout performed automatically. 
Set to 0 to have the routine leave all graphics settings as they are. }

\item{select}{Allows the  plot for a single model term to be selected for printing. e.g. if you just want the plot for the second smooth term set \code{select=2}. }

\item{scale}{ set to -1 (default) to have the same y-axis scale for each plot, and to 0 for a 
          different y axis for each plot. Ignored if \code{ylim} supplied.}

\item{n}{ number of points used for each 1-d plot - for a nice smooth plot this needs to be several times the estimated 
  degrees of freedom for the smooth. Default value 100.}

\item{n2}{Square root of number of points used to grid estimates of 2-d
  functions for contouring.}

\item{pers}{Set to \code{TRUE} if you want perspective plots for 2-d
  terms.}

\item{theta}{One of the perspective plot angles.}

\item{phi}{The other perspective plot angle.}

\item{jit}{Set to TRUE if you want rug plots for 1-d terms to be jittered.}

\item{xlab}{If supplied then this will be used as the x label for all plots.}

\item{ylab}{If supplied then this will be used as the y label for all plots.}

\item{main}{Used as title (or z axis label) for plots if supplied.}

\item{ylim}{If supplied then this pair of numbers are used as the y limits for each plot.}

\item{xlim}{If supplied then this pair of numbers are used as the x limits for each plot.}

\item{too.far}{If greater than 0 then this is used to determine when a location is too
far from data to be plotted when plotting 2-D smooths. This is useful since smooths tend to go wild away from data.
The data are scaled into the unit square before deciding what to exclude, and \code{too.far} is a distance 
within the unit square.}

\item{all.terms}{if set to \code{TRUE} then the partial effects of parametric
model components are also plotted, via a call to \code{\link{termplot}}. Only
terms of order 1 can be plotted in this way.}

\item{shade}{Set to \code{TRUE} to produce shaded regions as confidence bands
for smooths (not avaliable for parametric terms, which are plotted using \code{termplot}).}

\item{shade.col}{define the color used for shading confidence bands.}

\item{shift}{constant to add to each smooth (on the scale of the linear
predictor) before plotting. Can be useful for some diagnostics, or with \code{trans}.}

\item{trans}{function to apply to each smooth  (after any shift), before
plotting. \code{shift} and \code{trans} are occasionally useful as a means for
getting plots on the response scale, when the model consists only of a single smooth.}

\item{seWithMean}{if \code{TRUE} the component smooths are shown with confidence 
intervals that include the uncertainty about the overall mean. If \code{FALSE} then the 
uncertainty relates purely to the centred smooth itself. An extension of the argument
presented in Nychka (1988) suggests that \code{TRUE} results in better coverage performance,
and this is also suggested by simulation.}

\item{unconditional}{ if \code{TRUE} then the smoothing parameter uncertainty
          corrected covariance matrix is used to compute uncertainty
          bands, if available. Otherwise the bands treat the smoothing
          parameters as fixed.}

\item{by.resids}{Should partial residuals be plotted for terms with \code{by} variables? 
Usually the answer is no, they would be meaningless.}

\item{scheme}{Integer (0,1 or 2) or integer vector selecting a plotting scheme for each plot. 
  \code{scheme == 0} produces a smooth curve with dashed curves 
indicating 2 standard error bounds. \code{scheme == 1} illustrates the error bounds using a shaded
region. For \code{scheme==0}, contour plots are produced for 2-d smooths with the x-axes labelled with the first covariate
name and the y axis with the second covariate name. For 2-d smooths \code{scheme==1} produces a perspective plot, while \code{scheme==2} produces a heatmap, with overlaid contours. }

\item{...}{ other graphics parameters to pass on to plotting commands.}

}

% \details{ to be added...
%%  ~~ If necessary, more details than the description above ~~
% }
\value{The function generates plots.

}
\references{

Pya, N. and Wood, S.N. (2015) Shape constrained additive models. Statistics and Computing, 25(3), 543-559

Pya, N. (2010) Additive models with shape constraints. PhD thesis. University of Bath. Department of Mathematical Sciences

Wood S.N. (2006) Generalized Additive Models: An Introduction with R. Chapman and Hall/CRC Press.


}
\author{
Natalya Pya <nat.pya@gmail.com> based on the \code{plot.gam} of the \code{mgcv} by Simon Wood

}

\seealso{  \code{\link{scam}}}

\examples{
## simulating data...
n <- 200
set.seed(1)
x0 <- rep(1:4,50)
x1 <- runif(n)*6-3
f1 <- 3*exp(-x1^2) # unconstrained smooth term
x2 <- runif(n)*4-1;
f2 <- exp(4*x2)/(1+exp(4*x2)) # monotone increasing smooth
x3 <- runif(n)*5;
f3 <- -log(x3)/5  # monotone decreasing smooth
f <- f1+f2+f3
y <- 2*x0 + f + rnorm(n)*0.3
x0 <- factor(x0)

## fit the model and plot ...
b <- scam(y~x0+s(x1,k=15,bs="cr")+s(x2,k=30,bs="mpi")+s(x3,k=30,bs="mpd"))
plot(b,pages=1,residuals=TRUE,all.terms=TRUE,shade=TRUE,shade.col=3)    


\dontrun{
 ## example with 2-d plots...
 ## simulating data...
   set.seed(2)
   n <- 30
   x0 <- rep(1:9,100)
   x1 <- sort(runif(n)*4-1)
   x2 <- sort(runif(n))
   x3 <- runif(n*n, 0, 1)
   f <- matrix(0,n,n)
   for (i in 1:n) for (j in 1:n) 
       { f[i,j] <- -exp(4*x1[i])/(1+exp(4*x1[i]))+2*sin(pi*x2[j])}
   f1 <- as.vector(t(f))
   f2 <- x3*0
   e <- rnorm(length(f1))*0.1
   y <- 2*x0 + f1 + f2 + e
   x0 <- factor(x0)
   x11 <-  matrix(0,n,n)
   x11[,1:n] <- x1
   x11 <- as.vector(t(x11))
   x22 <- rep(x2,n)
   dat <- list(x0=x0,x1=x11,x2=x22,x3=x3,y=y)
## fit model  and plot ...
   b <- scam(y~x0+s(x1,x2,k=c(10,10),bs=c("tesmd1","ps"),m=2)+s(x3),data=dat)
   op <- par(mfrow=c(2,2))
   plot(b,all.terms=TRUE)
   plot(y,b$fitted.values,xlab="Simulated data",ylab="Fitted data",pch=19,cex=.3)
   par(op) 
   
## and use of schemes...
   op <- par(mfrow=c(2,2))
   plot(b,all.terms=TRUE,scheme=1)
   par(op)
   op <- par(mfrow=c(2,2))
   plot(b,all.terms=TRUE,scheme=c(2,1))
   par(op)

  }
}

\keyword{models} \keyword{smooth} \keyword{regression} \keyword{hplot} %-- one or more ...

