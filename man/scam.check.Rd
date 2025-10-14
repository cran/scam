\name{scam.check}
%\Rdversion{1.0}
\alias{scam.check}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Some diagnostics for a fitted scam object
}
\description{
The function takes a fitted \code{scam} object produced by \code{scam()} and provides diagnostic information about the model fitting procedure and results. It is analogous to the \code{gam.check()} function from the mgcv package. By default, the function produces three deviance residual plots and a plot of the response versus fitted values, and prints information about the convergence of the smoothness selection optimization. Alternatively, the function can generate plots for checking the normalized (randomized) quantile residuals of the fitted \code{scam} object.
}
\usage{
scam.check(b,type=c("deviance","rquantile","pearson","response"),old.style=FALSE, 
        pch=".", setseed=NULL, rep=0, level=.9, rl.col=3, rep.col="gray80",\dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{b}{ a fitted \code{scam} object as produced by \code{scam()}.}
 \item{old.style}{produces qq-norm plots as it was in scam versions < 1.2-15 when set to \code{TRUE}.}
  \item{type}{type of residuals, see \code{\link{residuals.scam}}, used in
    all plots.}
  \item{rep, level, rep.col}{
    arguments passed to \code{\link{qq.scam}()} when \code{old.style} is
    \code{FALSE} (default).} 
  \item{rl.col}{color for the reference line on the quantile-quantile plot.}
  \item{pch}{plot character to use for the quantile-quantile plot.}
  \item{setseed}{seed to be used in getting randomized quantile residuals.}
  \item{\dots}{extra graphics parameters to pass to plotting functions.}

}
\details{\code{scam.check} produces four standard residual plots, and some convergence diagnostics. By default, the residuals are deviance residuals. The QQ plot in this case is created by a call to \code{qq.scam()}, and plots the deviance residuals against approximate theoretical quantilies of the deviance residual distribution, according to the fitted model.

The type of the residuals can be changed to \code{"rquantile"}, (randomized) quantile residuals as described by Dunn and Smyth (1996). For models with non-normal response distribution, (randomized) quantile residuals may provide a better diagnostic tool. The QQ plot for these residuals is a standard normal QQ plot.
Randomization is applied only when the response variable is discrete.

The printed output includes information about the optimization process used for smoothing parameter selection.
}

\references{
Wood S.N. (2006) Generalized Additive Models: An Introduction with R. Chapman and Hall/CRC Press.

Dunn, P. K. and Smyth, G. K. (1996) Randomised quantile residuals,
    \emph{J. Comput. Graph. Statist.}, \bold{5}, 236--244
}
\author{

Natalya Pya \email{nat.pya@gmail.com} based partly on \code{mgcv} by Simon N Wood
}
\seealso{  \code{\link{residuals.scam}}
}

\examples{
library(scam)
set.seed(2)
n <- 200
x1 <- runif(n)*4-1; x2 <- runif(n)*3-1;
f1 <- function(x) exp(4*x)/(1+exp(4*x)) # monotone increasing smooth
f2 <- function(x) exp(-3*x)/15  # monotone decreasing and convex smooth
f <- f1(x1) + f2(x2)
y <- f+ rnorm(n)*0.2
dat <- data.frame(x1=x1,x2=x2,y=y)
b <- scam(y~ s(x1,bs="mpi")+s(x2,bs="mdcx"),data=dat)
plot(b,pages=1)
scam.check(b)
 
 \donttest{
## binomial model... 
set.seed(4)
n <- 400
scale <- 3.3
x0 <- runif(n, 0, 1)
f0 <- function(x) 2 * sin(pi * x)
x1 <- runif(n)*4-1
f1 <- function(x) exp(4*x)/(1+exp(4*x)) 
f <- f0(x0) + f1(x1) 
f <- (f-2)*scale
g <- binomial()$linkinv(f)
y <- rbinom(g,1,g)
dat <- data.frame(y=y,x0=x0,x1=x1)

b1 <- scam(y~s(x0)+s(x1,bs="mpi"),family=binomial,data=dat)
scam.check(b1)
## residuals plots using randomized quantile residuals...
scam.check(b1, type="rquantile")

## Residuals vs linear predictor for four realizations of the quantile residuas...
old.par <- par(mfrow=c(2,2))
plot(b1$linear.predictors,residuals(b1, type="rquantile"))
plot(b1$linear.predictors,residuals(b1, type="rquantile"))
plot(b1$linear.predictors,residuals(b1, type="rquantile"))
plot(b1$linear.predictors,residuals(b1, type="rquantile"))
par(old.par)
}
}

\keyword{models} 
\keyword{smooth} 
\keyword{regression} %-- one or more ..

