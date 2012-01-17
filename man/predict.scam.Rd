\name{predict.scam}
\alias{predict.scam}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Prediction from fitted SCAM model}
\description{ 
This function is a clone of the \code{mgcv} library code \code{\link{predict.gam}} with some modifications 
to adopt shape preserving smooth terms.
It takes a fitted \code{scam} object produced by \code{scam()} 
and produces predictions given a new set of values for the model covariates 
or the original values used for the model fit. Predictions can be accompanied
by standard errors, based on the posterior distribution of the model
coefficients. }

\usage{
\method{predict}{scam}(object,newdata,type="link",se.fit=FALSE,terms=NULL,
    block.size=1000,newdata.guaranteed=FALSE,na.action=na.pass,...)
}
%- maybe also `usage' for other objects documented here.
\arguments{ 
  
 \item{object}{ a fitted \code{scam} object as produced by \code{scam()}.
                }
 \item{newdata}{ A data frame or list containing the values of the model covariates at which predictions
  are required. If this is not provided then predictions corresponding to the
                original data are returned. If \code{newdata} is provided then
                it should contain all the variables needed for prediction: a
                warning is generated if not. } 
 
\item{type}{ When this has the value \code{"link"} (default) the linear predictor (possibly with
associated standard errors) is returned. When \code{type="terms"} each component of the 
linear predictor is returned seperately (possibly with standard errors): this includes 
parametric model components, followed by each smooth component, but excludes
any offset and any intercept. \code{type="iterms"} is the same, except that any standard errors 
returned for smooth components will include the uncertainty about the intercept/overall mean.  When 
\code{type="response"} predictions 
on the scale of the response are returned (possibly with approximate
standard errors). When \code{type="lpmatrix"} then a matrix is returned
which yields the values of the linear predictor (minus any offset) when
postmultiplied by the
parameter vector (in this case \code{se.fit} is ignored). The latter
option is most useful for getting variance estimates for quantities derived from
the model: for example integrated quantities, or derivatives of smooths. A
linear predictor matrix can also be used to implement approximate prediction
outside \code{R} (see example code, below). }

\item{se.fit}{ when this is TRUE (not default) standard error estimates are returned for each prediction.}

\item{terms}{if \code{type=="terms"} then only results for the terms given in this array
will be returned.}

\item{block.size}{maximum number of predictions to process per call to underlying
code: larger is quicker, but more memory intensive. Set to < 1 to use total number
of predictions as this.}

\item{newdata.guaranteed}{Set to \code{TRUE} to turn off all checking of
\code{newdata} except for sanity of factor levels: this can speed things up
for large prediction tasks, but \code{newdata} must be complete, with no
\code{NA} values for predictors required in the model. }

\item{na.action}{what to do about \code{NA} values in \code{newdata}. With the
default \code{na.pass}, any row of \code{newdata} containing \code{NA} values
for required predictors, gives rise to \code{NA} predictions (even if the term concerned has no
\code{NA} predictors). \code{na.exclude} or \code{na.omit} result in the
dropping of \code{newdata} rows, if they contain any \code{NA} values for
required predictors. If \code{newdata} is missing then \code{NA} handling is 
determined from \code{object$na.action}.}

\item{...}{ other arguments.}

}

\value{ If \code{type=="lpmatrix"} then a matrix is returned which will
  give a vector of linear predictor values (minus any offest) at the supplied covariate
  values, when applied to the model coefficient vector. 
Otherwise, if \code{se.fit} is \code{TRUE} then a 2 item list is returned with items (both arrays) \code{fit}
and \code{se.fit} containing predictions and associated standard error estimates, otherwise an 
array of predictions is returned. The dimensions of the returned arrays depends on whether 
\code{type} is \code{"terms"} or not: if it is then the array is 2 dimensional with each 
term in the linear predictor separate, otherwise the array is 1 dimensional and contains the 
linear predictor/predicted values (or corresponding s.e.s). The linear predictor returned termwise will 
not include the offset or the intercept.

\code{newdata} can be a data frame, list or model.frame: if it's a model frame
then all variables must be supplied.

}

\details{See \code{\link{predict.gam}} for details.
}

\references{

Chambers and Hastie (1993) Statistical Models in S. Chapman & Hall.

Wood S.N. (2006) \emph{Generalized Additive Models: An Introduction with R}. Chapman and Hall/CRC Press.

Pya, N. (2010) \emph{Additive models with shape constraints}. PhD thesis. University of Bath. Department of Mathematical Sciences

}
\author{ 
Natalya Pya <nat.pya@gmail.com> based partly on \code{mgcv} by Simon Wood

}

\seealso{  \code{\link{scam}}, \code{\link{plot.scam}}}

\examples{
\dontrun{
library(scam)
set.seed(2)
n <- 200
x1 <- runif(n)*6-3
f1 <- 3*exp(-x1^2) # unconstrained term
x2 <- runif(n)*4-1;
f2 <- exp(4*x2)/(1+exp(4*x2)) # monotone increasing smooth
f <- f1+f2
y <- f+rnorm(n)*0.2
dat <- data.frame(x1=x1,x2=x2,y=y)
b <- scam(y~s(x1,k=15,bs="cr",m=2)+s(x2,k=30,bs="mpi",m=2),
    family=gaussian(link="identity"),data=dat)

newd <- data.frame(x1=seq(-2.9,2.9,length.out=30),x2=seq(-0.9,2.9,length.out=30))
pred <- predict.scam(b,newd)
pred
predict(b,newd,type="terms",se=TRUE)
plot(b,se=TRUE,residuals=TRUE,pages=1)
 }
}

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..

