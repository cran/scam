\name{anova.scam}
\alias{anova.scam}
\alias{print.anova.scam}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Approximate hypothesis tests related to  SCAM fits}
\description{Performs hypothesis tests relating to one or more fitted \code{scam} objects.
 The function is a clone of \code{anova.gam} of the \code{mgcv} package. 

The documentation below is similar to that of object \code{\link{anova.gam}}.
}
\usage{
\method{anova}{scam}(object, ..., dispersion = NULL, test = NULL,
                    freq = FALSE,p.type=0)
\method{print}{anova.scam}(x, digits = max(3, getOption("digits") - 3),...)
}
%- maybe also `usage' for other objects documented here.
\arguments{ 

\item{object,...}{ fitted model objects of class \code{scam} as produced by \code{scam()}.}
\item{x}{an \code{anova.scam} object produced by a single model call to \code{anova.scam()}.} 
\item{dispersion}{ a value for the dispersion parameter: not normally used.}
\item{test}{what sort of test to perform for a multi-model call. One of
\code{"Chisq"}, \code{"F"} or \code{"Cp"}. }
\item{freq}{whether to use frequentist or Bayesian approximations for parametric term 
p-values. See \code{\link{summary.gam}} for details.}
\item{p.type}{selects exact test statistic to use for single smooth term p-values. See
 \code{\link{summary.scam}} for details.}
\item{digits}{number of digits to use when printing output.}
}

\details{ see \code{\link{anova.gam}} for details.
}

\value{In the multi-model case \code{anova.scam} produces output identical to
\code{\link{anova.glm}}, which it in fact uses.

In the single model case an object of class \code{anova.scam} is produced,
which is in fact an object returned from \code{\link{summary.scam}}.

\code{print.anova.scam} simply produces tabulated output.

}

\references{

Scheipl, F., Greven, S. and Kuchenhoff, H. (2008) Size and power of tests for a zero random effect variance or polynomial 
regression in additive and linear mixed models. Comp. Statist. Data Anal. 52, 3283-3299

Wood, S.N. (2013a) On p-values for smooth components of an extended generalized additive model. Biometrika 100:221-228

Wood, S.N. (2013b) A simple test for random effects in regression models. Biometrika 100:1005-1010
}

\author{ Simon N. Wood \email{simon.wood@r-project.org} }

\section{WARNING}{ If models 'a' and 'b' differ only in terms with no un-penalized components then 
p values from anova(a,b) are unreliable, and usually much too low.

Default P-values will usually be wrong for parametric terms penalized using `paraPen': use freq=TRUE
to obtain better p-values when the penalties are full rank and represent conventional random effects.

For a single model, interpretation is similar to drop1, not anova.lm.
} 

\seealso{  \code{\link{scam}}, \code{\link{predict.scam}},
\code{\link{scam.check}}, \code{\link{summary.scam}}, \code{\link{anova.gam}} 
}

\examples{
library(scam)
set.seed(0)
fac <- rep(1:4,20)
x1 <- runif(80)*5
x2 <- runif(80,-1,2)
x3 <- runif(80, 0, 1)
y <- fac+log(x1)/5
y <- y + exp(-1.3*x2) +rnorm(80)*0.1
fac <- factor(fac)
b <- scam(y ~ fac+s(x1,bs="mpi")+s(x2,bs="mpd")+s(x3))

b1 <- scam(y ~ fac+s(x1,bs="mpi")+s(x2,bs="mpd"))
anova(b,b1,test="F")

## b2 <- scam(y ~ fac +s(x1)+s(x2)+te(x1,x2))

}
\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ...






