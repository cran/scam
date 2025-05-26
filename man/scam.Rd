\name{scam}
% \Rdversion{1.2-18}
\alias{scam}
\title{Shape constrained additive models (SCAM) and integrated smoothness selection}

\description{This function fits a SCAM to data. Various shape constrained smooths (SCOP-splines), including univariate smooths subject to monotonicity, convexity, or monotonicity plus convexity, bivariate smooths with double or single monotonicity  are available as model terms. See \code{\link{shape.constrained.smooth.terms}} for a complete overview of what is available. Smoothness selection is estimated as part of the fitting. Confidence/credible intervals are available for each smooth term.    

  The shaped constrained smooths have been added to the \code{gam()} in the \code{mgcv} package setup using the \code{smooth.construct} function.  The routine calls a \code{gam()} function for the model set up, but there are separate functions for the model fitting, \code{\link{scam.fit}}, and smoothing parameter selection, \code{\link{bfgs_gcv.ubre}}. Any smooth available in the \code{mgcv} can be taken as a model term for SCAM. User-defined smooths can be included as well.
}

\usage{
scam(formula, family = gaussian(), data = list(), gamma = 1, 
      sp = NULL, weights = NULL, offset = NULL,optimizer=c("bfgs","newton"), 
      optim.method=c("Nelder-Mead","fd"),scale = 0, knots=NULL,
      not.exp=FALSE, start= NULL, etastart=NULL,mustart= NULL,
      control=list(),AR1.rho=0, AR.start=NULL,drop.unused.levels=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{A SCAM formula.  
 This is exactly like the formula for a GAM (see \code{formula.gam} of the \code{mgcv} library) except that shape constrained smooth terms,
 can be added in the expression of the form, e.g., 
 \code{s(x1,k=12,bs="mpi",by=z),}  where \code{bs} indicates the basis to use for the constrained smooth (increasing in this case): 
 the built in options for the shape constrained smooths are described in \cr 
 \code{\link{shape.constrained.smooth.terms}}, 
} 

  \item{family}{A family object specifying the distribution and link to use in
          fitting etc. See \code{\link{glm}} and \code{\link{family}} for more
         details.}

  \item{data}{A data frame or list containing the model response variable and 
      covariates required by the formula. By default the variables are taken 
        from 
       \code{environment(formula)}: typically the environment from 
         which \code{gam} is called.} 

  \item{gamma}{A constant multiplier to inflate the model degrees of freedom in the GCV or UBRE/AIC score.}
            
  \item{sp}{A vector of smoothing parameters can be provided here. Smoothing parameters must be supplied in the order that 
       the smooth terms appear in the model formula. The default \code{sp=NULL} indicates that smoothing parameters should be estimated.
       If \code{length(sp)} does not correspond to the number of underlying smoothing parameters or negative values supplied then
       the vector is ignored and all the smoothing parameters will be estimated.}

  \item{weights}{Prior weights on the data.}

  \item{offset}{Used to supply a model offset for use in fitting. Note that this offset will always be completely ignored when predicting, unlike an offset 
                included in \code{formula}. This conforms to the behaviour of \code{lm}, \code{glm} and \code{gam}.}

  \item{optimizer}{An array specifying the numerical optimization methods 
      to optimize the smoothing parameter estimation criterion (specified in the first 
       element  of \code{optimizer}) and to use to estimate the model coefficients
        (specified in the second element of \code{optimizer}). For the model coefficients
       estimation there are two alternatives: "newton" (default) and "bfgs" methods.
            For the smoothing parameter selection the available methods are "bfgs" (default)
       for the built in to \code{scam} package routine \code{\link{bfgs_gcv.ubre}}, "optim", "nlm", 
           "nlm.fd" (based on finite-difference approximation of the derivatives), "efs". 
     "efs" for the extended Fellner Schall method of Wood and Fasiolo (2017) 
      (rather than minimizing REML as in \code{gam(mgcv)} this minimizes the GCV/UBRE criterion)
     Note that 'bfgs' method for the coefficient estimation works only with 'efs'.}

   \item{optim.method}{In case of \code{optimizer="optim"} this specifies the numerical method to be used in \code{optim}
                      in the first element, the second element of \code{optim.method} indicates whether the finite difference
              approximation should be used ("fd") or analytical gradient ("grad"). The default is \code{optim.method=c("Nelder-Mead","fd")}.}

  \item{scale}{If this is positive then it is taken as the known scale parameter of the exponential family distribution.
     Negative value indicates that the scale paraemter is unknown. 0 indicates that the scale parameter is 1  for Poisson and binomial
    and unknown otherwise. This conforms to the behaviour of \code{gam}.} 
  
  \item{knots}{An optional list containing user specified knot values to be used for basis construction. Different terms can use different numbers of knots.}
  
  \item{not.exp}{if \code{TRUE} then \code{notExp(x,b,threshold)} re-parameterization function will be used in place of \code{exp()} (default value) to ensure positivity in the model coefficients (scop-splines coefficients). \code{notExp()} is a softplus function, \code{1/b*log(1+exp(b*x)}, as implemented in PyTorch, it reverts to the linear function when x*b > threshold, for numerical stability.}

  \item{start}{Initial values for the model coefficients.}

  \item{etastart}{Initial values for the linear predictor.}

  \item{mustart}{Initial values for the expected values.}

  \item{control}{A list of fit control parameters to replace defaults returned by \code{\link{scam.control}}.
                 Values not set assume default values.}

  \item{AR1.rho}{The AR1 correlation parameter. An AR1 error model can be used for the residuals of Gaussian-identity link models. Standardized residuals (approximately uncorrelated under correct model) returned in
\code{std.rsd} if non-zero. }

  \item{AR.start}{logical variable of same length as data, \code{TRUE} at first observation of an independent section of AR1 correlation. Very first observation in data frame does not need this. If \code{NULL} (default) then there are no breaks in AR1 correlaion.}

  \item{drop.unused.levels}{as with \code{gam} by default unused levels are dropped from factors before fitting.}

}

\details{
A shape constrained additive model (SCAM) is a generalized linear model (GLM) 
in which the linear predictor is given by strictly parametric components plus a sum of smooth functions of the covariates where
some of the functions are assumed to be shape constrained. For example,
\deqn{\log(E(Y_i)) = X_i^*b+f_1(x_{1i})+m_2(x_{2i})+f_3(x_{3i})}{log(E(Y_i))=X_i*b+f_1(x_1i)+m_2(x_2i)+f_3(x_3i)}
where the independent response variables \eqn{Y_i}{Y_i} follow Poisson distribution with \code{log} link function,
\eqn{f_1}{f_1}, \eqn{m_2}{m_2}, and \eqn{f_3}{f_3} are smooth functions of the corresponding covariates, and \eqn{m_2}{m_2} 
is subject to monotone increasing constraint.  

Available shape constrained smooths are described in \code{\link{shape.constrained.smooth.terms}}. The shape constrained terms are added to the \code{formula} via \code{s()} terms that can have more that one argument. Unconstrained terms of the \code{mgcv} package, specified using \code{s}, \code{te},
\code{ti}, \code{t2}, can be added as well. There is an option to indicate whether the smooth term is a fixed degrees of freedom regression spline (unpenalized) or a penalized regression spline, by setting the argument \code{fx} of \code{s} term to \code{TRUE} or \code{FALSE} (default). See an example below.

Residual auto-correlation with a simple AR1 correlation structure can be dealt with, for Gaussian models with identity link. Currently, the AR1 correlation parameter should be supplied (rather than estimated) in \code{AR1.rho}. \code{AR.start} input argument (logical) allows to set independent sections of AR1 correlation. Standardized residuals (approximately uncorrelated under correct model) are returned in \code{std.rsd} if \code{AR1.rho} is non zero. Use \code{acf(model$std.rsd)} for computing and plotting estimates of the autocorrelation function to check correlation.

}
\value{The function returns an object of class \code{"scam"} with the following elements (this agrees with \code{gamObject}):

\item{aic}{AIC of the fitted model: the degrees of freedom used to calculate this are the effective degrees of freedom of the model, and the likelihood is evaluated at the maximum of the penalized likelihood, not at the MLE.}

\item{assign}{Array whose elements indicate which model term (listed in
           \code{pterms}) each parameter relates to: applies only to non-smooth terms.}

\item{bfgs.info}{If \code{optimizer[1]="bfgs"}, a list of convergence diagnostics relating to the BFGS method of smoothing parameter selection. The items are: 
              \code{conv}, indicates why the BFGS algorithm of the smoothness selection terminated; \code{iter}, number of iterations of the BFGS taken to get convergence; \code{grad}, the gradient of the GCV/UBRE score at 
              convergence; \code{score.hist}, the succesive values of the score up until convergence.} 

\item{call}{the matched call.}

\item{coefficients}{the coefficients of the fitted model. Parametric
          coefficients are  first, followed  by coefficients for each
          spline term in turn.}

\item{coefficients.t}{the parametrized coefficients of the fitted model (exponentiated for the monotonic smooths).}
\item{conv}{indicates whether or not the iterative fitting method converged.} 
\item{CPU.time}{indicates the real and CPU time (in seconds) taken by the fitting process in case of unknown smoothing parameters}
 
\item{data}{the original supplied data argument. Only included if the \code{scam} argument
    \code{keepData} is set to \code{TRUE} (default is \code{FALSE}).} 

\item{deviance}{model deviance (not penalized deviance).}

\item{df.null}{null degrees of freedom.}

\item{df.residual}{effective residual degrees of freedom of the model.}
           
\item{edf}{estimated degrees of freedom for each model parameter. Penalization
       means that many of these are less than 1.}
\item{edf1}{alternative estimate of edf.}

\item{efs.info}{If \code{optimizer[1]="efs"}, a list of convergence diagnostics relating to the extended Fellner Schall method fot smoothing parameter selection. The items are: 
              \code{conv}, indicates why the efs algorithm of the smoothness selection terminated; \code{iter}, number of iterations of the efs taken to get convergence; \code{score.hist}, the succesive values of the score up until convergence. }

\item{family}{family object specifying distribution and link used.}

\item{fitted.values}{fitted model predictions of expected value for each
          datum.}

\item{formula}{the model formula.}

\item{gcv.ubre}{the minimized GCV or UBRE score.}

\item{dgcv.ubre}{the gradient of the GCV or UBRE score.}

\item{iter}{number of iterations of the Newton-Raphson method taken to get convergence.}

\item{linear.predictors}{fitted model prediction of link function of
expected value for  each datum.}

\item{method}{\code{"GCV"} or \code{"UBRE"}, depending on the fitting criterion used.}

\item{min.edf}{Minimum possible degrees of freedom for whole model.}

\item{model}{model frame containing all variables needed in original model fit.}

\item{nlm.info}{If \code{optimizer[1]="nlm"} or \code{optimizer[1]="nlm.fd"}, a list of convergence diagnostics relating to the BFGS method of smoothing parameter selection. The items are: 
              \code{conv}, indicates why the BFGS algorithm of the smoothness selection terminated; \code{iter}, number of iterations of
              BFGS taken to get convergence; \code{grad}, the gradient of the GCV/UBRE score at 
              convergence.} 

\item{not.exp}{if \code{TRUE} then \code{notExp()} function will be used in place of \code{exp} (default value) in positivity ensuring model parameters re-parameterization.}

\item{nsdf}{number of parametric, non-smooth, model terms including the
          intercept.}

\item{null.deviance}{deviance for single parameter model.} 

\item{offset}{model offset.}

\item{optim.info}{If \code{optimizer[1]="optim"}, a list of convergence diagnostics relating to the BFGS method of smoothing parameter selection. The items are: 
              \code{conv}, indicates why the BFGS algorithm of the smoothness selection terminated; \code{iter}, number of iterations of
              BFGS taken to get convergence; \code{optim.method}, the numerical optimization method used.} 

\item{prior.weights}{prior weights on observations.} 

\item{pterms}{\code{terms} object for strictly parametric part of model.}

\item{R}{Factor R from QR decomposition of weighted model matrix, unpivoted to be in same column order as model matrix.}

\item{residuals}{the working residuals for the fitted model.}

\item{scale.estimated}{\code{TRUE} if the scale parameter was estimated, \code{FALSE} otherwise.}

\item{sig2}{estimated or supplied variance/scale parameter.}

\item{smooth}{list of smooth objects, containing the basis information for each term in the 
model formula in the order in which they appear. These smooth objects are returned by
the \code{\link[mgcv]{smooth.construct}} objects.}

\item{sp}{estimated smoothing parameters for the model. These are the underlying smoothing
parameters, subject to optimization.}

\item{std.rsd}{Standardized residuals (approximately uncorrelated under correct model) if \code{AR1.rho} non zero }

\item{termcode}{an integer indicating why the optimization process of the smoothness selection
          terminated (see \code{\link{bfgs_gcv.ubre}}).}

\item{terms}{\code{terms} object of \code{model} model frame.}

\item{trA}{trace of the influence matrix, total number of the estimated degrees of freedom (\code{sum(edf)}).}

\item{var.summary}{A named list of summary information on the predictor variables. See \code{\link[mgcv]{gamObject}}.}

\item{Ve}{frequentist estimated covariance matrix for the parameter
         estimators.}


\item{Vp}{estimated covariance matrix for the parameters. This is a Bayesian
            posterior covariance matrix that results from adopting a particular Bayesian
            model of the smoothing process.}

\item{Ve.t}{frequentist estimated covariance matrix for the reparametrized parameter
         estimators obtained using the delta method. Particularly useful for testing whether terms are zero. Not so
         useful for CI's as smooths are usually biased.}


\item{Vp.t}{estimated covariance matrix for the reparametrized parameters obtained using the delta method. 
           Paricularly useful for creating credible/confidence intervals.}

\item{weights}{final weights used in the Newton-Raphson iteration.}

%\item{X}{model matrix.}

\item{cmX}{column means of the model matrix (with elements corresponding to smooths set to zero).}

\item{contrasts}{contrasts associated with a factor.}
\item{xlevels}{levels of a factor variable used in the model.}

\item{y}{response data.}


}

\references{
Pya, N. and Wood, S.N. (2015) Shape constrained additive models. Statistics and Computing, 25(3), 543-559

Pya, N. (2010) Additive models with shape constraints. PhD thesis. University of Bath. Department of Mathematical Sciences

Wood, S.N. (2011) Fast stable restricted maximum likelihood 
and marginal likelihood estimation of semiparametric generalized linear 
models. Journal of the Royal Statistical Society (B) 73(1):3-36

Wood S.N. (2017) Generalized Additive Models: An Introduction with R (2nd edition). Chapman
and Hall/CRC Press

Wood, S.N. (2006) On confidence intervals for generalized additive models based on penalized regression splines. Australian and New Zealand Journal of Statistics. 48(4): 445-464.

Wood, S.N. and M. Fasiolo (2017) A generalized Fellner-Schall method for smoothing
parameter optimization with application to Tweedie location, scale and shape models.
Biometrics 73 (4), 1071-1081
}
\author{
Natalya Pya <nat.pya@gmail.com> based partly on \code{mgcv} by Simon Wood

}
\seealso{\code{\link{scam-package}}, \code{\link{shape.constrained.smooth.terms}}, 
    \code{\link[mgcv]{gam}}, \code{\link[mgcv]{s}},
    \code{\link{plot.scam}}, \code{\link{summary.scam}}, \code{\link{scam.check}}, \code{\link{predict.scam}}


}
\examples{
## Gaussian model, two smooth terms: unconstrained and increasing...
## simulating data...
require(scam)
set.seed(4)
n <- 200
x1 <- runif(n)*6-3
f1 <- 3*exp(-x1^2) # unconstrained term
x2 <- runif(n)*4-1;
f2 <- exp(4*x2)/(1+exp(4*x2)) # monotone increasing smooth
y <- f1+f2 +rnorm(n)*.5
dat <- data.frame(x1=x1,x2=x2,y=y)
## fit model, get results, and plot...
b <- scam(y~s(x1,bs="cr")+s(x2,bs="mpi"),data=dat)
b
summary(b)
plot(b,pages=1,shade=TRUE)


## Gaussian model, two smooth terms: increasing and mixed (decreasing and convex)...
## simulating data...
set.seed(4)
n <- 200
x1 <- runif(n)*4-1;
f1 <- exp(4*x1)/(1+exp(4*x1)) # increasing smooth
x2 <- runif(n)*3-1;
f2 <- exp(-3*x2)/15  # decreasing and convex smooth
y <- f1+f2 + rnorm(n)*.4
dat <- data.frame(x1=x1,x2=x2,y=y)
## fit model, results, and plot...
b <- scam(y~ s(x1,bs="mpi")+s(x2, bs="mdcx"),data=dat)
summary(b)
plot(b,pages=1,scale=0,shade=TRUE)


## using the extended Fellner-Schall method for smoothing parameter selection...
b0 <- scam(y~ s(x1,bs="mpi")+s(x2,bs="mdcx"),data=dat,optimizer="efs")
summary(b0)

## using the extended Fellner-Schall method for smoothing parameter selection, 
## and BFGS for model coefficient estimation...
b0 <- scam(y~ s(x1,bs="mpi")+s(x2,bs="mdcx"),data=dat,optimizer=c("efs","bfgs"))
summary(b0)

## using optim() for smoothing parameter selection...
b1 <- scam(y~ s(x1,bs="mpi")+s(x2,bs="mdcx"),data=dat,optimizer="optim")
summary(b1)

b2 <- scam(y~ s(x1,bs="mpi")+s(x2,bs="mdcx"),data=dat,optimizer="optim",
           optim.method=c("BFGS","fd"))
summary(b2)

## using nlm()...
b3 <- scam(y~ s(x1,bs="mpi")+s(x2,bs="mdcx"),data=dat,optimizer="nlm")
summary(b3)


## Poisson model ....
## simulating data...
set.seed(2)
n <- 200
x1 <- runif(n)*6-3
f1 <- 3*exp(-x1^2) # unconstrained term
x2 <- runif(n)*4-1;
f2 <- exp(4*x2)/(1+exp(4*x2)) # monotone increasing smooth
f <- f1+f2
y <- rpois(n,exp(f))
dat <- data.frame(x1=x1,x2=x2,y=y)
## fit model, get results, and plot...
b <- scam(y~s(x1,bs="cr")+s(x2,bs="mpi"),
      family=poisson(link="log"),data=dat,optimizer=c("efs","bfgs"))
summary(b)
plot(b,pages=1,shade=TRUE)
scam.check(b)

## Gamma model...
## simulating data...
set.seed(6)
n <- 300
x1 <- runif(n)*6-3
f1 <- 1.5*sin(1.5*x1) # unconstrained term
x2 <- runif(n)*4-1;
f2 <- 1.5/(1+exp(-10*(x2+.75)))+1.5/(1+exp(-5*(x2-.75))) # increasing smooth
x3 <- runif(n)*6-3;
f3 <- 3*exp(-x3^2)  # unconstrained term
f <- f1+f2+f3
y <- rgamma(n,shape=1,scale=exp(f))
dat <- data.frame(x1=x1,x2=x2,x3=x3,y=y)
## fit model, get results, and plot...
b <- scam(y~s(x1,bs="ps")+s(x2,k=15,bs="mpi")+s(x3,bs="ps"),
          family=Gamma(link="log"),data=dat,optimizer=c("efs","bfgs"))
b
summary(b)
old.par <- par(mfrow=c(2,2))
plot(b,shade=TRUE)
par(old.par) 
  ## run with unpenalized terms...
b1 <- scam(y~s(x1,bs="ps")+s(x2,k=15,bs="mpi",fx=TRUE)+s(x3,bs="ps",fx=TRUE),
          family=Gamma(link="log"),data=dat,optimizer=c("efs","bfgs"))  
summary(b1)


## example with random effect smoother...
set.seed(2)
n <- 200
x1 <- runif(n)*6-3
f1 <- 3*exp(-x1^2) # unconstrained term
x2 <- runif(n)*4-1;
f2 <- exp(4*x2)/(1+exp(4*x2)) # increasing smooth
f <- f1+f2
a <- factor(sample(1:10,200,replace=TRUE))   
Xa <- model.matrix(~a-1)    # random main effects
y <- f + Xa\%*\%rnorm(length(levels(a)))*.5 + rnorm(n)*.4    
dat <- data.frame(x1=x1,x2=x2,y=y,a=a)
## fit model and plot...
b <- scam(y~s(x1,bs="cr")+s(x2,bs="mpi")+s(a,bs="re"), data=dat)
summary(b)
scam.check(b)
plot(b,pages=1,shade=TRUE)

## example with AR1 errors...
set.seed(8)
n <- 500
x1 <- runif(n)*6-3
f1 <- 3*exp(-x1^2) # unconstrained term
x2 <- runif(n)*4-1;
f2 <- exp(4*x2)/(1+exp(4*x2)) # increasing smooth
f <- f1+f2
e <- rnorm(n,0,sd=2)
for (i in 2:n) e[i] <- .6*e[i-1] + e[i]
y <- f + e
dat <- data.frame(x1=x1,x2=x2,y=y)  
b <- scam(y~s(x1,bs="cr")+s(x2,k=25,bs="mpi"),
            data=dat, AR1.rho=.6, optimizer="efs")
b
## raw residuals still show correlation...
acf(residuals(b)) 
## but standardized are now fine...
acf(b$std.rsd)
}

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..
\concept{Varying coefficient model}
\concept{Functional linear model}
\concept{Penalized GLM}
\concept{Generalized Additive Model}
\concept{Penalized regression}
\concept{Spline smoothing}
\concept{Penalized regression spline}
\concept{Generalized Cross Validation}
\concept{Smoothing parameter selection}
\concept{tensor product smoothing}
\concept{P-spline}


