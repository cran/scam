\name{scam}
% \Rdversion{1.1-2}
\alias{scam}
\title{Shape constrained additive models (SCAM) and integrated smoothness selection}

\description{This function fits a SCAM to data. Univariate smooths subject to monotonicity, convexity, or monotonicity plus convexity are available as model terms, as well as
   bivariate smooths with double or single monotonicity. Smoothness selection is estimated as part of the fitting. 
   Confidence/credible intervals are available for each smooth term.    

  All the shaped constrained smooths have been added to the \code{mgcv(gam)} setup using the
\code{smooth.construct} function. The routine calls a \code{mgcv{gam}} function for the model set up, but there are 
  separate functions for the model fitting, \code{\link{scam.fit}}, and smoothing parameter selection, 
  \code{\link{bfgs_gcv.ubre}}. Any unconstrained smooth available in \code{gam} can be taken 
  as model terms. 

 A function \code{\link{extrapolate.uni.scam}}  to predict future values of the response variable in case of a single univariate shape constrained term has been added. 

}
\usage{
scam(formula, family = gaussian(), data = list(), gamma = 1, 
      sp = NULL, weights = NULL, offset = NULL, 
      optimizer="bfgs", optim.method=c("Nelder-Mead","fd"), 
      scale = 0, epsilon = 1e-08, check.analytical=FALSE,
     del=1e-4, start= NULL, etastart, mustart)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{A SCAM formula.  
 This is exactly like the formula for a GAM (see \code{formula.gam} of the \code{mgcv} library) except that monotone smooth terms,
 can be added in the expression of the form \cr
 \code{s(x1,k=12,bs="mpi",by=z),}\cr  where \code{bs} indicates the basis to use for the constrained smooth: 
 the built in options for the monotonic smooths are described in \code{\link{shape.constrained.smooth.terms}}, 
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

  \item{offset}{Can be used to supply a model offset for use in fitting.}

  \item{optimizer}{The numerical optimization method to use to optimize the smoothing 
                   parameter estimation criterion. "bfgs" for the built in to \code{scam} package routine \code{\link{bfgs_gcv.ubre}},
                   "optim", "nlm", "nlm.fd" (based on finite-difference approximation of the derivatives).}

   \item{optim.method}{In case of \code{optimizer="optim"} this specifies the numerical method to be used in \code{optim}
                      in the first element, the second element of \code{optim.method} indicates whether the finite difference
              approximation should be used ("fd") or analytical gradient ("grad"). The default is \code{optim.method=c("Nelder-Mead","fd")}.}

  \item{scale}{If this is positive then it is taken as the known scale parameter of the exponential family distribution.
     Negative value indicates that the scale paraemter is unknown. 0 indicates that the scale parameter is 1  for Poisson and binomial
    and unknown otherwise. This conforms to the behaviour of \code{gam}.} 

  \item{epsilon}{A positive scalar giving the convergence control for the model fitting algorithm.}

  \item{check.analytical}{If this is \code{TRUE} then finite difference derivatives of GCV/UBRE score will be calculated.} 

  \item{del}{A positive scalar (default is 1e-4) giving an increment for finite difference approximation when 
             \code{check.analytical=TRUE}.}

  \item{start}{Initial values for the model coefficients}

  \item{etastart}{Initial values for the linear predictor}

  \item{mustart}{Initial values for the expected values}
}

\details{
A shape constrained additive model (SCAM) is a generalized linear model (GLM) 
in which the linear predictor is given by strictly parametric components plus a sum of smooth functions of the covariates where
some of the functions are assumed to be shape constrained. For example,
\deqn{\log(E(Y_i)) = X_i^*b+f_1(x_{1i})+m_2(x_{2i})+f_3(x_{3i})}{log(E(Y_i))=X_i*b+f_1(x_1i)+m_2(x_2i)+f_3(x_3i)}
where the independent response variables \eqn{Y_i}{Y_i} follow Poisson distribution with \code{log} link function,
\eqn{f_1}{f_1}, \eqn{m_2}{m_2}, and \eqn{f_3}{f_3} are smooth functions of the corresponding covariates, and \eqn{m_2}{m_2} 
is subject to monotone increasing constraint.  

All available shape constrained smooths are decsribed in \code{\link{shape.constrained.smooth.terms}}.


}
\value{The function returns an object of class \code{"scam"} with the following elements (this agrees with \code{gamObject}):

\item{aic}{AIC of the fitted model: the degrees of freedom used to calculate this are the effective degrees of freedom of the model, and the likelihood is evaluated at the maximum of the penalized likelihood, not at the MLE.}

\item{assign}{Array whose elements indicate which model term (listed in
           \code{pterms}) each parameter relates to: applies only to non-smooth terms.}

\item{bfgs.info}{If \code{optimizer="bfgs"}, a list of convergence diagnostics relating to the BFGS method of smoothing parameter selection. The items are: 
              \code{conv}, indicates why the BFGS algorithm of the smoothness selection terminated; \code{iter}, number of iterations of
              BFGS taken to get convergence; \code{grad}, the gradient of the GCV/UBRE score at 
              convergence.} 

\item{optim.info}{If \code{optimizer="optim"}, a list of convergence diagnostics relating to the BFGS method of smoothing parameter selection. The items are: 
              \code{conv}, indicates why the BFGS algorithm of the smoothness selection terminated; \code{iter}, number of iterations of
              BFGS taken to get convergence; \code{optim.method}, the numerical optimization method used.} 

\item{nlm.info}{If \code{optimizer="nlm"} or \code{optimizer="nlm.fd"}, a list of convergence diagnostics relating to the BFGS method of smoothing parameter selection. The items are: 
              \code{conv}, indicates why the BFGS algorithm of the smoothness selection terminated; \code{iter}, number of iterations of
              BFGS taken to get convergence; \code{grad}, the gradient of the GCV/UBRE score at 
              convergence.} 

\item{coefficients}{the coefficients of the fitted model. Parametric
          coefficients are  first, followed  by coefficients for each
          spline term in turn.}

\item{coefficients.t}{the parametrized coefficients of the fitted model (exponentiated for the monotonic smooths).}
\item{conv}{indicates whether or not the iterative fitting method converged.} 
\item{CPU.time}{indicates the real and CPU time (in seconds) taken by the fitting process in case of unknown smoothing parameters}
 
\item{data}{the original supplied data argument.} 
   
\item{deviance}{model deviance (not penalized deviance).}
        
\item{edf}{estimated degrees of freedom for each model parameter. Penalization
       means that many of these are less than 1.}

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

\item{model}{model frame containing all variables needed in original model fit.}

\item{nsdf}{number of parametric, non-smooth, model terms including the
          intercept.}

\item{null.deviance}{deviance for single parameter model.} 

\item{offset}{model offset.}

\item{prior.weights}{prior weights on observations.} 

\item{pterms}{\code{terms} object for strictly parametric part of model.}

\item{residuals}{the working residuals for the fitted model.}

\item{scale.known}{ \code{FALSE} if the scale parameter was estimated, \code{TRUE} otherwise.}

\item{sig2}{estimated or supplied variance/scale parameter.}

\item{smooth}{list of smooth objects, containing the basis information for each term in the 
model formula in the order in which they appear. These smooth objects are returned by
the \code{\link{smooth.construct}} objects.}

\item{sp}{estimated smoothing parameters for the model. These are the underlying smoothing
parameters, subject to optimization.}

\item{termcode}{an integer indicating why the optimization process of the smoothness selection
          terminated (see \code{\link{bfgs_gcv.ubre}}).}

\item{terms}{\code{terms} object of \code{model} model frame.}

\item{trA}{trace of the influence matrix, total number of the estimated degrees of freedom (\code{sum(edf)}).}

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

\item{X}{model matrix.}

\item{cmX}{column means of the model matrix (with elements corresponding to smooths set to zero).}

\item{y}{response data.}


}

\references{
Pya, N. (2010) \emph{Additive models with shape constraints}. PhD thesis. University of Bath. Department of Mathematical Sciences

Wood, S.N. (2008) Fast stable direct fitting and smoothness selection for generalized
additive models. J.R.Statist.Soc.B 70(3):495-518. [Generalized additive model methods]

Wood S.N. (2006a) \emph{Generalized Additive Models: An Introduction with R}. Chapman
and Hall/CRC Press.

Wood, S.N. (2006b) On confidence intervals for generalized additive models based on penalized regression splines. Australian and New Zealand Journal of Statistics. 48(4): 445-464.
}
\author{
Natalya Pya <nat.pya@gmail.com> based partly on \code{mgcv} by Simon Wood

}
\seealso{\code{\link{scam-package}}, \code{\link{shape.constrained.smooth.terms}}, 
    \code{\link{gam}}, \code{\link{s}},
    \code{\link{plot.scam}}, \code{\link{summary.scam}}, \code{\link{scam.check}}, \code{\link{predict.scam}}, \code{\link{extrapolate.uni.scam}}


}


\examples{
##**********************************
## Gaussian model ....
   ## simulating data...

set.seed(2)
n <- 200
x1 <- runif(n)*6-3
f1 <- 3*exp(-x1^2) # unconstrained term
f1 <- (f1-min(f1))/(max(f1)-min(f1)) # function scaled to have range [0,1]
x2 <- runif(n)*4-1;
f2 <- exp(4*x2)/(1+exp(4*x2)) # monotone increasing smooth
f2 <- (f2-min(f2))/(max(f2)-min(f2)) # function scaled to have range [0,1]
f <- f1+f2
y <- f+rnorm(n)*0.1
dat <- data.frame(x1=x1,x2=x2,y=y)
  ## fit model, results, and plot...
b <- scam(y~s(x1,k=15,bs="cr",m=2)+s(x2,k=25,bs="mpi",m=2),
    family=gaussian(link="identity"),data=dat)
print(b)
summary(b)
plot(b,pages=1)


##************************************
## Gaussian model ....
   ## simulating data...

set.seed(2)
n <- 200
x1 <- runif(n)*4-1;
f1 <- exp(4*x1)/(1+exp(4*x1)) # monotone increasing smooth
x2 <- runif(n)*3-1;
f2 <- exp(-3*x2)/15  # monotone decreasing and convex smooth
f <- f1+f2
y <- f+ rnorm(n)*0.2
dat <- data.frame(x1=x1,x2=x2,y=y)
  ## fit model, results, and plot...
b <- scam(y~ s(x1,k=25,bs="mpi",m=2)+s(x2,k=25,bs="mdcx",m=2),
    family=gaussian(link="identity"),data=dat)
print(b)
summary(b)
plot(b,pages=1,scale=0)

##***********************************
\dontrun{
## using optim() for smoothing parameter selection...
b1 <- scam(y~ s(x1,k=25,bs="mpi",m=2)+s(x2,k=25,bs="mdcx",m=2),
    family=gaussian(link="identity"),data=dat, optimizer="optim")
summary(b1)

b2 <- scam(y~ s(x1,k=25,bs="mpi",m=2)+s(x2,k=25,bs="mdcx",m=2),
    family=gaussian(link="identity"),data=dat, optimizer="optim",
    optim.method=c("BFGS","fd"))
summary(b2)

## using nlm()...
b3 <- scam(y~ s(x1,k=25,bs="mpi",m=2)+s(x2,k=25,bs="mdcx",m=2),
    family=gaussian(link="identity"),data=dat, optimizer="nlm")
summary(b3)
}


##************************************
\dontrun{
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
  ## fit model, results, and plot...
b <- scam(y~s(x1,k=15,bs="cr",m=2)+s(x2,k=30,bs="mpi",m=2),
      family=poisson(link="log"),data=dat,optimizer="nlm.fd")
print(b)
summary(b)
plot(b,pages=1)
scam.check(b)

## Gamma model...
   ## simulating data...
set.seed(3)
n <- 200
x1 <- runif(n)*6-3
f1 <- 1.5*sin(1.5*x1) # unconstrained term
x2 <- runif(n)*4-1;
f2 <- 1.5/(1+exp(-10*(x2+0.75)))+1.5/(1+exp(-5*(x2-0.75))) # monotone increasing smooth
x3 <- runif(n)*6-3;
f3 <- 3*exp(-x3^2)  # unconstrained term
f <- f1+f2+f3
y <- rgamma(n,shape=1,scale=exp(f))
dat <- data.frame(x1=x1,x2=x2,x3=x3,y=y)
   ## fit model, results, and plot...
b <- scam(y~s(x1,k=15,bs="ps",m=2)+s(x2,k=30,bs="mpi",m=2)+
            s(x3,k=15,bs="ps",m=2),family=Gamma(link="log"),
            data=dat,optimizer="nlm")
print(b)
summary(b)
par(mfrow=c(2,2))
plot(b)

## bivariate example...
 ## simulating data...
   set.seed(2)
   n <- 30
   x1 <- sort(runif(n)*4-1)
   x2 <- sort(runif(n))
   f1 <- matrix(0,n,n)
   for (i in 1:n) for (j in 1:n) 
       { f1[i,j] <- -exp(4*x1[i])/(1+exp(4*x1[i]))+2*sin(pi*x2[j])}
   f <- as.vector(t(f1))
   y <- f+rnorm(length(f))*0.1
   x11 <-  matrix(0,n,n)
   x11[,1:n] <- x1
   x11 <- as.vector(t(x11))
   x22 <- rep(x2,n)
   dat <- list(x1=x11,x2=x22,y=y)
## fit model  and plot ...
   b <- scam(y~s(x1,x2,k=c(10,10),bs=c("tesmd1","ps"),m=2),
            family=gaussian(link="identity"), data=dat,sp=NULL)
   summary(b)
   par(mfrow=c(2,2),mar=c(4,4,2,2))
   plot(b,se=TRUE)
   plot(b,pers=TRUE,theta = 30, phi = 40)
   plot(y,b$fitted.values,xlab="Simulated data",ylab="Fitted data")

 }
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


