## (c) Natalya Pya (2012-2025). Provided under GPL 2.
## similar to residuals.gam() (c) Simon N Wood....
## (2025) added quantile residuals in version 1.2-20...

residuals.scam <-function(object, type = c("deviance", "pearson","scaled.pearson", "working", "response", "rquantile"),setseed=NULL, ...)
# calculates residuals for scam object the same as residulas.gam(), but from version 1.2-20 includes (randomized) quantile residuals...
{ type <- match.arg(type)
  y <- object$y
  mu <- object$fitted.values
  y.mu <- y-mu
##  family <- object$family
  wts <- object$prior.weights
  res<- switch(type, working = object$residuals,
                     scaled.pearson = y.mu*wts^.5/sqrt(object$sig2*object$family$variance(mu)),
                     pearson = y.mu*wts^.5/sqrt(object$family$variance(mu)),
                     deviance = { d.res <- sqrt(pmax(object$family$dev.resids(y,mu,wts),0))
                                  ifelse(y>mu , d.res, -d.res)             
                                 },
                     response = y.mu,
                     rquantile = rquantile.residuals(y=y,fv=mu, wt=wts, scale=object$sig2, family=object$family,setseed=setseed)
              )
  res <- naresid(object$na.action,res)
  res
}


## based on residuals routines (c) Gavin L. Simpson r package gratia and (c) Mikis Stasinopoulos et al r package gamlss...

rquantile.residuals <- function(y, fv, wt, scale, family,setseed=NULL)
## calculates non-randomized quantile residuals for continuous distributions 
## and randomized quantile residuals for discrete distributions (Dunn and Smyth, 1996; Feng et al. 2020)...
## y - observed response, fv - fitted values, wt - prior weights, scale - dispersion or scale
## setseed allows to set seeds
{
  if (is.null(family$cdf)) {
    family <- add.family.cdf(family)
  }
  # if it's still NULL...
  if (is.null(family$cdf)) {
    stop("Quantile residuals are not available for this family.")
  }
  cdf <- family$cdf(q=y, mu=fv, wt=wt, scale=scale, log.p=FALSE)
  fam.name <- family[1]
  discrete <- FALSE
  if (fam.name == "binomial" || fam.name == "poisson") 
           discrete <- TRUE 
  if(!discrete) ## for continuous distributions get non-randomized quantile residuals
     rquantile <- qnorm(cdf)
  else { ## for discrete distributions -> randomized quantile residuals
      aval <- family$cdf(q = y - 1L, mu = fv, wt = wt, scale = scale, log.p = FALSE) 
      # generate a random value between aval and cdf...
      if (!is.null(setseed)) set.seed(setseed)  
      uval <- runif(n = length(y), aval, cdf)  
      uval <- ifelse(uval>0.999999,uval-.1e-15,uval) 
      uval <- ifelse(uval<0.000001,uval+ .1e-15,uval)
      rquantile <- qnorm(uval)  
  }  
  rquantile  
}


## adding cdf to the families, needed for calculating the residual quantiles...
## based on cdf routines (c) Gavin L. Simpson, r package gratia...

add.family.cdf <- function(family) 
## adding cdf to family
{  if (!is.null(family$cdf)) {
    return(family)
  }
  fam.name <- family[1]$family
  qfun <- switch(EXPR = fam.name,
             poisson = cdf.poisson,
             binomial = cdf.binom,
             gaussian = cdf.norm,
             Gamma  = cdf.gamma,
             inverse.gaussian = cdf.invgaussian
           )
  # add the CDF fun to the family
  family$cdf <- qfun
  family
}


cdf.poisson <- function(q, mu, wt, scale, log.p = FALSE) {
  ppois(q, lambda = mu, log.p = log.p)
}

cdf.norm <- function(q, mu, wt, scale, log.p = FALSE) {
  pnorm(q, mean = mu, sd = sqrt(scale / wt), log.p = log.p)
}

cdf.binom <- function(q, mu, wt, scale, log.p = FALSE) {
  pbinom(q * (wt + as.numeric(wt == 0)), size = wt, prob = mu, log.p = log.p)
}

cdf.gamma <- function(q, mu, wt, scale, log.p = FALSE) {
   pgamma(q, shape = 1 / scale, scale = mu * scale, log.p = log.p)
}

## inverse gaussian taken from (c) Mikis Stasinopoulos r package gamlss.dist....

cdf.invgaussian <- function(q, mu = 1, wt, scale = 1,  log.p = FALSE) 
## lower.tail = TRUE, only for lower tail of CDF
{    sigma <- scale
    if (any(mu < 0))  stop(paste("mu must be positive", "\n", "")) 
    if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
  ## if (any(q < 0))  stop(paste("y must be positive", "\n", ""))  
      lq <- length(q)                                                                    
   sigma <- rep(sigma, length = lq)
      mu <- rep(mu, length = lq)           
    cdf1 <- pnorm(((q/mu)-1)/(sigma*sqrt(q))) 
   lcdf2 <- (2/(mu*sigma^2))+pnorm((-((q/mu)+1))/(sigma*sqrt(q)),log.p=TRUE)
     cdf <- cdf1+ exp(lcdf2)
  ##  if(lower.tail==TRUE) cdf  <- cdf else  cdf <- 1-cdf 
    if(log.p==FALSE) cdf  <- cdf else  cdf <- log(cdf) 
    cdf <- ifelse(q <= 0, 0, cdf) 
    cdf
}




##########################################################
## residual.scam() function used before version 1.2-20 ...

residuals.scam_19 <-function(object, type = c("deviance", "pearson","scaled.pearson", "working", "response"),...)
# calculates residuals for scam object the same as residulas.gam()...
{ type <- match.arg(type)
  y <- object$y
  mu <- object$fitted.values
  y.mu <- y-mu
##  family <- object$family
  wts <- object$prior.weights
  res<- switch(type,working = object$residuals,
         scaled.pearson = y.mu*wts^.5/sqrt(object$sig2*object$family$variance(mu)),
              pearson = y.mu*wts^.5/sqrt(object$family$variance(mu)),
              deviance = { d.res <- sqrt(pmax(object$family$dev.resids(y,mu,wts),0))
                           ifelse(y>mu , d.res, -d.res)             
                         },
              response = y.mu)
  res <- naresid(object$na.action,res)
  res
}

