## incorporating the emmeans package (copyright (c) 2012-2025 Russell V. Lenth ) support in the scam from version 1.2-21
## (c) Natalya Pya (2026). Provided under GPL 2

###################################################################################
## from the 'Extending *emmeans* vignette...                                     ##
## need to write two S3 methods, 'recover_data' and 'emm_basis', for the class   ##
## of object that your model-fitting function returns.                           ##
## i) The recover_data method is needed to recreate the dataset so that          ##
## the reference grid can be identified.                                         ##
## ii) The emm_basis method then determines the linear functions needed          ##
## to evaluate each point in the reference grid and to obtain associated         ##  
## information such as the variance-covariance matrix needed to do estimation    ## 
## and testing.                                                                  ##
## iii) 'emm_basis'calls a utility function 'smooth.is.random- to sort out which ##
## smoothers are associated with random effects, needed because emmeans focuses  ##
## on fixed effects                                                              ##
###################################################################################


### for scam::scam objects
### similar to recover_data gam...
recover_data.scam = function(object, ...) 
## recreates the dataset so that the reference grid can be identified.    
{
    if (length(object$smooth) > 0) { # get rid of random terms
        fixnm = unlist(lapply(object$smooth, function(s) {
                              if(.smooth.is.random(s)) ""
                              else c(s$term, s$by)
            }))
        fixnm = union(emmeans::.all.vars(delete.response(object$pterms)), fixnm)
        fixnm = setdiff(fixnm, c("1", "", "NA"))
        object$terms = terms(.reformulate(fixnm, env = environment(terms(object)))) 
    }
    recover_data.lm(object, ...)
}


### emm_basis method for scam::scam objects
### extra arg `freq` and `untransformed' as in `vcov.scam`          
emm_basis.scam = function(object, trms, xlev, grid,
                         freq = FALSE, untransformed = FALSE, ...) 
## Obtains six things and return them in a named list:
## 1) matrix 'X' of linear functions for each point in the reference grid,
## 2) regression coefficients 'bhat',
## 3) variance-covariance matrix 'V',
## 4) matrix 'nbasis' for non-estimable functions,
## 5) function 'dffun(k,dfargs)' for computing degrees of freedom for the linear function sum(k*bhat),
## 6) a list 'dfargs' of arguments to pass to 'dffun'.                       
{
    if (length(object$smooth) > 0) { # get rid of random terms 
        rand = sapply(object$smooth, function(s) {ifelse(.smooth.is.random(s), s$label, NA)})
        rand = if (all(is.na(rand))) NULL else rand[!is.na(rand)]
    }
    else
        rand = NULL
    X = predict.scam(object, newdata = grid, type = "lpmatrix",
                          exclude = rand, newdata.guaranteed = TRUE)
    keep = if (is.null(rand)) rep(TRUE, ncol(X)) else apply(X, 2, function(x) !all(x == 0))
    X = X[, keep, drop = FALSE]
    bhat = as.numeric(object$coefficients.t[keep])
    V = emmeans::.my.vcov(object, freq = freq, untransformed = untransformed, ...)[keep, keep]
  
    nbasis = estimability::all.estble
    dfargs = list(df = object$df.residual)
    dffun = function(k, dfargs) dfargs$df
    list(X=X, bhat=bhat, nbasis=nbasis, V=V, dffun=dffun, dfargs=dfargs)
}


# Local utility for identifying random smooths
.smooth.is.random = function(s) {
    rcls = c("random.effect", "fs.interaction") ### what to look for
    cls = c(class(s), unname(unlist(lapply(s$margin, class))))
    any(cls %in% rcls)
}

###############################################################################
# below are internal (unexported) functions of the emmeans package, needed    #
# for running recover_data.scam()                                             #
#    Copyright (c) 2012-2019 Russell V. Lenth                                 #
###############################################################################
### Contributed by Jonathon Love, https://github.com/jonathon-love          ###
### and adapted by RVL to exclude terms like df$trt or df[["trt"]]          ###
###############################################################################
# reformulate internally in emmeans
# same as stats::reformulate, except it surrounds term labels with backsticks

.reformulate <- function (termlabels, response = NULL, intercept = TRUE, env = parent.frame())
{
    if (!is.character(termlabels) || !length(termlabels))
        stop("'termlabels' must be a character vector of length at least one")
    has.resp = !is.null(response)
    termlabels = sapply(trimws(termlabels), function(x) 
        if (length(grep("\\$|\\[\\[", x)) > 0) x
        else paste0("`", x, "`"))
    termtext = paste(if (has.resp) "response", "~", 
                     paste(termlabels, collapse = "+"), collapse = "")
# prev version:                     paste0("`", trimws(termlabels), "`", collapse = "+"), collapse = "")
    if (!intercept)
        termtext = paste(termtext, "- 1")
    rval = eval(parse(text = termtext, keep.source = FALSE)[[1L]])
    if (has.resp)
        rval[[2L]] = if (is.character(response))
            as.symbol(response)
    else response
    environment(rval) = env
    rval
}

recover_data.lm = function(object, frame = object$model, ...) {
        fcall = object$call
    emmeans::recover_data(fcall, delete.response(terms(object)), object$na.action, 
                 frame = frame, pwts = weights(object), ...)
}





