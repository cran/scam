## (c) Natalya Pya (2012-2024). Provided under GPL 2.
## based on (c) Simon N Wood predict.gam(mgcv) ...


predict.scam <- function(object,newdata,type="link",se.fit=FALSE,terms=NULL,exclude=NULL,
                       block.size=NULL,newdata.guaranteed=FALSE,na.action=na.pass,...) 
{
# This function is used for predicting from a SCAM. object is a scam object, newdata a dataframe to
# be used in prediction......
#
# Type == "link"     - for linear predictor
#      == "response" - for fitted values
#      == "terms"    - for individual terms on scale of linear predictor 
#      == "iterms"   - exactly as "terms" except that se's include uncertainty about mean for unconstrained smooths
#      == "lpmatrix" - for matrix mapping parameters to l.p.
# Steps are:
#  1. Set newdata to object$model if no newdata supplied
#  2. split up newdata into manageable blocks if too large
#  3. Obtain parametric model matrix (safely)
#  4. Work through smooths calling prediction.matrix constructors for each term
#  5. Work out required quantities
# 
# The splitting into blocks enables blocks of compiled code to be called efficiently
# using smooth class specific prediction matrix constructors, without having to 
# build up potentially enormous prediction matrices.
# if newdata.guaranteed == TRUE then the data.frame is assumed complete and
# ready to go, so that only factor levels are checked for sanity.
# 
# if `terms' is non null then it should be a list of terms to be returned 
# when type=="terms". 
# if `object' has an attribute `para.only' then only parametric terms of order
# 1 are returned for type=="terms": i.e. only what termplot can handle.
#
# if no new data is supplied then na.action does nothing, otherwise 
# if na.action == "na.pass" then NA predictors result in NA predictions (as lm
#                   or glm)
#              == "na.omit" or "na.exclude" then NA predictors result in
#                       dropping
# if GC is TRUE then gc() is called after each block is processed
# exclude: if type=="terms" or type="iterms" then terms (smooth or parametric) named
# in this array will not be returned. Otherwise any smooth terms named in this
# array will be set to zero. If NULL then no terms are excluded.

  if (type!="link"&&type!="terms"&&type!="iterms"&&type!="response"&&type!="lpmatrix")  
  { warning("Unknown type, reset to terms.")
    type<-"terms"
  }
  if (!inherits(object,"scam")) stop("predict.scam can only be used to predict from scam objects")

#  if (ncol(attr(object$terms,"factors")) == 1)
#      { if (max(newdata) > max(object$data[,attr(b$terms,"term.labels")])) 
#           stop("predict.scam can only be used for data within #the range of observed values, please use extrapolate.scam #otherwise")  }

  ## to mimic behaviour of predict.lm, some resetting is required ...
  if (missing(newdata)) na.act <- object$na.action else {
    if (is.null(na.action)) na.act <- NULL 
    else {
      na.txt <- if (is.character(na.action)||is.function(na.action)) get.na.action(na.action) else "na.pass"
      #na.txt <- "na.pass"
      #if (is.character(na.action))
      #na.txt <- substitute(na.action) else
      #if (is.function(na.action)) na.txt <- deparse(substitute(na.action))
      if (na.txt=="na.pass") na.act <- "na.exclude" else
      if (na.txt=="na.exclude") na.act <- "na.omit" else
      na.act <- na.action
    }
  } ## ... done

 # get data from which to predict.....  
  nd.is.mf <- FALSE # need to flag if supplied newdata is already a model frame 
  ## get name of response...
  # yname <- all.vars(object$terms)[attr(object$terms,"response")]
  yname <- attr(attr(object$terms,"dataClasses"),"names")[attr(object$terms,"response")]
  if (newdata.guaranteed==FALSE) {
    if (missing(newdata)) { # then "fake" an object suitable for prediction 
      newdata <- object$model
      new.data.ok <- FALSE
      nd.is.mf <- TRUE
      response <- newdata[[yname]] ## ok even with "cbind(foo,bar)" as yname 
    } else {  # do an R ``standard'' evaluation to pick up data
      new.data.ok <- TRUE
      if (is.data.frame(newdata)&&!is.null(attr(newdata,"terms"))) { # it's a model frame
        if (sum(!(names(object$model)%in%names(newdata)))) stop(
        "newdata is a model.frame: it should contain all required variables\n")
         nd.is.mf <- TRUE
      } else {
        ## Following is non-standard to allow convenient splitting into blocks
        ## below, and to allow checking that all variables are in newdata ...

        ## get names of required variables, less response, but including offset variable
        ## see ?terms.object and ?terms for more information on terms objects
        # yname <- all.vars(object$terms)[attr(object$terms,"response")] ## redundant
        resp <- get.var(yname,newdata,FALSE)
        naresp <- FALSE
        #if (!is.null(object$family$predict)&&!is.null(newdata[[yname]])) {
	if (!is.null(object$family$predict)&&!is.null(resp)) {
          ## response provided, and potentially needed for prediction (e.g. Cox PH family): not in scam yet
          if (!is.null(object$pred.formula)) object$pred.formula <- attr(object$pred.formula,"full")
          response <- TRUE
          Terms <- terms(object)
          #resp <- newdata[[yname]]
	  if (is.matrix(resp)) {
            if (sum(is.na(rowSums(resp)))>0) stop("no NAs allowed in response data for this model")
          } else { ## vector response
            if (sum(is.na(resp))>0) {
              naresp <- TRUE ## there are NAs in supplied response
              ## replace them with a numeric code, so that rows are not dropped below
              rar <- range(resp,na.rm=TRUE)
              thresh <- rar[1]*1.01-rar[2]*.01
              resp[is.na(resp)] <- thresh
              newdata[[yname]] <- thresh 
            }
	  }  
        } else { ## response not provided
          response <- FALSE 
          Terms <- delete.response(terms(object))
        }
        allNames <- if (is.null(object$pred.formula)) all.vars(Terms) else all.vars(object$pred.formula)
        if (length(allNames) > 0) { 
          ff <- if (is.null(object$pred.formula)) reformulate(allNames) else  object$pred.formula
          if (sum(!(allNames%in%names(newdata)))) { 
            warning("not all required variables have been supplied in  newdata!\n")
          }
          ## note that `xlev' argument not used here, otherwise `as.factor' in 
          ## formula can cause a problem ... levels reset later.
          newdata <- eval(model.frame(ff,data=newdata,na.action=na.act),parent.frame())
          if (naresp) newdata[[yname]][newdata[[yname]]<=thresh] <- NA ## reinstate as NA  
        } ## otherwise it's intercept only and newdata can be left alone
        na.act <- attr(newdata,"na.action")
        #response <- if (response) newdata[[yname]] else NULL
	response <- if (response) get.var(yname,newdata,FALSE) else NULL
      }
    }
  } else { ## newdata.guaranteed == TRUE
    na.act <- NULL
    new.data.ok=TRUE ## it's guaranteed!
    if (!is.null(attr(newdata,"terms"))) nd.is.mf <- TRUE
    #response <- newdata[[yname]]
    response <- get.var(yname,newdata,FALSE)
  }


  ## now check the factor levels and split into blocks...

  if (new.data.ok){
    ## check factor levels are right ...
    names(newdata)->nn # new data names
    colnames(object$model)->mn # original names
    for (i in 1:length(newdata)) 
    if (nn[i]%in%mn && is.factor(object$model[,nn[i]])){ # then so should newdata[[i]] be 
       ## newdata[[i]]<-factor(newdata[[i]],levels=levels(object$model[,nn[i]])) # set prediction levels to fit levels
      levm <- levels(object$model[,nn[i]]) ## original levels
      levn <- levels(factor(newdata[[i]])) ## new levels
      if (sum(!levn%in%levm)>0) { ## check not trying to sneak in new levels 
        msg <- paste("factor levels",paste(levn[!levn%in%levm],collapse=", "),"not in original fit",collapse="")
        warning(msg)
      }
      ## set prediction levels to fit levels...
      if (is.matrix(newdata[[i]])) {
        dum <- factor(newdata[[i]],levels=levm)
	dim(dum) <- dim(newdata[[i]])
	newdata[[i]] <- dum
      } else newdata[[i]] <- factor(newdata[[i]],levels=levm)
    }
    if (type=="newdata") return(newdata)

    # split prediction into blocks, to avoid running out of memory
    if (length(newdata)==1) newdata[[2]]<-newdata[[1]] # avoids data frame losing its labels and dimensions below!
    if (is.null(dim(newdata[[1]]))) np<-length(newdata[[1]]) 
    else np <- dim(newdata[[1]])[1] 
    nb <- length(object$coefficients)
    if (is.null(block.size)) block.size <- 1000
    if (block.size < 1) block.size <- np
  } else { # no new data, just use object$model
    np <- nrow(object$model)
    nb <- length(object$coefficients)
  }

  if (type=="lpmatrix") block.size <- NULL ## nothing gained by blocking in this case - and offset handling easier this way

  ## split prediction into blocks, to avoid running out of memory
  if (is.null(block.size)) { 
    ## use one block as predicting using model frame
    ## and no block size supplied... 
    n.blocks <- 1
    b.size <- array(np,1)
  } else {
    n.blocks <- np %/% block.size
    b.size <- rep(block.size,n.blocks)
    last.block <- np-sum(b.size)
    if (last.block>0) {
      n.blocks <- n.blocks+1  
      b.size[n.blocks] <- last.block
    }
  }

  # setup prediction arrays...
  n.smooth<-length(object$smooth)
  if (type=="lpmatrix"){
   H<-matrix(0,np,nb)
  } else if (type=="terms"||type=="iterms"){
    term.labels<-attr(object$pterms,"term.labels")
    para.only <- attr(object,"para.only")
    if (is.null(para.only)) para.only <- FALSE  # if TRUE then only return information on parametric part
    n.pterms <- length(term.labels)
    fit<-array(0,c(np,n.pterms+as.numeric(para.only==0)*n.smooth))
    if (se.fit) se<-fit
    ColNames<-term.labels
  } else { ## "response" or "link"
    fit<-array(0,np)
    if (se.fit) se<-fit
  }
  stop<-0

  ####################################
  ## Actual prediction starts here...
  ####################################

 # Terms <- list(delete.response(object$pterms))
  Terms <- delete.response(object$pterms)

  s.offset <- NULL # to accumulate any smooth term specific offset
  any.soff <- FALSE # indicator of term specific offset existence
  if (n.blocks > 0) for (b in 1:n.blocks){  # work through prediction blocks
    start<-stop+1
    stop<-start+b.size[b]-1
    if (n.blocks==1) data <- newdata  else data<-newdata[start:stop,]
    X <- matrix(0,b.size[b],nb)
    Xoff <- matrix(0,b.size[b],n.smooth) ## term specific offsets 
    ## implements safe prediction for parametric part as described in
    ## http://developer.r-project.org/model-fitting-functions.txt
    if (new.data.ok){
      if (nd.is.mf) mf <- model.frame(data,xlev=object$xlevels) 
      else {
        mf <- model.frame(Terms,data,xlev=object$xlevels)
        if (!is.null(cl <- attr(object$pterms,"dataClasses"))) .checkMFClasses(cl,mf)
      } 
       ## next line is just a work around to prevent a spurious warning (e.g. R 3.6) from
       ## model.matrix if contrast relates to a term in mf which is not
       ## part of Terms (model.matrix doc actually defines contrast w.r.t. mf,
       ## not Terms)...
        # Xp <- model.matrix(Terms,mf,contrasts=object$contrasts) 
      oc <- if (length(object$contrasts)==0) object$contrasts else
	      object$contrasts[names(object$contrasts)%in%attr(Terms,"term.labels")]
      Xp <- model.matrix(Terms,mf,contrasts=oc)
    } else {
       Xp <- model.matrix(Terms,object$model)
       mf <- newdata # needed in case of offset, below
    }
    
    if (object$nsdf) X[,1:object$nsdf]<-Xp

    if (n.smooth) for (k in 1:n.smooth) { ## loop through smooths
        klab <- object$smooth[[k]]$label
        if ((is.null(terms)||(klab%in%terms))&&(is.null(exclude)||!(klab%in%exclude))) {
          Xfrag <- PredictMat(object$smooth[[k]],data)	
          if (!is.matrix(Xfrag)) Xfrag <- matrix(Xfrag,nrow=nrow(data))
          ## added code specific for scam....
          if (inherits(object$smooth[[k]], c("mpi.smooth","mpic.smooth","mpd.smooth", "cv.smooth", "cx.smooth",
                 "mdcv.smooth","mdcx.smooth","micv.smooth","micx.smooth","po.smooth",  
                 "ipo.smooth","cpopspline.smooth",
                 "tedmi.smooth","tedmd.smooth","temicx.smooth","temicv.smooth", "tedecx.smooth",
                 "tedecv.smooth","tecvcv.smooth","tecxcx.smooth","tecxcv.smooth")))
               X[,object$smooth[[k]]$first.para:object$smooth[[k]]$last.para] <- Xfrag[,2:ncol(Xfrag)]
            else if (inherits(object$smooth[[k]], c("dpo.smooth"))) 
               X[,object$smooth[[k]]$first.para:object$smooth[[k]]$last.para] <- Xfrag[,1:(ncol(Xfrag)-1)]
            else if (inherits(object$smooth[[k]], c("miso.smooth","mifo.smooth"))) ## 'zero start' and  'zero end' increasing constraints
               X[,object$smooth[[k]]$first.para:object$smooth[[k]]$last.para] <- Xfrag[,-object$smooth[[k]]$n.zero]
          #  else if (inherits(object$smooth[[k]], c("mipoc.smooth"))) ## 'pss zero ' increasing constraint ...
           #    X[,object$smooth[[k]]$first.para:object$smooth[[k]]$last.para] <- Xfrag
            else if (inherits(object$smooth[[k]], c("tesmi1.smooth", "tesmi2.smooth", "tesmd1.smooth",
                  "tesmd2.smooth","tescx.smooth","tescv.smooth"))) { ## for single monotonicity/ convexity             
                  X[,object$smooth[[k]]$first.para:object$smooth[[k]]$last.para] <- Xfrag%*%object$smooth[[k]]$Zc               
                ##  X[,object$smooth[[k]]$first.para:object$smooth[[k]]$last.para] <- sweep(X[,object$smooth[[k]]$first.para:object$smooth[[k]]$last.para],2,object$smooth[[k]]$cmX)
               } else if (inherits(object$smooth[[k]], c("tismi.smooth","tismd.smooth"))) {## for smooth interaction             
                   X[,object$smooth[[k]]$first.para:object$smooth[[k]]$last.para] <- Xfrag %*%object$smooth[[k]]$Zc 
                } else if (inherits(object$smooth[[k]], c("lipl.smooth"))) {## local scop-spline increasing up to given point and a plateau from it               
                   n.extra.col <- object$smooth[[k]]$n.zero.col ## ncol(Xfrag) - (object$smooth[[k]]$last.para-object$smooth[[k]]$first.para+1)
                   ind.extra <- c((object$smooth[[k]]$last.para+1):(ncol(Xfrag)+1))
                   X <- cbind(X,matrix(0,nrow(X), n.extra.col))
                   X[,object$smooth[[k]]$first.para:(object$smooth[[k]]$first.para+ncol(Xfrag)-1)] <- Xfrag                                   
               } else ## unconstrainded smooths, 'by' constrained and local scop smooths...
                    X[,object$smooth[[k]]$first.para:object$smooth[[k]]$last.para] <- Xfrag
          Xfrag.off <- attr(Xfrag,"offset") ## any term specific offsets?
          if (!is.null(Xfrag.off)) { Xoff[,k] <- Xfrag.off; any.soff <- TRUE }
        }
        if (type=="terms"||type=="iterms") ColNames[n.pterms+k] <- klab
    } ## smooths done

    # Now have prediction matrix for this block, now do something with it
    if (type=="lpmatrix") { 
      if (inherits(object$smooth[[k]], c("lipl.smooth"))) {
          H <- cbind(H,matrix(0,nrow(H), n.extra.col))
          H[start:stop,] <- X
      } else H[start:stop,] <- X
      if (any.soff) s.offset <- rbind(s.offset,Xoff)
    } else if (type=="terms" ||type=="iterms") {
      ind <- 1:length(object$assign)
      if (n.pterms)  # work through parametric part
      for (i in 1:n.pterms){
        ii <- ind[object$assign==i]
        ## CORRECTIONS FOR SCAM....
        fit[start:stop,i] <- as.matrix(X[,ii,drop=FALSE])%*%object$coefficients.t[ii]
        if (se.fit) se[start:stop,i] <-
               sqrt(rowSums((as.matrix(X[,ii,drop=FALSE])%*%object$Vp.t[ii,ii])*as.matrix(X[,ii,drop=FALSE])))
      }

      if (n.smooth&&!para.only) { 
        for (k in 1:n.smooth) # work through the smooth terms 
        { first<-object$smooth[[k]]$first.para;last<-object$smooth[[k]]$last.para

          ## CORRECTED for SCAM ...
          fit[start:stop,n.pterms+k]<-X[,first:last,drop=FALSE]%*%object$coefficients.t[first:last] + Xoff[,k]
          if (se.fit) { # diag(Z%*%V%*%t(Z))^0.5; Z=X[,first:last]; V is sub-matrix of Vp
               if (inherits(object$smooth[[k]], c("mpi.smooth","mpd.smooth","cv.smooth", "cx.smooth",
                               "mdcv.smooth","mdcx.smooth","micv.smooth","micx.smooth","po.smooth","cpopspline.smooth"))){
                      if (nrow(X)==1) # prediction vector if prediction is made for only one value of covariates
                              X1 <- c(1,t(X[,first:last]))
                          else 
                              X1 <- cbind(rep(1,nrow(X)),X[,first:last]) # prediction matrix
                      Vp <- object$Vp.t[c(1,first:last),c(1,first:last)] 
                      Vp[,1] <- rep(0,nrow(Vp))
                      Vp[1,] <- rep(0,ncol(Vp))
                      se[start:stop,n.pterms+k] <- sqrt(rowSums((X1%*%Vp)*X1))
               } else if (inherits(object$smooth[[k]], c("ipo.smooth","dpo.smooth"))){
                           if (nrow(X)==1) # prediction vector if prediction is made for only one value of covariates
                              X1 <- if  (inherits(object$smooth[[k]], c("dpo.smooth"))) c(t(X[,first:last]),1)
                                    else c(1,t(X[,first:last]))
                           else 
                              X1 <- if  (inherits(object$smooth[[k]], c("dpo.smooth"))) cbind(X[,first:last],rep(1,nrow(X))) 
                                    else cbind(rep(1,nrow(X)),X[,first:last]) # prediction matrix
                           # X0 - model matrix of the original data....
                           object.X <- model.matrix(object)
                           X0 <- if  (inherits(object$smooth[[k]], c("dpo.smooth"))) cbind(object.X[,first:last], rep(1,nrow(object.X))) 
                                 else cbind(rep(1,nrow(object.X)),object.X[,first:last]) 
                           q <- ncol(X0)
                           onet <- matrix(rep(1,nrow(X0)),1,nrow(X0))
                           A <- onet%*%X0
                           qrX <- qr(X0)                
                           R <- qr.R(qrX) 
                           qrA <- qr(t(A))
                           if (inherits(object$smooth[[k]], c("dpo.smooth"))) {
                                      R <- R[-q,]
                                      RZa <- t(qr.qty(qrA,t(R)))[,1:(q-1)] 
                            } else {  R <- R[-1,]
                                      RZa <- t(qr.qty(qrA,t(R)))[,2:q] 
                                   }   
                            RZa.inv <- solve(RZa)
                            RZaR <- RZa.inv%*%R
                            if (nrow(X)==1)  
                                 XZa <- if  (inherits(object$smooth[[k]], c("dpo.smooth"))) t(qr.qty(qrA,X1))[,1:(q-1)]
                                        else  t(qr.qty(qrA,X1))[,2:q]          
                               else   
                                   XZa <- if  (inherits(object$smooth[[k]], c("dpo.smooth"))) t(qr.qty(qrA,t(X1)))[,1:(q-1)]
                                           else  t(qr.qty(qrA,t(X1)))[,2:q] 
                                   
                            Ga <- XZa%*%RZaR
                            if (inherits(object$smooth[[k]], c("dpo.smooth"))) {
                                  Vp <- object$Vp.t[c(first:last,1),c(first:last,1)] 
                                  Vp[,q] <- rep(0,nrow(Vp))
                                  Vp[q,] <- rep(0,ncol(Vp))
                             } else {
                                  Vp <- object$Vp.t[c(1,first:last),c(1,first:last)] 
                                  Vp[,1] <- rep(0,nrow(Vp))
                                  Vp[1,] <- rep(0,ncol(Vp))
                               }
                             se[start:stop,n.pterms+k] <- sqrt(pmax(0,rowSums((Ga%*%Vp)*Ga)))                                                       
               } else if (inherits(object$smooth[[k]], c("tedmi.smooth","tedmd.smooth","tesmi1.smooth",
                                    "tismi.smooth","tismd.smooth","tesmi2.smooth", "tesmd1.smooth", "tesmd2.smooth","temicx.smooth", "temicv.smooth",
                                   "tedecx.smooth", "tedecv.smooth", "tescx.smooth", 
                                   "tescv.smooth","tecvcv.smooth","tecxcx.smooth","tecxcv.smooth"))) {
                      # X0 - model matrix of the original data....
                      object.X <- model.matrix(object)
                      X0 <- cbind(rep(1,nrow(object.X)),object.X[,first:last]) 
                      onet <- matrix(rep(1,nrow(X0)),1,nrow(X0))
                      if (nrow(X)==1) # prediction vector if prediction is made for only one value of covariates
                            X1 <- c(1,t(X[,first:last]))
                        else 
                            X1 <- cbind(rep(1,nrow(X)),X[,first:last]) # prediction matrix
                      A <- onet%*%X0
                      qrX <- qr(X0)
                      R <- qr.R(qrX)
                      qrA <- qr(t(A))
                      q <- ncol(X0)
                      if (inherits(object$smooth[[k]], c("tedmi.smooth", "tedmd.smooth","tedmdc.smooth", "temicx.smooth",
                            "temicv.smooth", "tedecx.smooth","tedecv.smooth", "tecvcv.smooth",
                            "tecxcx.smooth","tecxcv.smooth"))) { # get RZaR for double monotonicity...
                            R <- R[-1,]
                            RZa <- t(qr.qty(qrA,t(R)))[,2:q] 
                            RZa.inv <- solve(RZa)
                            RZaR <- RZa.inv%*%R
                          }
                       else { # get RZaR for single monotonicity...
                            RZa <- t(qr.qty(qrA,t(R)))[,2:q]
                            RZatRZa.inv <- solve(crossprod(RZa)) ## solve(t(RZa)%*%RZa) 
                            Q <- qr.Q(qrX)
                            B1 <-  tcrossprod(RZatRZa.inv,RZa) ## RZatRZa.inv%*%t(RZa)
                            RZaR <- B1%*%R
                           }
                       if (nrow(X)==1)
                              XZa <- t(qr.qty(qrA,X1))[,2:ncol(X1)]
                       else
                              XZa <- t(qr.qty(qrA,t(X1)))[,2:ncol(X1)]
                       Ga <- XZa%*%RZaR%*%object$smooth[[k]]$Zc
                       Vp <- object$Vp.t[first:last,first:last]
                       se[start:stop,n.pterms+k] <- sqrt(rowSums((Ga%*%Vp)*Ga))
               } else if (inherits(object$smooth[[k]], c("miso.smooth", "mifo.smooth"))) { ## 'zero start' and passing through zero increasing constraint ...
                       se[start:stop,n.pterms+k] <- sqrt(rowSums((X[,first:last,drop=FALSE]%*%object$Vp.t[first:last,first:last])*X[,first:last,drop=FALSE]))
               } else { ## for local scop smooth terms, scop with'by' and unconstrained smooth terms..... 
                    if (type=="iterms"&& attr(object$smooth[[k]],"nCons")>0) { ## termwise se to "carry the intercept
                      X1 <- matrix(object$cmX,nrow(X),ncol(X),byrow=TRUE)
                      meanL1 <- object$smooth[[k]]$meanL1
                      if (!is.null(meanL1)) X1 <- X1 / meanL1              
                      X1[,first:last] <- X[,first:last,drop=FALSE]
                      se[start:stop,n.pterms+k] <- sqrt(pmax(0,rowSums((X1%*%object$Vp)*X1)))
                    } else  se[start:stop,n.pterms+k] <- ## terms strictly centred
                    sqrt(rowSums((X[,first:last,drop=FALSE]%*%object$Vp.t[first:last,first:last])*X[,first:last,drop=FALSE]))
                 }         
          } ## end if (se.fit)
        } ## end # work through the smooth terms 
        colnames(fit) <- ColNames
        if (se.fit) colnames(se) <- ColNames
      } else { # para.only
        if (para.only&&is.list(object$pterms)) { 
            ## have to use term labels that match original data, or termplot fails 
            ## to plot. This only applies for 'para.only==1' calls which are 
            ## designed for use from termplot called from plot.gam
            term.labels <- unlist(lapply(object$pterms,attr,"term.labels"))
          }
        colnames(fit) <- term.labels
        if (se.fit) colnames(se) <- term.labels
        if (para.only) { 
           # retain only terms of order 1 - this is to make termplot work
           order <- if (is.list(object$pterms)) unlist(lapply(object$pterms,attr,"order")) else attr(object$pterms,"order") #attr(object$pterms,"order")
           term.labels <- term.labels[order==1]
           ## fit <- as.matrix(as.matrix(fit)[,order==1])
           fit <- fit[,order==1,drop=FALSE]
           colnames(fit) <- term.labels
           if (se.fit) { ## se <- as.matrix(as.matrix(se)[,order==1])
              se <- se[,order==1,drop=FALSE] 
              colnames(se) <- term.labels } 
           }
        }
       ##if (!is.null(terms)) { # return only terms requested via `terms'
       ## if (sum(!(terms %in%colnames(fit)))) 
       ## warning("non-existent terms requested - ignoring")
       ## else { names(term.labels) <- term.labels
       ##   term.labels <- term.labels[terms]  # names lost if only one col
       ##   fit <- as.matrix(as.matrix(fit)[,terms])
       ##   colnames(fit) <- term.labels
       ##   if (se.fit) {se <- as.matrix(as.matrix(se)[,terms])
       ##   colnames(se) <- term.labels}
       ## }
       ##}
      
    } else {# "link" or "response"
        k<-attr(attr(object$model,"terms"),"offset")
        ## CORRECTED for SCAM...
        fit[start:stop]<-X%*%object$coefficients.t + rowSums(Xoff)
        if (!is.null(k)) fit[start:stop]<-fit[start:stop]+model.offset(mf) + rowSums(Xoff)
        if (se.fit) se[start:stop]<-sqrt(rowSums((X%*%object$Vp.t)*X))
        if (type=="response") {# transform    
           fam<-object$family;linkinv <- fam$linkinv;dmu.deta<-fam$mu.eta  
           if (se.fit) se[start:stop] <- se[start:stop]*abs(dmu.deta(fit[start:stop])) 
           fit[start:stop] <- linkinv(fit[start:stop])
        }
      }
    rm(X)   
  } ## end of prediction block loop

  if ((type=="terms"||type=="iterms")&&(!is.null(terms)||!is.null(exclude))) { # return only terms requested via `terms'
    cnames <- colnames(fit)
    if (!is.null(terms)) {
      if (sum(!(terms %in%cnames))) 
        warning("non-existent terms requested - ignoring")
      else { 
        fit <- fit[,terms,drop=FALSE]
        if (se.fit) {
           se <- se[,terms,drop=FALSE]
        }
      }
    }
    if (!is.null(exclude)) {
      if (sum(!(exclude %in%cnames))) 
        warning("non-existent exclude terms requested - ignoring")
      else { 
        exclude <- which(cnames%in%exclude) ## convert to numeric column index
        fit <- fit[,-exclude,drop=FALSE]
        if (se.fit) {
          se <- se[,-exclude,drop=FALSE]
        }
      }
    }
  }

  #if (type=="response"&&!is.null(fit1)) {
  #  fit <- fit1
  #  if (se.fit) se <- se1
  #}


  rn <- rownames(newdata)
  if (type=="lpmatrix") { 
    if (inherits(object$smooth[[k]], c("lipl.smooth"))) 
             colnames(H[,-ind.extra]) <- names(object$coefficients)
        else colnames(H) <- names(object$coefficients)
    rownames(H)<-rn
    if (!is.null(s.offset)) { 
      s.offset <- napredict(na.act,s.offset)
      attr(H,"offset") <- s.offset ## term specific offsets...
    }
    H <- napredict(na.act,H)
  } else { 
    if (se.fit) { 
      if (is.null(nrow(fit))) {
        names(fit) <- rn
        names(se) <- rn
        fit <- napredict(na.act,fit)
        se <- napredict(na.act,se) 
      } else { 
        rownames(fit)<-rn
        rownames(se)<-rn
        fit <- napredict(na.act,fit)
        se <- napredict(na.act,se)
      }
      H<-list(fit=fit,se.fit=se) 
    } else { 
      H<-fit
      if (is.null(nrow(H))) names(H) <- rn else
      rownames(H)<-rn
      H <- napredict(na.act,H)
    }
  }
  if ((type=="terms"||type=="iterms")&&attr(object$terms,"intercept")==1) attr(H,"constant") <- object$coefficients[1]
  H # ... and return
} ## end predict.scam


################################
## below is from mgcv r ...
get.na.action <- function(na.action) {
## get the name of the na.action whether function or text string.
## avoids deparse(substitute(na.action)) which is easily broken by 
## nested calls.
  if (is.character(na.action)) {
    if (na.action%in%c("na.omit","na.exclude","na.pass","na.fail")) 
      return(na.action) else stop("unrecognised na.action")
  } 
  if (!is.function(na.action)) stop("na.action not character or function")
  a <- try(na.action(c(0,NA)),silent=TRUE)
  if (inherits(a,"try-error")) return("na.fail")
  if (inherits((attr(a,"na.action")),"omit")) return("na.omit")
  if (inherits((attr(a,"na.action")),"exclude")) return("na.exclude")
  return("na.pass")
} ## get.na.action

###############################################################################
### ISSUES.....
#
#* predict.scam "terms" and "iterms" don't deal properly with 
#  shape constrained smooths: se.fit retured are the same for both types

