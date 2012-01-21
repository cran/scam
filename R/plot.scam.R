
## PLOT FUNCTION based on plot.gam(mgcv)....

plot.scam <- function (x, residuals = FALSE, rug = TRUE, se = TRUE, 
    pages = 0, 
    select = NULL, scale = -1, n = 100, n2 = 40, pers = FALSE, 
    theta = 30, phi = 30, jit = FALSE, xlab = NULL, ylab = NULL, 
    main = NULL, ylim = NULL, xlim = NULL, too.far = 0.1, all.terms = FALSE, 
    shade = FALSE, shade.col = "gray80", shift = 0, trans = I, 
    seWithMean = FALSE, by.resids = FALSE, ...) 
{
    sub.edf <- function(lab, edf) {
        pos <- regexpr(":", lab)[1]
        if (pos < 0) {
            pos <- nchar(lab) - 1
            lab <- paste(substr(lab, start = 1, stop = pos), 
                ",", round(edf, digits = 2), ")", sep = "")
        }
        else {
            lab1 <- substr(lab, start = 1, stop = pos - 2)
            lab2 <- substr(lab, start = pos - 1, stop = nchar(lab))
            lab <- paste(lab1, ",", round(edf, digits = 2), lab2, 
                sep = "")
        }
        lab
    }
    sp.contour <- function(x, y, z, zse, xlab = "", ylab = "", 
        zlab = "", titleOnly = FALSE, se.plot = TRUE, se.mult = 1, 
        trans = I, shift = 0, ...) {
        gap <- median(zse, na.rm = TRUE)
        zr <- max(trans(z + zse + shift), na.rm = TRUE) - min(trans(z - 
            zse + shift), na.rm = TRUE)
        n <- 10
        while (n > 1 && zr/n < 2.5 * gap) n <- n - 1
        zrange <- c(min(trans(z - zse + shift), na.rm = TRUE), 
            max(trans(z + zse + shift), na.rm = TRUE))
        zlev <- pretty(zrange, n)
        yrange <- range(y)
        yr <- yrange[2] - yrange[1]
        xrange <- range(x)
        xr <- xrange[2] - xrange[1]
        ypos <- yrange[2] + yr/10
        args <- as.list(substitute(list(...)))[-1]
        args$x <- substitute(x)
        args$y <- substitute(y)
        args$type = "n"
        args$xlab <- args$ylab <- ""
        args$axes <- FALSE
        do.call("plot", args)
        cs <- (yr/10)/strheight(zlab)
        if (cs > 1) 
            cs <- 1
        tl <- strwidth(zlab)
        if (tl * cs > 3 * xr/10) 
            cs <- (3 * xr/10)/tl
        args <- as.list(substitute(list(...)))[-1]
        n.args <- names(args)
        zz <- trans(z + shift)
        args$x <- substitute(x)
        args$y <- substitute(y)
        args$z <- substitute(zz)
        if (!"levels" %in% n.args) 
            args$levels <- substitute(zlev)
        if (!"lwd" %in% n.args) 
            args$lwd <- 2
        if (!"labcex" %in% n.args) 
            args$labcex <- cs * 0.65
        if (!"axes" %in% n.args) 
            args$axes <- FALSE
        if (!"add" %in% n.args) 
            args$add <- TRUE
        do.call("contour", args)
        if (is.null(args$cex.main)) 
            cm <- 1
        else cm <- args$cex.main
        if (titleOnly) 
            title(zlab, cex.main = cm)
        else {
            xpos <- xrange[1] + 3 * xr/10
            xl <- c(xpos, xpos + xr/10)
            yl <- c(ypos, ypos)
            lines(xl, yl, xpd = TRUE, lwd = args$lwd)
            text(xpos + xr/10, ypos, zlab, xpd = TRUE, pos = 4, 
                cex = cs * cm, off = 0.5 * cs * cm)
        }
        if (is.null(args$cex.axis)) 
            cma <- 1
        else cma <- args$cex.axis
        axis(1, cex.axis = cs * cma)
        axis(2, cex.axis = cs * cma)
        box()
        if (is.null(args$cex.lab)) 
            cma <- 1
        else cma <- args$cex.lab
        mtext(xlab, 1, 2.5, cex = cs * cma)
        mtext(ylab, 2, 2.5, cex = cs * cma)
        if (!"lwd" %in% n.args) 
            args$lwd <- 1
        if (!"lty" %in% n.args) 
            args$lty <- 2
        if (!"col" %in% n.args) 
            args$col <- 2
        if (!"labcex" %in% n.args) 
            args$labcex <- cs * 0.5
        zz <- trans(z + zse + shift)
        args$z <- substitute(zz)
        do.call("contour", args)
        if (!titleOnly) {
            xpos <- xrange[1]
            xl <- c(xpos, xpos + xr/10)
            lines(xl, yl, xpd = TRUE, lty = args$lty, col = args$col)
            text(xpos + xr/10, ypos, paste("-", round(se.mult), 
                "se", sep = ""), xpd = TRUE, pos = 4, cex = cs * 
                cm, off = 0.5 * cs * cm)
        }
        if (!"lty" %in% n.args) 
            args$lty <- 3
        if (!"col" %in% n.args) 
            args$col <- 3
        zz <- trans(z - zse + shift)
        args$z <- substitute(zz)
        do.call("contour", args)
        if (!titleOnly) {
            xpos <- xrange[2] - xr/5
            xl <- c(xpos, xpos + xr/10)
            lines(xl, yl, xpd = TRUE, lty = args$lty, col = args$col)
            text(xpos + xr/10, ypos, paste("+", round(se.mult), 
                "se", sep = ""), xpd = TRUE, pos = 4, cex = cs * 
                cm, off = 0.5 * cs * cm)
        }
    }
    w.resid <- NULL
    if (length(residuals) > 1) {
        if (length(residuals) == length(x$residuals)) 
            w.resid <- residuals
        else warning("residuals argument to plot.gam is wrong length: ignored")
        partial.resids <- TRUE
    }
    else partial.resids <- residuals
    m <- length(x$smooth)
    order <- attr(x$pterms, "order")
    if (all.terms) 
        n.para <- sum(order == 1)
    else n.para <- 0
    if (m + n.para == 0) 
        stop("No terms to plot - nothing for plot.gam() to do.")
    if (se) {
        if (is.numeric(se)) 
            se2.mult <- se1.mult <- se
        else {
            se1.mult <- 2
            se2.mult <- 2 ############################### 1
        }
        if (se1.mult < 0) 
            se1.mult <- 0
        if (se2.mult < 0) 
            se2.mult <- 0
    }
    else se1.mult <- se2.mult <- 1
    if (se && x$Vp[1, 1] <= 0) {
        se <- FALSE
        warning("No variance estimates available")
    }
    n.plots <- m + n.para
    if (pages > n.plots) 
        pages <- n.plots
    if (pages < 0) 
        pages <- 0
    if (pages != 0) {
        ppp <- n.plots%/%pages
        if (n.plots%%pages != 0) {
            ppp <- ppp + 1
            while (ppp * (pages - 1) >= n.plots) pages <- pages - 
                1
        }
        c <- trunc(sqrt(ppp))
        if (c < 1) 
            c <- 1
        r <- ppp%/%c
        if (r < 1) 
            r <- 1
        while (r * c < ppp) r <- r + 1
        while (r * c - ppp > c && r > 1) r <- r - 1
        while (r * c - ppp > r && c > 1) c <- c - 1
        oldpar <- par(mfrow = c(r, c))
    }
    else {
        ppp <- 1
        oldpar <- par()
    }
    if ((pages == 0 && prod(par("mfcol")) < n.plots && dev.interactive()) || 
        pages > 1 && dev.interactive()) 
        ask <- TRUE
    else ask <- FALSE
    if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }
    if (partial.resids) {
        fv.terms <- predict(x, type = "terms")
        if (is.null(w.resid)) 
            w.resid <- x$residuals * sqrt(x$weights)
    }
    pd <- list()
    i <- 1
    if (m > 0) 
        for (i in 1:m) {
            if (x$smooth[[i]]$dim == 1) {
                raw <- x$model[x$smooth[[i]]$term]
                xx <- seq(min(raw), max(raw), length = n)
                if (x$smooth[[i]]$by != "NA") {
                  by <- rep(1, n)
                  dat <- data.frame(x = xx, by = by)
                  names(dat) <- c(x$smooth[[i]]$term, x$smooth[[i]]$by)
                }
                else {
                  dat <- data.frame(x = xx)
                  names(dat) <- x$smooth[[i]]$term
                }
                X <- PredictMat(x$smooth[[i]], dat)
                first <- x$smooth[[i]]$first.para
                last <- x$smooth[[i]]$last.para
                p <- x$coefficients[first:last]
                intercept <- 0
### additions for shape-constrained smooths ...
                if (inherits(x$smooth[[i]], c("mpi.smooth","mpd.smooth", 
                   "mdcv.smooth","mdcx.smooth","micv.smooth","micx.smooth"))){
                        q <- ncol(X)
                        beta.c <- c(0,exp(p))
                        fit.c <- X%*%beta.c # fitted values for the mono P-splines identifiability constraints
                          # get an intercept - difference between the fit with monoP-splne constraints and 
                            # fit with the centering constraint...
                        intercept <- -sum(fit.c)/n 
                        onet <- matrix(rep(1,n),1,n)
                        A <- onet%*%X 
                        qrX <- qr(X)
                        R <- qr.R(qrX) # matrix R from QR decomposition
                        qrA <- qr(t(A))
                     #   RZa <- t(qr.qty(qrA,t(R)))[,2:q] 
                     #   RZatRZa.inv <- solve(t(RZa)%*%RZa)
                     #   beta.a <- RZatRZa.inv%*%t(RZa)%*%R%*%beta.c
                        R <- R[-1,]
                        RZa <- t(qr.qty(qrA,t(R)))[,2:q] 
                        RZa.inv <- solve(RZa)
                        RZaR <- RZa.inv%*%R
                        beta.a <- RZaR%*%beta.c
                        p <- c(0,beta.a)
                        p <- qr.qy(qrA,p)
                   #    fit <- X%*%p
                }
                offset <- attr(X, "offset")
                if (is.null(offset)) 
                  fit <- X %*% p - intercept
                else fit <- X %*% p + offset - intercept
                if (se) {
                  if (inherits(x$smooth[[i]], c("mpi.smooth","mpd.smooth", 
                    "mdcv.smooth","mdcx.smooth","micv.smooth","micx.smooth"))){
                        XZa <- t(qr.qty(qrA,t(X)))[,2:q]
                        Ga <- XZa%*%RZaR
                        Vp <- x$Vp.t[c(1,first:last),c(1,first:last)] 
                        Vp.c <- Vp
                        Vp.c[,1] <- rep(0,nrow(Vp))
                        Vp.c[1,] <- rep(0,ncol(Vp))
                        se.fit <- sqrt(rowSums((Ga%*%Vp.c)*Ga))
                  }
                  else {
                      if (seWithMean && inherits(attr(x$smooth[[i]], 
                            "qrc"), "qr")) {
                           X1 <- matrix(x$cmX, nrow(X), ncol(x$Vp), 
                           byrow = TRUE)
                           meanL1 <- x$smooth[[i]]$meanL1
                           if (!is.null(meanL1)) 
                                X1 <- X1/meanL1
                           X1[, first:last] <- X
                           se.fit <- sqrt(rowSums((X1 %*% x$Vp.t) * X1))
                      }
                      else se.fit <- sqrt(rowSums((X %*% x$Vp.t[first:last, 
                        first:last]) * X))
                  }
                }
                edf <- sum(x$edf[first:last])
                xterm <- x$smooth[[i]]$term
                if (is.null(xlab)) 
                  xlabel <- xterm
                else xlabel <- xlab
                if (is.null(ylab)) 
                  ylabel <- sub.edf(x$smooth[[i]]$label, edf)
                else ylabel <- ylab
                pd.item <- list(fit = fit, dim = 1, x = xx, ylab = ylabel, 
                  xlab = xlabel, raw = raw[[1]])
                if (partial.resids) {
                  pd.item$p.resid <- fv.terms[, length(order) + 
                    i] + w.resid
                }
                if (se) {
                          pd.item$se = se.fit * se1.mult
                }
                pd[[i]] <- pd.item
                rm(pd.item)
            }
            else if (x$smooth[[i]]$dim == 2) {
                xterm <- x$smooth[[i]]$term[1]
                if (is.null(xlab)) 
                  xlabel <- xterm
                else xlabel <- xlab
                yterm <- x$smooth[[i]]$term[2]
                if (is.null(ylab)) 
                  ylabel <- yterm
                else ylabel <- ylab
                raw <- data.frame(x = as.numeric(x$model[xterm][[1]]), 
                  y = as.numeric(x$model[yterm][[1]]))
                n2 <- max(10, n2)
                xm <- seq(min(raw$x), max(raw$x), length = n2)
                ym <- seq(min(raw$y), max(raw$y), length = n2)
                xx <- rep(xm, n2)
                yy <- rep(ym, rep(n2, n2))
                if (too.far > 0) 
                  exclude <- exclude.too.far(xx, yy, raw$x, raw$y, 
                    dist = too.far)
                else exclude <- rep(FALSE, n2 * n2)
                if (x$smooth[[i]]$by != "NA") {
                  by <- rep(1, n2^2)
                  dat <- data.frame(x = xx, y = yy, by = by)
                  names(dat) <- c(xterm, yterm, x$smooth[[i]]$by)
                }
                else {
                  dat <- data.frame(x = xx, y = yy)
                  names(dat) <- c(xterm, yterm)
                }
                X <- PredictMat(x$smooth[[i]], dat)
                first <- x$smooth[[i]]$first.para
                last <- x$smooth[[i]]$last.para
                p <- x$coefficients[first:last]
                intercept <- 0
### addition for double and single monotonicity...
                if (inherits(x$smooth[[i]], c("tedmi.smooth","tedmd.smooth", 
                    "tesmi1.smooth","tesmi2.smooth","tesmd1.smooth","tesmd2.smooth")))
                   { p.ident <- x$p.ident[first:last]
                     count<-0 # count the number of exponentiated parameters 
                     for (j in 1:length(p))
                          {if (p.ident[j]==1) count<-count+1}
                     iv<-array(0, dim=c(count,1)) # index vector for the exponent. parameters
                     k<-1
                     for (j in 1:length(p))
                        {if (p.ident[j]==1) {iv[k]<-j; k<-k+1}}
                     p[iv] <- exp(p[iv])
                     beta <- x$smooth[[i]]$Zc%*%p
                     fit.c <- X%*%beta # fitted values for the mono P-splines identifiability constraints
                          # get an fit1intercept - difference between the fit with monoP-splne constraints and 
                            # fit with the centering constraint...
                     intercept <- -sum(fit.c)/(length(fit.c)) 
                     onet <- matrix(rep(1,nrow(X)),1,nrow(X))
                     A <- onet%*%X
                     qrX <- qr(X)
                     R <- qr.R(qrX)
                     qrA <- qr(t(A))
                     q <- ncol(X)
                     if (inherits(x$smooth[[i]], c("tedmi.smooth","tedmd.smooth")))
                     { # get `beta.a' for double monotonicity...
                         R <- R[-1,]
                         RZa <- t(qr.qty(qrA,t(R)))[,2:q] 
                         RZa.inv <- solve(RZa)
                         RZaR <- RZa.inv%*%R
                         beta.a <- RZaR%*%beta
                     }
                     else # get `beta.a' for single monotonicity...
                     { RZa <- t(qr.qty(qrA,t(R)))[,2:q]
                       RZatRZa.inv <- solve(t(RZa)%*%RZa) 
                       Q <- qr.Q(qrX)
                       B1 <- RZatRZa.inv%*%t(RZa)
                       RZaR <- B1%*%R
                       beta.a <- RZaR%*%beta + B1%*%(intercept*colSums(Q))
                     }
                     p <- c(0,beta.a)
                     p <- qr.qy(qrA,p)
                }
                offset <- attr(X, "offset")
                if (is.null(offset)) 
                  fit <- X %*% p - intercept
                else fit <- X %*% p + offset - intercept
                fit[exclude] <- NA
                if (se) {
                   if (inherits(x$smooth[[i]], c("tedmi.smooth","tedmd.smooth", 
                    "tesmi1.smooth","tesmi2.smooth","tesmd1.smooth","tesmd2.smooth"))) {
                          XZa <- t(qr.qty(qrA,t(X)))[,2:ncol(X)]
                          Ga <- XZa%*%RZaR%*%x$smooth[[i]]$Zc
                          Vp <- x$Vp.t[first:last,first:last] 
                          se.fit <- sqrt(rowSums((Ga%*%Vp)*Ga))
                   }
                   else {
                      if (seWithMean && inherits(attr(x$smooth[[i]], 
                           "qrc"), "qr")) {
                          X1 <- matrix(x$cmX, nrow(X), ncol(x$Vp), 
                                  byrow = TRUE)
                          meanL1 <- x$smooth[[i]]$meanL1
                          if (!is.null(meanL1)) 
                                    X1 <- X1/meanL1
                          X1[, first:last] <- X
                          se.fit <- sqrt(rowSums((X1 %*% x$Vp.t) * X1))
                      }
                      else se.fit <- sqrt(rowSums((X %*% x$Vp.t[first:last, 
                                    first:last]) * X))
                   }
                   se.fit[exclude] <- NA
                }
                edf <- sum(x$edf[first:last])
                if (is.null(main)) {
                  title <- sub.edf(x$smooth[[i]]$label, edf)
                }
                else title <- main
                pd.item <- list(fit = fit, dim = 2, xm = xm, 
                  ym = ym, ylab = ylabel, xlab = xlabel, title = title, 
                  raw = raw)
                if (is.null(ylim)) 
                  pd.item$ylim <- range(ym)
                else pd.item$ylim <- ylim
                if (is.null(xlim)) 
                  pd.item$xlim <- range(xm)
                else pd.item$xlim <- xlim
                if (se) { 
                       pd.item$se = se.fit * se2.mult
                }
                pd[[i]] <- pd.item
                rm(pd.item)
            }
            else {
                pd[[i]] <- list(dim = x$smooth[[i]]$dim)
            }
        }
    if (se) {
        k <- 0
        if (scale == -1 && is.null(ylim)) 
            if (m > 0) 
                for (i in 1:m) {
                  if (pd[[i]]$dim == 1) {
                       ul <- pd[[i]]$fit + pd[[i]]$se
                       ll <- pd[[i]]$fit - pd[[i]]$se
                       
                      if (k == 0) {
                      ylim <- c(min(ll), max(ul))
                      k <- 1
                      }
                     else {
                      if (min(ll) < ylim[1]) 
                        ylim[1] <- min(ll)
                      if (max(ul) > ylim[2]) 
                        ylim[2] <- max(ul)
                    }
                    if (partial.resids) {
                      ul <- max(pd[[i]]$p.resid, na.rm = TRUE)
                      if (ul > ylim[2]) 
                        ylim[2] <- ul
                      ll <- min(pd[[i]]$p.resid, na.rm = TRUE)
                      if (ll < ylim[1]) 
                        ylim[1] <- ll
                    }
                  }
                }
        j <- 1
        if (m > 0) 
            for (i in 1:m) {
                if (is.null(select) || i == select) {
                  if (pd[[i]]$dim == 1) {
                       ul <- pd[[i]]$fit + pd[[i]]$se
                       ll <- pd[[i]]$fit - pd[[i]]$se
                   
                    if (scale == 0 && is.null(ylim)) {
                      ylimit <- c(min(ll), max(ul))
                      if (partial.resids) {
                        max.r <- max(pd[[i]]$p.resid, na.rm = TRUE)
                        if (max.r > ylimit[2]) 
                          ylimit[2] <- max.r
                        min.r <- min(pd[[i]]$p.resid, na.rm = TRUE)
                        if (min.r < ylimit[1]) 
                          ylimit[1] <- min.r
                      }
                    }
                    if (!is.null(ylim)) 
                      ylimit <- ylim
                    if (shade) {
                      plot(pd[[i]]$x, trans(pd[[i]]$fit + shift), 
                        type = "n", xlab = pd[[i]]$xlab, ylim = trans(ylimit + 
                          shift), xlim = xlim, ylab = pd[[i]]$ylab, 
                        main = main, ...)
                      polygon(c(pd[[i]]$x, pd[[i]]$x[n:1], pd[[i]]$x[1]), 
                        trans(c(ul, ll[n:1], ul[1]) + shift), 
                        col = shade.col, border = NA)
                      lines(pd[[i]]$x, trans(pd[[i]]$fit + shift))
                    }
                    else {
                      plot(pd[[i]]$x, trans(pd[[i]]$fit + shift), 
                        type = "l", xlab = pd[[i]]$xlab, ylim = trans(ylimit + 
                          shift), xlim = xlim, ylab = pd[[i]]$ylab, 
                        main = main, ...)
                      if (is.null(list(...)[["lty"]])) {
                        lines(pd[[i]]$x, trans(ul + shift), lty = 2, 
                          ...)
                        lines(pd[[i]]$x, trans(ll + shift), lty = 2, 
                          ...)
                      }
                      else {
                        lines(pd[[i]]$x, trans(ul + shift), ...)
                        lines(pd[[i]]$x, trans(ll + shift), ...)
                      }
                    }
                    if (partial.resids && (by.resids || x$smooth[[i]]$by == 
                      "NA")) {
                      if (length(pd[[i]]$raw) == length(pd[[i]]$p.resid)) {
                        if (is.null(list(...)[["pch"]])) 
                          points(pd[[i]]$raw, trans(pd[[i]]$p.resid + 
                            shift), pch = ".", ...)
                        else points(pd[[i]]$raw, trans(pd[[i]]$p.resid + 
                          shift), ...)
                      }
                      else {
                        warning("Partial residuals do not have a natural x-axis location for linear functional terms")
                      }
                    }
                    if (rug) {
                      if (jit) 
                        rug(jitter(as.numeric(pd[[i]]$raw)), 
                          ...)
                      else rug(as.numeric(pd[[i]]$raw), ...)
                    }
                  }
                  else if (pd[[i]]$dim == 2) {
                    if (pers) {
                      if (!is.null(main)) 
                        pd[[i]]$title <- main
                      persp(pd[[i]]$xm, pd[[i]]$ym, matrix(trans(pd[[i]]$fit + 
                        shift), n2, n2), xlab = pd[[i]]$xlab, 
                        ylab = pd[[i]]$ylab, zlab = pd[[i]]$title, 
                        ylim = pd[[i]]$ylim, xlim = pd[[i]]$xlim, 
                        theta = theta, phi = phi, ...)
                    }
                    else {
                      sp.contour(pd[[i]]$xm, pd[[i]]$ym, matrix(pd[[i]]$fit, 
                        n2, n2), matrix(pd[[i]]$se, n2, n2), 
                        xlab = pd[[i]]$xlab, ylab = pd[[i]]$ylab, 
                        zlab = pd[[i]]$title, titleOnly = !is.null(main), 
                        se.mult = se2.mult, trans = trans, shift = shift, 
                        ...)
                      if (rug) {
                        if (is.null(list(...)[["pch"]])) 
                          points(pd[[i]]$raw$x, pd[[i]]$raw$y, 
                            pch = ".", ...)
                        else points(pd[[i]]$raw$x, pd[[i]]$raw$y, 
                          ...)
                      }
                    }
                  }
                  else {
                    warning("no automatic plotting for smooths of more than two variables")
                  }
                }
                j <- j + pd[[i]]$dim
            }
    }
    else {
        k <- 0
        if (scale == -1 && is.null(ylim)) 
            if (m > 0) 
                for (i in 1:m) {
                  if (pd[[i]]$dim == 1) {
                    if (k == 0) {
                      if (partial.resids) 
                        ylim <- range(pd[[i]]$p.resid, na.rm = TRUE)
                      else ylim <- range(pd[[i]]$fit)
                      k <- 1
                    }
                    else {
                      if (partial.resids) {
                        if (min(pd[[i]]$p.resid) < ylim[1]) 
                          ylim[1] <- min(pd[[i]]$p.resid, na.rm = TRUE)
                        if (max(pd[[i]]$p.resid) > ylim[2]) 
                          ylim[2] <- max(pd[[i]]$p.resid, na.rm = TRUE)
                      }
                      else {
                        if (min(pd[[i]]$fit) < ylim[1]) 
                          ylim[1] <- min(pd[[i]]$fit)
                        if (max(pd[[i]]$fit) > ylim[2]) 
                          ylim[2] <- max(pd[[i]]$fit)
                      }
                    }
                  }
                }
        j <- 1
        if (m > 0) 
            for (i in 1:m) {
                if (is.null(select) || i == select) {
                  if (pd[[i]]$dim == 1) {
                    if (scale == 0 && is.null(ylim)) {
                      if (partial.resids) 
                        ylimit <- range(pd[[i]]$p.resid, na.rm = TRUE)
                      else ylimit <- range(pd[[i]]$fit)
                    }
                    if (!is.null(ylim)) 
                      ylimit <- ylim
                    plot(pd[[i]]$x, trans(pd[[i]]$fit + shift), 
                      type = "l", , xlab = pd[[i]]$xlab, ylab = pd[[i]]$ylab, 
                      ylim = trans(ylimit + shift), xlim = xlim, 
                      main = main, ...)
                    if (rug) {
                      if (jit) 
                        rug(jitter(as.numeric(pd[[i]]$raw)), 
                          ...)
                      else rug(as.numeric(pd[[i]]$raw), ...)
                    }
                    if (partial.resids && (by.resids || x$smooth[[i]]$by == 
                      "NA")) {
                      if (is.null(list(...)[["pch"]])) 
                        points(pd[[i]]$raw, trans(pd[[i]]$p.resid + 
                          shift), pch = ".", ...)
                      else points(pd[[i]]$raw, trans(pd[[i]]$p.resid + 
                        shift), ...)
                    }
                  }
                  else if (pd[[i]]$dim == 2) {
                    if (!is.null(main)) 
                      pd[[i]]$title <- main
                    if (pers) {
                      persp(pd[[i]]$xm, pd[[i]]$ym, matrix(trans(pd[[i]]$fit + 
                        shift), n2, n2), xlab = pd[[i]]$xlab, 
                        ylab = pd[[i]]$ylab, zlab = pd[[i]]$title, 
                        theta = theta, phi = phi, xlim = pd[[i]]$xlim, 
                        ylim = pd[[i]]$ylim, ...)
                    }
                    else {
                      contour(pd[[i]]$xm, pd[[i]]$ym, matrix(trans(pd[[i]]$fit + 
                        shift), n2, n2), xlab = pd[[i]]$xlab, 
                        ylab = pd[[i]]$ylab, main = pd[[i]]$title, 
                        xlim = pd[[i]]$xlim, ylim = pd[[i]]$ylim, 
                        ...)
                      if (rug) {
                        if (is.null(list(...)[["pch"]])) 
                          points(pd[[i]]$raw$x, pd[[i]]$raw$y, 
                            pch = ".", ...)
                        else points(pd[[i]]$raw$x, pd[[i]]$raw$y, 
                          ...)
                      }
                    }
                  }
                  else {
                    warning("no automatic plotting for smooths of more than one variable")
                  }
                }
                j <- j + pd[[i]]$dim
            }
    }
    if (n.para > 0) {
        class(x) <- c("gam", "glm", "lm")
        if (is.null(select)) {
            attr(x, "para.only") <- TRUE
            termplot(x, se = se, rug = rug, col.se = 1, col.term = 1)
        }
        else {
            if (select > m) {
                select <- select - m
                term.labels <- attr(x$pterms, "term.labels")
                term.labels <- term.labels[order == 1]
                if (select <= length(term.labels)) {
                  if (interactive() && m && i%%ppp == 0) 
                    termplot(x, terms = term.labels[select], 
                      se = se, rug = rug, col.se = 1, col.term = 1)
                }
            }
        }
    }
    if (pages > 0) 
        par(oldpar)
}


