
## summary function for scam() (clone of summary.gam())...
##### function to get all the summary information....

summary.scam <- function (object,dispersion = NULL,freq = FALSE,alpha = 0,...) 
{
    pinv <- function(V, M, rank.tol = 1e-06) {
        D <- La.svd(V)
        M1 <- length(D$d[D$d > rank.tol * D$d[1]])
        if (M > M1) 
            M <- M1
        if (M + 1 <= length(D$d)) 
            D$d[(M + 1):length(D$d)] <- 1
        D$d <- 1/D$d
        if (M + 1 <= length(D$d)) 
            D$d[(M + 1):length(D$d)] <- 0
        res <- D$u %*% diag(D$d) %*% D$v
        attr(res, "rank") <- M
        res
    }
    p.table <- pTerms.table <- s.table <- NULL
    if (freq) 
        covmat <- object$Ve.t
    else covmat <- object$Vp.t
    name <- names(object$coefficients.t)
    dimnames(covmat) <- list(name, name)
    covmat.unscaled <- covmat/object$sig2
    est.disp <- TRUE
    if (object$method == "UBRE") 
        est.disp <- FALSE
    if (!is.null(dispersion)) {
        covmat <- dispersion * covmat.unscaled
        est.disp <- FALSE
    }
    else dispersion <- object$sig2

    se <- 0
    for (i in 1:length(object$coefficients)) se[i] <- covmat[i, 
        i]^0.5
    residual.df <- length(object$y) - object$trA
    if (object$nsdf > 0) {
        p.coeff <- object$coefficients.t[1:object$nsdf]
        p.se <- se[1:object$nsdf]
        p.t <- p.coeff/p.se
        if (!est.disp) {
            p.pv <- 2 * pnorm(abs(p.t), lower.tail = FALSE)
            p.table <- cbind(p.coeff, p.se, p.t, p.pv)
            dimnames(p.table) <- list(names(p.coeff), c("Estimate", 
                "Std. Error", "z value", "Pr(>|z|)"))
        }
        else {
            p.pv <- 2 * pt(abs(p.t), df = residual.df, lower.tail = FALSE)
            p.table <- cbind(p.coeff, p.se, p.t, p.pv)
            dimnames(p.table) <- list(names(p.coeff), c("Estimate", 
                "Std. Error", "t value", "Pr(>|t|)"))
        }
    }
    else {
        p.coeff <- p.t <- p.pv <- array(0, 0)
    }
    term.labels <- attr(object$pterms, "term.labels")
    nt <- length(term.labels)
    if (nt > 0) {
        np <- length(object$assign)
        Vb <- matrix(covmat[1:np, 1:np], np, np)
        bp <- array(object$coefficients[1:np], np)
        pTerms.pv <- array(0, nt)
        attr(pTerms.pv, "names") <- term.labels
        pTerms.df <- pTerms.chi.sq <- pTerms.pv
        for (i in 1:nt) {
            ind <- object$assign == i
            b <- bp[ind]
            V <- Vb[ind, ind]
            if (length(b) == 1) {
                V <- 1/V
                pTerms.df[i] <- nb <- 1
                pTerms.chi.sq[i] <- V * b * b
            }
            else {
                V <- pinv(V, length(b), rank.tol = .Machine$double.eps^0.5)
                pTerms.df[i] <- nb <- attr(V, "rank")
                pTerms.chi.sq[i] <- t(b) %*% V %*% b
            }
            if (!est.disp) 
                pTerms.pv[i] <- pchisq(pTerms.chi.sq[i], df = nb, 
                  lower.tail = FALSE)
            else pTerms.pv[i] <- pf(pTerms.chi.sq[i]/nb, df1 = nb, 
                df2 = residual.df, lower.tail = FALSE)
        }
        if (!est.disp) {
            pTerms.table <- cbind(pTerms.df, pTerms.chi.sq, pTerms.pv)
            dimnames(pTerms.table) <- list(term.labels, c("df", 
                "Chi.sq", "p-value"))
        }
        else {
            pTerms.table <- cbind(pTerms.df, pTerms.chi.sq/pTerms.df, 
                pTerms.pv)
            dimnames(pTerms.table) <- list(term.labels, c("df", 
                "F", "p-value"))
        }
    }
    else {
        pTerms.df <- pTerms.chi.sq <- pTerms.pv <- array(0, 0)
    }
    m <- length(object$smooth)
    df <- edf <- s.pv <- chi.sq <- array(0, m)
    if (m > 0) {
        if (!freq) {
            X <- object$X
            X <- X[!is.na(rowSums(X)), ]
        }
        for (i in 1:m) {
            start <- object$smooth[[i]]$first.para
            stop <- object$smooth[[i]]$last.para
            V <- covmat[start:stop, start:stop] # covmat of gamma (transposed parameters)
            p <- object$coefficients.t[start:stop] # transposed parameters of a smooth
            edf[i] <- sum(object$edf[start:stop])
            if (freq) {
                M1 <- object$smooth[[i]]$df
                M <- min(M1, ceiling(2 * sum(object$edf[start:stop])))
                V <- pinv(V, M)
                chi.sq[i] <- t(p) %*% V %*% p
                df[i] <- attr(V, "rank")
            }
            else {
                Xt <- X[, start:stop]
                ft <- Xt %*% p
                if (FALSE) {
                  trial.rank <- ceiling(edf[i])
                  if (edf[i] - trial.rank > 0) 
                    trial.rank <- trial.rank + 1
                  ed <- eigXVX(Xt, V, trial.rank)
                  iv <- 1/ed$values
                  chi.sq[i] <- sum(((t(ed$vectors) %*% ft) * 
                    sqrt(iv))^2)
                  df[i] <- edf[i] + 0.5
                }
                else {
                  df[i] <- min(ncol(Xt), edf[i])
                  D <- pinvXVX(Xt, V, df[i])
                  df[i] <- df[i] + alpha * sum(object$smooth[[i]]$sp < 
                    0)
                  chi.sq[i] <- sum((t(D) %*% ft)^2)
                  if (FALSE) {
                    if (inherits(attr(object$smooth[[i]], "qrc"), 
                      "qr")) {
                      X1 <- matrix(object$cmX, nrow(Xt), ncol(object$Vp.t), 
                        byrow = TRUE)
                      meanL1 <- object$smooth[[i]]$meanL1
                      if (!is.null(meanL1)) 
                        X1 <- X1/meanL1
                      X1[, start:stop] <- Xt
                      df[i] <- edf[i]
                      D <- pinvXVX(X1, object$Vp.t, df[i] + 1)
                    }
                    else {
                      df[i] <- edf[i]
                      D <- pinvXVX(Xt, V, df[i])
                    }
                    fm <- sum(D %*% (t(D) %*% ft))/sum(colSums(D)^2)
                    chi.sq[i] <- sum((t(D) %*% (ft - fm))^2)
                  }
                }
            }
            names(chi.sq)[i] <- object$smooth[[i]]$label
            if (!est.disp) 
                s.pv[i] <- pchisq(chi.sq[i], df = df[i], lower.tail = FALSE)
            else s.pv[i] <- pf(chi.sq[i]/df[i], df1 = df[i], 
                df2 = residual.df, lower.tail = FALSE)
            if (df[i] < 0.5) 
                s.pv[i] <- NA
        }
        if (!est.disp) {
            if (freq) {
                s.table <- cbind(edf, df, chi.sq, s.pv)
                dimnames(s.table) <- list(names(chi.sq), c("edf", 
                  "Est.rank", "Chi.sq", "p-value"))
            }
            else {
                s.table <- cbind(edf, df, chi.sq, s.pv)
                dimnames(s.table) <- list(names(chi.sq), c("edf", 
                  "Ref.df", "Chi.sq", "p-value"))
            }
        }
        else {
            if (freq) {
                s.table <- cbind(edf, df, chi.sq/df, s.pv)
                dimnames(s.table) <- list(names(chi.sq), c("edf", 
                  "Est.rank", "F", "p-value"))
            }
            else {
                s.table <- cbind(edf, df, chi.sq/df, s.pv)
                dimnames(s.table) <- list(names(chi.sq), c("edf", 
                  "Ref.df", "F", "p-value"))
            }
        }
    }
    w <- object$prior.weights
    nobs <- nrow(object$X)
    r.sq <- 1 - sum(w^2 * (object$y - object$fitted.values)^2)/
                   (var(w * object$y) * residual.df)
    dev.expl <- (object$null.deviance - object$deviance)/object$null.deviance
    ret <- list(p.coeff = p.coeff, se = se, p.t = p.t, p.pv = p.pv, 
        residual.df = residual.df, m = m, chi.sq = chi.sq, s.pv = s.pv, 
        scale = dispersion, r.sq = r.sq, family = object$family, 
        formula = object$formula, n = nobs, dev.expl = dev.expl, 
        edf = edf, dispersion = dispersion, pTerms.pv = pTerms.pv, 
        pTerms.chi.sq = pTerms.chi.sq, pTerms.df = pTerms.df, 
        cov.unscaled = covmat.unscaled, cov.scaled = covmat, 
        p.table = p.table, pTerms.table = pTerms.table, s.table = s.table, 
        method = object$method,sp.criterion = object$gcv.ubre,sp=object$sp,
        dgcv.ubre=object$dgcv.ubre,termcode=object$termcode,
         gcv.ubre=object$gcv.ubre,optimizer=object$optimizer)
    class(ret) <- "summary.scam"
    ret
}


################## mgcv::: pinvXVX

pinvXVX <- function (X, V, rank = NULL) 
{
    k <- floor(rank)
    nu <- rank - k
    if (nu > 0) 
        k1 <- k + 1
    else k1 <- k
    qrx <- qr(X)
    R <- qr.R(qrx)
    V <- R %*% V[qrx$pivot, qrx$pivot] %*% t(R)
    V <- (V + t(V))/2
    ed <- eigen(V, symmetric = TRUE)
    vec <- qr.qy(qrx, rbind(ed$vectors, matrix(0, nrow(X) - ncol(X), 
        ncol(X))))
    if (k1 < ncol(vec)) 
        vec <- vec[, 1:k1, drop = FALSE]
    if (k == 0) {
        vec <- t(t(vec) * sqrt(nu/ed$val[1]))
        return(vec)
    }
    if (nu > 0) {
        if (k > 1) 
            vec[, 1:(k - 1)] <- t(t(vec[, 1:(k - 1)])/sqrt(ed$val[1:(k - 
                1)]))
        b12 <- 0.5 * nu * (1 - nu)
        if (b12 < 0) 
            b12 <- 0
        b12 <- sqrt(b12)
        B <- matrix(c(1, b12, b12, nu), 2, 2)
        ev <- diag(ed$values[k:k1]^-0.5)
        B <- ev %*% B %*% ev
        eb <- eigen(B, symmetric = TRUE)
        rB <- eb$vectors %*% diag(sqrt(eb$values)) %*% t(eb$vectors)
        vec[, k:k1] <- t(rB %*% t(vec[, k:k1]))
    }
    else {
        vec <- t(t(vec)/sqrt(ed$val[1:k]))
    }
    vec
}
 
################ mgcv::: eigXVX

eigXVX <- function (X, V, rank = NULL, tol = .Machine$double.eps^0.5) 
{
    qrx <- qr(X)
    R <- qr.R(qrx)
    V <- R %*% V[qrx$pivot, qrx$pivot] %*% t(R)
    V <- (V + t(V))/2
    ed <- eigen(V, symmetric = TRUE)
    ind <- abs(ed$values) > max(abs(ed$values)) * tol
    erank <- sum(ind)
    if (is.null(rank)) {
        rank <- erank
    }
    else {
        if (rank < erank) 
            ind <- 1:rank
        else rank <- erank
    }
    vec <- qr.qy(qrx, rbind(ed$vectors, matrix(0, nrow(X) - ncol(X), 
        ncol(X))))
    list(values = ed$values[ind], vectors = vec[, ind], rank = rank)
}


##### print.summary.scam .....
print.summary.scam <- function (x, digits = max(3, getOption("digits") - 3), signif.stars = getOption("show.signif.stars"), 
    ...) 
{
    print(x$family)
    cat("Formula:\n")
    print(x$formula)
    if (length(x$p.coeff) > 0) {
        cat("\nParametric coefficients:\n")
        printCoefmat(x$p.table, digits = digits, signif.stars = signif.stars, 
            na.print = "NA", ...)
    }
    cat("\n")
    if (x$m > 0) {
        cat("Approximate significance of smooth terms:\n")
        printCoefmat(x$s.table, digits = digits, signif.stars = signif.stars, 
            has.Pvalue = TRUE, na.print = "NA", cs.ind = 1, ...)
    }
    cat("\nR-sq.(adj) = ", formatC(x$r.sq, digits = 3, width = 5))
    if (length(x$dev.expl) > 0) 
        cat("   Deviance explained = ", formatC(x$dev.expl * 
            100, digits = 3, width = 4), "%\n", sep = "")
    cat( x$method," score = ", formatC(x$sp.criterion, digits = 5), 
            sep = "")
    cat("  Scale est. = ", formatC(x$scale, digits = 5, width = 8, 
        flag = "-"), "  n = ", x$n, "\n", sep = "")
    if (x$optimizer == "bfgs"){
               if (x$termcode!= 1) {
                   dgcv.ubre <- max(abs(x$dgcv.ubre)*max(abs(log(x$sp)),1)/max(abs(x$gcv.ubre),1))
                  cat("\nBFGS termination condition:\n", dgcv.ubre,"\n",sep = "")
               }   
    }
    cat("\n")
    invisible(x)
}



