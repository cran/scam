
## printing the results of the scam...

print.scam <- function (x,...) 
   {
    print(x$family)
    cat("Formula:\n")
    print(x$formula)
    n.smooth <- x$n.smooth
    if (n.smooth == 0) 
        cat("Total model degrees of freedom", sum(x$edf), "\n")
    else {
        edf <- 0
        cat("\nEstimated degrees of freedom:\n")
        for (i in 1:n.smooth) edf[i] <-  sum(x$edf[x$smooth[[i]]$first.para:x$smooth[[i]]$last.para])
        edf.str <- format(edf, digits = 5)
        for (i in 1:n.smooth) {
            cat(edf.str[i], " ", sep = "")
            if (i%%7 == 0) 
                cat("\n")
        }
        cat(" total =", x$trA, "\n")
    }
    cat("\n",x$method," score: ", x$gcv.ubre, "\n", sep = "")
#    dgcv <- 0
#    cat("\nGCV gradient:\n")
#    for (i in 1:length(x$sp)) dgcv[i] <- x$gcv.g[i]
#    dgcv.str <- format(dgcv, digits = 5)
#    for (i in 1:length(x$sp)) {
#       cat(dgcv.str[i], " ", sep = "")
#       if (i%%7 == 0) 
#           cat("\n")
#    }
    cat("\n")
    invisible(x)
}

