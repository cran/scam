###################################################################
## printing the results of the scam (clone of print.gam())...    ##
###################################################################

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
        for (i in 1:n.smooth) 
           edf[i] <-  sum(x$edf[x$smooth[[i]]$first.para:x$smooth[[i]]$last.para])
        edf.str <- format(round(edf,digits=4),digits=3,scientific=FALSE) 
        for (i in 1:n.smooth) {
            cat(edf.str[i], " ", sep = "")
            if (i%%7 == 0) 
                cat("\n")
        }
        cat(" total =",round(sum(x$edf),digits=2),"\n")    ##  (" total =", x$trA, "\n")
    }
    if (!is.null(x$gcv.ubre))
         cat("\n",x$method," score: ", formatC(x$gcv.ubre, digits = 5), "\n", sep = "")
    if (!is.null(x$rank) && x$rank< length(x$coefficients)) cat("rank: ",x$rank,"/",length(x$coefficients),sep="") 
    cat("\n")
    invisible(x)
}

