
print.scam.version <- function()
{ library(help=scam)$info[[1]] -> version
  version <- version[pmatch("Version",version)]
  um <- strsplit(version," ")[[1]]
  version <- um[nchar(um)>0][2]
  hello <- paste("This is scam ",version,".",sep="")
  packageStartupMessage(hello)
}

.onAttach <- function(...) { 
  print.scam.version()
}

.onLoad <- function(libname, pkgname) {
    if (requireNamespace("emmeans", quietly = TRUE))
        emmeans::.emm_register(c("scam"), pkgname)
}



##.onUnload <- function(libpath) library.dynam.unload("scam", libpath)


