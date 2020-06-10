\name{vis.scam}
\alias{vis.scam}
\alias{persp}
\title{Visualization of SCAM objects}
\description{ Produces perspective or contour plot views of \code{scam} model
predictions.
The code is a clone of \code{vis.gam} of the \code{mgcv} package.
 
}
\usage{
vis.scam(x,view=NULL,cond=list(),n.grid=30,too.far=0,col=NA,
        color="heat",contour.col=NULL,se=-1,type="link",
        plot.type="persp",zlim=NULL,nCol=50,...)
}

\arguments{
        The documentation below is the same as in documentation object \code{\link{vis.gam}}.

  \item{x}{a \code{scam} object, produced by \code{scam()}}

  \item{view}{an array containing the names of the two main effect terms to be displayed on the 
   x and y dimensions of the plot. If omitted the first two suitable terms
   will be used.}

  \item{cond}{a named list of the values to use for the other predictor terms
  (not in \code{view}). Variables omitted from this list will have the closest observed value to the median 
  for continuous variables, or the most commonly occuring level for factors. Parametric matrix variables have 
  all the entries in each column set to the observed column entry closest to the column median.
}
  \item{n.grid}{The number of grid nodes in each direction used for calculating the 
   plotted surface.} 

 \item{too.far}{ plot grid nodes that are too far from the points defined by the variables given in \code{view} 
can be excluded from the plot. \code{too.far} determines what is too far. The grid is scaled into the unit 
square along with the \code{view} variables and then grid nodes more than \code{too.far} from the predictor variables 
are excluded.}
\item{col}{The colours for the facets of the plot. If this is \code{NA} then if \code{se}>0 the facets are transparent, 
otherwise the colour scheme specified in \code{color} is used. If \code{col} is not \code{NA} then it is used as the facet 
colour.}
 \item{color}{ the colour scheme to use for plots when \code{se}<=0. One of \code{"topo"}, \code{"heat"}, \code{"cm"}, 
\code{"terrain"}, \code{"gray"} or \code{"bw"}. Schemes \code{"gray"} and
\code{"bw"} also modify the colors used when \code{se}>0.}
\item{contour.col}{sets the colour of contours when using \code{plot.type="contour"}. Default scheme used if \code{NULL}.}
\item{se}{if less than or equal to zero then only the predicted surface is plotted, but if greater than zero, then 3 
surfaces are plotted, one at the predicted values minus \code{se} standard errors, one at the predicted values and one at
the predicted values plus \code{se} standard errors.}
\item{type}{\code{"link"} to plot on linear predictor scale and \code{"response"} to plot on the response scale.}
\item{plot.type}{one of \code{"contour"} or \code{"persp"}.}
\item{zlim}{a two item array giving the lower and upper limits for the z-axis
scale. \code{NULL} to choose automatically.}
\item{nCol}{The number of colors to use in color schemes.}

\item{...}{other options to pass on to \code{\link{persp}},
\code{\link{image}} or \code{\link{contour}}.}


}
\value{Simply produces a plot.}

\author{Simon Wood \email{simon.wood@r-project.org}

}


\seealso{
\code{\link{persp}}, \code{\link{vis.gam}},  and \code{\link{scam}}.
}
\examples{
library(scam)

# Example with factor variable
set.seed(0)
fac<-rep(1:4,20)
x <- runif(80)*5;
y <- fac+log(x)/5+rnorm(80)*0.1
fac <- factor(fac)
b <- scam(y~fac+s(x,bs="mpi"))

vis.scam(b,theta=-35,color="heat") # factor example

# Example with "by" variable

z<-rnorm(80)*0.4   
y<-as.numeric(fac)+log(x)*z+rnorm(80)*0.1
b<-scam(y~fac+s(x,by=z))
g <- gam(y~fac+s(x,by=z))

vis.scam(b,theta=-35,color="terrain",cond=list(z=1)) # by variable example
vis.scam(b,view=c("z","x"),theta= 65) # plot against by variable
## compare with gam(mgcv)...
vis.gam(g,theta=-35,color="terrain",cond=list(z=1)) # by variable example
vis.gam(g,view=c("z","x"),theta= 65) # plot against by variable

## all three smooths are increasing...
set.seed(2)
n <- 400
x0 <- runif(n, 0, 1)
x <- runif(n, 0, 1)
f1 <- log(x *5)
f2 <-  exp(2 * x) - 4
f3 <-  5* sin(x)
e <- rnorm(n, 0, 2)
fac <- as.factor(sample(1:3,n,replace=TRUE))
fac.1 <- as.numeric(fac==1)
fac.2 <- as.numeric(fac==2)
fac.3 <- as.numeric(fac==3)
y <- f1*fac.1 + f2*fac.2 + f3*fac.3 + e 
dat <- data.frame(y=y,x=x,fac=fac,f1=f1,f2=f2,f3=f3)

b1 <- scam(y ~ s(x,by=fac,bs="mpi"),data=dat,optimizer="efs")
plot(b1,pages=1,scale=0)
summary(b1)
vis.scam(b1,theta=-40,color="terrain",cond=list(z=1))

## note that the preceding, b1, fit is the same as....
b2 <- scam(y ~ s(x,by=as.numeric(fac==1),bs="mpi")+s(x,by=as.numeric(fac==2),bs="mpi")+
    s(x,by=as.numeric(fac==3),bs="mpi"),data=dat,optimizer="efs")
summary(b2)

## Also note that as in gam() when using factor 'by' variables, centering
## constraints are applied to the smooths, which usually means that the 'by'
## variable should be included as a parametric term, as well. That is not 
## the case with scam(), since here a 'zero intercept' constraint is 
## applied (although scam's fitted smooths are centred for plotting).
g1 <- gam(y ~ fac+s(x,by=fac),data=dat)
g2 <- gam(y ~ s(x,by=fac),data=dat)
summary(g1)
summary(g2)
}
\keyword{hplot} \keyword{models} \keyword{smooth} \keyword{regression}
