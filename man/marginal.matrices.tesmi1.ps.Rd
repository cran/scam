\name{marginal.matrices.tesmi1.ps}
% \Rdversion{1.0}
\alias{marginal.matrices.tesmi1.ps}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Constructs marginal model matrices for "tesmi1" and "tesmd1" bivariate smooths in case of B-splines basis functions 
      for both unconstrained marginal smooths 
}
\description{This function returns the marginal model matrices and the list of penalty matrices
for the tensor product bivariate smooth with the single monotone increasing or decreasing restriction along the first covariate. The 
marginal smooth functions of both covariates are constructed using the B-spline basis functions. 
}
\usage{
marginal.matrices.tesmi1.ps(x, z, xk, zk, m, q1, q2)
}
\arguments{
  \item{x}{A numeric vector of the values of the first covariate at which to evaluate the B-spline marginal functions.
    The values in \code{x} must be between \code{xk[m[1]+2]} and
          \code{xk[length(xk) - m[1] - 1]}.}

  \item{z}{A numeric vector of the values of the second covariate at which to evaluate the B-spline marginal functions.
    The values in \code{z} must be between \code{zk[m[2]+2]} and
          \code{zk[length(zk) - m[2] - 1]}.}

  \item{xk}{A numeric vector of knot positions for the first covariate, \code{x},  with non-decreasing values.}

  \item{zk}{A numeric vector of knot positions for the second covariate,\code{z}, with non-decreasing values.}

  \item{m}{A pair of two numbers where \code{m[i]+1} denotes the order of the basis of the \eqn{i^{th}}{i^{th}} marginal smooth
   (e.g. \code{m[i] = 2} for a cubic spline.)} 

  \item{q1}{A number denoting the basis dimension of the first marginal smooth.}

  \item{q2}{A number denoting the basis dimension of the second marginal smooth.}

}
\details{The function is not called directly, but is rather used internally by 
the constructor 

\code{\link{smooth.construct.tesmi1.smooth.spec}} and \code{\link{smooth.construct.tesmd1.smooth.spec}} . 

%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{X1}{Marginal model matrix for the first monotonic marginal smooth.}
  \item{X2}{Marginal model matrix for the second unconstrained marginal smooth.}
  \item{S}{A list of penalty matrices for this tensor product smooth.}
}
\references{

Pya, N. and Wood, S.N. (2014) Shape constrained additive models. Statistics and Computing, DOI \url{10.1007/s11222-013-9448-7}

Pya, N. (2010) Additive models with shape constraints. PhD thesis. University of Bath. Department of Mathematical Sciences

Wood S.N. (2006) Generalized Additive Models: An Introduction with R. Chapman and Hall/CRC Press.

}
\author{
Natalya Pya <nat.pya@gmail.com>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{smooth.construct.tesmi1.smooth.spec}}, 
\code{\link{smooth.construct.tesmi2.smooth.spec}}, 
\code{\link{marginal.matrices.tesmi2.ps}},
\code{\link{smooth.construct.tesmd1.smooth.spec}}, 
\code{\link{smooth.construct.tesmd2.smooth.spec}}
}
\keyword{models} \keyword{regression} %-- one or more ..

