% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/design_matrix.R
\name{design_matrix}
\alias{design_matrix}
\title{Generic function for creating a design matrix H}
\usage{
design_matrix(x, ...)
}
\arguments{
\item{x}{A basis function object.}

\item{...}{Additional parameters.}
}
\description{
This is a generic function which calles the appropriate methods depending
on the class of the object \code{x}.
}
\examples{
obj <- polynomial.object(M=2)
obs <- c(0,.2,.5)
polyn <- design_matrix(obj, obs)

#----------------

obj <- rbf.object(M=2)
obs <- c(0,.2,.5)
rbf <- design_matrix(obj, obs)

}
\seealso{
\code{\link{polynomial_basis}}, \code{\link{rbf_basis}},
 \code{\link{design_matrix.polynomial}}, \code{\link{design_matrix.rbf}}
}

