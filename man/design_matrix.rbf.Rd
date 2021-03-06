% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/design_matrix.R
\name{design_matrix.rbf}
\alias{design_matrix.rbf}
\title{Creates an RBF design matrix H}
\usage{
\method{design_matrix}{rbf}(x, obs, ...)
}
\arguments{
\item{x}{A basis object.}

\item{obs}{A vector of observations.}

\item{...}{Additional parameters}
}
\value{
A list containing the design matrix \code{H} and the basis object.
The dimensions of the matrix H are Nx(M+1), where N is the length of the
observations, and M is the number of radial basis functions. The updated
\code{basis} object contains also the updated centers of RBFs.
}
\description{
\code{design_matrix.rbf} creates a design matrix H using radial
basis functions of degree M.
}
\examples{
obj <- rbf.object(M=3)
obs <- c(0,.2,.5, 0.3)
des_mat <- design_matrix(obj, obs)

}
\seealso{
\code{\link{design_matrix}}, \code{\link{rbf_basis}}
}

