% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bpr_optimize.R
\name{bpr_optim.list}
\alias{bpr_optim.list}
\title{Optimization method for the BPR NLL function using list x}
\usage{
\method{bpr_optim}{list}(x, w = NULL, basis = NULL, fit_feature = NULL,
  cpg_dens_feat = FALSE, opt_method = "CG", opt_itnmax = 100,
  is_parallel = TRUE, no_cores = NULL, ...)
}
\arguments{
\item{x}{A list of elements of length N, where each element is an L x 3
matrix of observations, where 1st column contains the locations. The 2nd
and 3rd columns contain the total trials and number of successes at the
corresponding locations, repsectively.}

\item{w}{A vector of parameters (i.e. coefficients of the basis functions)}

\item{basis}{A 'basis' object. See \code{\link{polynomial.object}}}

\item{fit_feature}{Additional feature on how well the profile fits the
methylation data.}

\item{cpg_dens_feat}{Additional feature for the CpG density across the
promoter region.}

\item{opt_method}{The optimization method to be used. See
\code{\link[stats]{optim}} for possible methods. Default is 'CG'.}

\item{opt_itnmax}{Optional argument giving the maximum number of iterations
for the corresponding method. See \code{\link[stats]{optim}} for details.}

\item{is_parallel}{Logical, indicating if code should be run in parallel.}

\item{no_cores}{Number of cores to be used, default is max_no_cores - 1.}

\item{...}{Additional parameters}
}
\value{
A list containing the following elements:
\itemize{
 \item{ \code{W_opt}: An Nx(M+1) matrix with the optimized parameter values.
   Each row of the matrix corresponds to each element of the list x. The
   columns are of the same length as the parameter vector w (i.e. number
   of basis functions).
 }
 \item{ \code{Mus}: An N x M matrix with the RBF centers if basis object is
  \code{\link{rbf.object}}, otherwise NULL.}
 \item{ \code{basis}: The basis object.
 }
 \item{ \code{w}: The initial values of the parameters w.
 }
 \item{ \code{x_extrema}: The min and max values of each promoter region.
 }
}
}
\description{
\code{bpr_optim.list} minimizes the negative log likelihood of the BPR
function. Since it cannot be evaluated analytically, an optimization
procedure is used. The \code{\link[stats]{optim}} packages is used for
performing optimization. This method calls \code{\link{bpr_optim.matrix}}
to process each element of the list.
}
\examples{
ex_data <- bpr_data
basis <- rbf.object(M=3)
out_opt <- bpr_optim(x = ex_data, is_parallel = TRUE, basis = basis, opt_method = "CG")

}
\seealso{
\code{\link{bpr_optim}}, \code{\link{bpr_optim.matrix}}
}

