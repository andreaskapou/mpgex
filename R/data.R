#' Synthetic data for mpgex package
#'
#' A synthetic dataset containinig 600 entries.
#'
#' @format A list with 600 elements, where each element element is an L x 3
#'  matrix of observations, where:
#' \describe{
#'   \item{1st column}{locations of obseravtions}
#'   \item{2nd column}{total trials at corresponding locations}
#'   \item{3rd column}{number of successes at corresponding locations}
#' }
"bpr_data"

#' Synthetic control data for mpgex package
#'
#' A synthetic dataset containinig 600 entries for control data.
#'
#' @format A list with 600 elements, where each element element is an L x 3
#'  matrix of observations, where:
#' \describe{
#'   \item{1st column}{locations of obseravtions}
#'   \item{2nd column}{total trials at corresponding locations}
#'   \item{3rd column}{number of successes at corresponding locations}
#' }
"bpr_control_data"

#' Synthetic treatment data for mpgex package
#'
#' A synthetic dataset containinig 600 entries for treatment data.
#' @format A list with 600 elements, where each element element is an L x 3
#'  matrix of observations, where:
#' \describe{
#'   \item{1st column}{locations of obseravtions}
#'   \item{2nd column}{total trials at corresponding locations}
#'   \item{3rd column}{number of successes at corresponding locations}
#' }
"bpr_treatment_data"


#' Synthetic data for mpgex package
#'
#' Corresponding gene expression data for the 'bpr_data'
#'
#' @format A vector of length 600
"gex_data"


#' Synthetic control data for mpgex package
#'
#' Corresponding gene expression data for the 'bpr_control_data'
#'
#' @format A vector of length 600
"gex_control_data"


#' Synthetic treatment data for mpgex package
#'
#' Corresponding gene expression data for the 'bpr_treatment_data'
#'
#' @format A vector of length 600
"gex_treatment_data"
