#' Illumina 450K experimental data
#'
#' A dataset containing beta values for 2 experimental condition in triplicats.
#'
#' @format A matrix with beta values:
#' \describe{
#'   \item{Rows}{CpG}
#'   \item{Columns}{Samples}
#'   ...
#' }
#' @source \url{http://www.diamondse.info/}
"betas"


#' Illumina 450K experimental setup
#'
#' A dataset containing beta values for 2 experimental condition in triplicats.
#'
#' @format A dataframe describing samples as in the Minfi Package:
#' \describe{
#'    \item{Sample_Name}{names of samples}
#'    \item{Sample_Group}{Experimental condition, must be a factor with first level corresponding to the control}
#'    \item{Pool_ID}{Index for paired analysis}
#'    }
#' @source \url{http://www.diamondse.info/}
"expSamples"


