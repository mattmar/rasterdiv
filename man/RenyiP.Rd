% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RenyiP.R
\name{RenyiP}
\alias{RenyiP}
\title{Parallel Computation of Renyi's Diversity Index}
\usage{
RenyiP(
  x,
  window = 1,
  alpha = 1,
  base = exp(1),
  na.tolerance = 1,
  debugging = FALSE,
  np = 1
)
}
\arguments{
\item{x}{A numeric matrix representing the data on which the index is to be calculated.}

\item{window}{The width of the moving window to consider for each cell. 
The actual window size will be `(2 * window + 1) x (2 * window + 1)`. Default is 1.}

\item{alpha}{The alpha parameter for Renyi's index, influencing sensitivity 
to species abundance. Default is 1.}

\item{base}{The base of the logarithm used in Renyi's formula. Default is `exp(1)` 
(natural logarithm).}

\item{na.tolerance}{The tolerance level for missing data within the moving window. 
A window will be processed only if the proportion of non-missing data is above this threshold. 
Value should be between 0 and 1. Default is 1.}

\item{debugging}{Boolean flag to enable or disable debugging messages. Default is FALSE.}

\item{np}{Number of processes for parallel computation.#'}
}
\value{
A matrix of the same dimensions as `x`, where each cell contains the 
        Renyi's diversity index calculated for the window around the cell.
}
\description{
This function computes Renyi's diversity index for each cell of a matrix, 
using a parallelized approach and considering a specified moving window.
}
\examples{
data <- matrix(runif(100), nrow = 10)
renyi_index <- RenyiP(data, window = 1, np = 1)

}
