% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ShannonP.R
\name{ShannonP}
\alias{ShannonP}
\title{Calculate Shannon-Wiener Index on a Matrix}
\usage{
ShannonP(x, window = 1, na.tolerance = 1, debugging = FALSE, np = 1)
}
\arguments{
\item{x}{A numeric matrix representing the data on which the index is to abe calculated.}

\item{window}{The width of the moving window to consider for each cell. 
The actual window size will be `(2 * window + 1) x (2 * window + 1)`. Default is 1.}

\item{na.tolerance}{The tolerance level for missing data within the moving window. 
A window will be processed only if the proportion of non-missing data is above this threshold. 
Value should be between 0 and 1. Default is 1.}

\item{debugging}{Boolean flag to enable or disable debugging messages. Default is FALSE.}

\item{np}{Number of processes for parallel computation.}
}
\value{
A matrix of the same dimensions as `x`, where each cell contains the 
        Shannon-Wiener Index calculated for the window around the cell.
}
\description{
This function computes Shannon-Wiener Index for each cell of a matrix, 
using a parallelized approach and considering a specified moving window.
}
\examples{
data <- matrix(runif(100), nrow = 10)
shannon_index <- ShannonP(data, window = 1, np = 1 )

}
