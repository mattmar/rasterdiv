% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/paRaoP.R
\name{paRaoP}
\alias{paRaoP}
\title{Parallelized Parametric Rao's index of quadratic entropy (Q)}
\usage{
paRaoP(
  x,
  alpha,
  window,
  dist_m,
  na.tolerance,
  diag,
  debugging,
  isfloat,
  mfactor,
  np,
  progBar
)
}
\arguments{
\item{x}{Matrix or data frame; the input data over which the index calculation is performed.}

\item{alpha}{Numeric; specifies the alpha value for the order of diversity in Hill's Index.}

\item{window}{Numeric; half of the side length of the square moving window used in the calculation.}

\item{dist_m}{Character; specifies the type of distance metric used in calculations.}

\item{na.tolerance}{Numeric; the threshold proportion of NA values allowed in the moving window. 
If exceeded, the calculation for that window is skipped. Values range from 0.0 (no tolerance) to 1.0.}

\item{diag}{Logical; indicates whether the diagonal of the distance matrix should be included in the 
computation. Typically set to FALSE.}

\item{debugging}{Logical; set to FALSE by default. If TRUE, additional console messages will be 
displayed for debugging purposes.}

\item{isfloat}{Logical; indicates whether the input data values are floating-point numbers.}

\item{mfactor}{Integer; multiplication factor in case of input data as float numbers.}

\item{np}{Number of processes for parallel computation.}

\item{progBar}{logical. If TRUE a progress bar is shown.}
}
\value{
A list of matrices corresponding to the computed Rao's index values. Each matrix in the list 
represents the calculations performed over the moving window, with dimensions equal to \code{dim(x)}.
}
\description{
This function computes the parametric Rao's index of quadratic entropy (Q), a measure of biodiversity 
that considers the evolutionary distances between species, utilizing parallel computing for enhanced 
performance. The computation is applied over a moving window across the input data.
}
\seealso{
\code{\link{paRao}} for the related non-parallelized function.
}
\author{
Duccio Rocchini \email{duccio.rocchini@unibo.it},
Matteo Marcantonio \email{marcantoniomatteo@gmail.com}
}
