\name{paRaoP}
\alias{paRaoP}
\title{
Parallelized Parametric Rao's index of quadratic entropy (Q)
}
\description{
Parametric Rao's index of quadratic entropy (Q).
}
\usage{
paRaoS(rasterm,alpha,w,dist_m,na.tolerance,
	diag,debugging,isfloat,mfactor)
}
\arguments{
\item{rasterm}{input data.}
\item{alpha}{alpha value for order of diversity in Hill's Index.}
\item{w}{half of the side of the square moving window.}
\item{dist_m}{Type of distance used.}
\item{na.tolerance}{a numeric value \eqn{(0.0-1.0)} which indicates the proportion of NA values that will be tolerated to calculate Rao's index in each moving window over \emph{x}. If the relative proportion of NA's in a moving window is bigger than na.tolerance, then the value of the window will be set as NA, otherwise Rao's index will be calculated considering the non-NA values. Default values is 0.0 (i.e., no tolerance for NA's).}
\item{diag}{Boolean. Diagonal of the distance matrix.}
\item{debugging}{a boolean variable set to FALSE by default. If TRUE, additional messages will be printed. For de-bugging only.}
\item{isfloat}{Are the input data floats?}
\item{debugging}{Debugging mode, TRUE or FALSE.}
}
\value{
A list of matrices of dimension \code{dim(x)} with length equal to the length of \code{alpha}.
}
\author{
Duccio Rocchini \email{duccio.rocchini@unibo.it} \cr
Marcantonio Matteo \email{marcantoniomatteo@gmail.com} \cr
}
\seealso{
\code{\link{paRao}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}
