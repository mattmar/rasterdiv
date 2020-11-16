\name{RenyiP}
\alias{RenyiP}
\title{
Parallelised Renyi's diversity index
}
\description{
Parallelised Renyi's diversity index.
}
\usage{
RenyiP(rasterm, w, alpha, base, na.tolerance, debugging)
}
\arguments{
\item{rasterm}{input data.}
\item{w}{half of the side of the square moving window.}
\item{alpha}{alpha value for order of diversity in Renyi's Index.}
\item{base}{base of the logarithm to calculate the index.}
\item{na.tolerance}{a numeric value \eqn{(0.0-1.0)} which indicates the proportion of NA values that will be tolerated to calculate Rao's index in each moving window over \emph{x}. If the relative proportion of NA's in a moving window is bigger than na.tolerance, then the value of the window will be set as NA, otherwise Rao's index will be calculated considering the non-NA values. Default values is 0.0 (i.e., no tolerance for NA's).}
\item{debugging}{a boolean variable set to FALSE by default. If TRUE, additional messages will be printed. For de-bugging only.}
}
\value{
Matrix or a list of matrixes with the Renyi index computed through moving window of the given size.
}
\author{
Marcantonio Matteo \email{marcantoniomatteo@gmail.com} \cr
Martina Iannacito \email{martina.iannacito@inria.fr} \cr
Duccio Rocchini \email{duccio.rocchini@unibo.it} \cr
}
\seealso{
\code{\link{Renyi}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}
