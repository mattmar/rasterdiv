\name{mpaRaoAreaS}
\alias{mpaRaoAreaS}
\title{
Area-Based Sequential Parametric Rao's index of quadratic entropy (Q)
}
\description{
Area-Based Sequential Parametric Rao's index of quadratic entropy (Q)
}
\usage{
mpaRaoAreaS(rasterm, area, alpha, simplify, dist_m, rescale, lambda, window)
}
\arguments{
\item{rasterm}{input raster data.}
\item{area}{input vector data.}
\item{alpha}{alpha value for order of diversity in Hill's Index.}
\item{simplify}{Rounding parameter.}
\item{dist_m}{Type of distance used.}
\item{rescale}{Scale and centre values in each of the element of x.}
\item{lambda}{Lambda value for Minkowski distance.}
\item{window}{moving window.}

}
\value{
A vector as the input with Rao's index pasted in additional columns.
}
\author{
Marcantonio Matteo \email{marcantoniomatteo@gmail.com} \cr
Duccio Rocchini \email{duccio.rocchini@unibo.it} \cr
Michele Torresani \email{michele.torresani@unibo.it} \cr
}
\seealso{
\code{\link{paRao}}
}
\keyword{internal}
