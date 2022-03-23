\name{paRaoAreaS}
\alias{paRaoAreaS}
\title{
Area-Based Sequential Parametric Rao's index of quadratic entropy (Q)
}
\description{
Area-Based Sequential Parametric Rao's index of quadratic entropy (Q)
}
\usage{
paRaoAreaS(rasterm, area, alpha, simplify)
}
\arguments{
\item{rasterm}{input raster data.}
\item{area}{input vector data.}
\item{alpha}{alpha value for order of diversity in Hill's Index.}
\item{simplify}{Rounding parameter.}
}
\value{
A vector as the input with Rao's index pasted additional columns.
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
