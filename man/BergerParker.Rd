\name{BergerParker}
\alias{BergerParker}
\title{
Berger-Parker's diversity index
}
\description{
Computes Berger-Parker's diversity index on different classes of numeric matrices using a moving window algorithm.
}
\usage{
BergerParker(x, window=3, rasterOut=TRUE, np=1, 
	na.tolerance=1.0, cluster.type="SOCK", debugging=FALSE)
}
\arguments{
\item{x}{input data may be a matrix, a Spatial Grid Data Frame, a RasterLayer or a list of these objects. In the latter case, only the first element of the list will be considered.}
\item{window}{the side of the square moving window, it must be a odd numeric value greater than 1 to ensure that the target pixel is in the centre of the moving window. Default value is 3.}
\item{rasterOut}{Boolean, if TRUE output will be in RasterLayer format with \emph{x} as template.}
\item{np}{the number of processes (cores) which will be spawned. Default value is 1.}
\item{na.tolerance}{a numeric value \eqn{(0.0-1.0)} which indicates the proportion of NA values that will be tolerated to calculate Berger-Parker's index in each moving window over \emph{x}. If the relative proportion of NA's in a moving window is bigger than na.tolerance, then the value of the window will be set as NA, otherwise Rao's index will be calculated considering the non-NA values. Default values is 1.0 (i.e., no tolerance for NA's).}
\item{cluster.type}{the type of cluster which will be created. The options are \code{"MPI"} (calls "makeMPIcluster"), \code{"FORK"} and \code{"SOCK"} (call "makeCluster"). Default type is \code{"SOCK"}.}
\item{debugging}{a boolean variable set to FALSE by default. If TRUE, additional messages will be printed. For de-bugging only.}
}
\details{
Berger-Parker's index is the relative abundance of the most abundant category (i.e., unique numerical values in the considered numerical matrix). Berger-Parker's index equals the logarithm of the inverse Renyi's index of order infinity, \eqn{log(1/{}^\infty H)} or the inverse of Hill's index of order infinity, \eqn{1/{}^\infty D}.
}
\value{
A numerical matrix with dimension as \code{dim(x)}.
}
\references{
Berger, W.H., Parker, F.L. (1970). Diversity of planktonic foraminifera in deep-sea sediments". Science, 168: 1345-1347.
}
\author{
Marcantonio Matteo \email{marcantoniomatteo@gmail.com} \cr
Martina Iannacito \email{martina.iannacito@inria.fr} \cr
Duccio Rocchini \email{duccio.rocchini@unibo.it} \cr
}
\note{
Linux users need to install libopenmpi for MPI parallel computing. Linux Ubuntu users may try:
apt-get update; apt-get upgrade; apt-get install mpi; apt-get install libopenmpi-dev; apt-get install r-cran-rmpi

Microsoft Windows users may need some additional work to use "MPI", see:\cr
https://bioinfomagician.wordpress.com/2013/11/18/installing-rmpi-mpi-for-r-on-mac-and-windows/
}

\examples{
#Minimal example; compute Renyi's index with alpha 1:5 
a <- matrix(c(10,10,10,20,20,20,20,30,30),ncol=3,nrow=3)
berpar <- BergerParker(x=a,window=3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}% use one of  RShowDoc("KEYWORDS")