\name{Shannon}
\alias{Shannon}
\title{
Shannon's diversity index (H')
}
\description{
Computes Shannon's diversity index (H') on different classes of numeric matrices using a moving window algorithm.
}
\usage{
Shannon(x, window=3, rasterOut=TRUE, np=1, na.tolerance=1, 
	cluster.type="SOCK", debugging=FALSE)
}
\arguments{
\item{x}{input data may be a matrix, a Spatial Grid Data Frame, a RasterLayer or a list of these objects. In the latter case, only the first element of the list will be considered.}
\item{window}{the side of the square moving window, it must be a odd numeric value greater than 1 to ensure that the target pixel is in the center of the moving window. Default value is 3.}
 \item{rasterOut}{Boolean, if TRUE the output will be in RasterLayer format with \emph{x} as template.}
\item{np}{the number of processes (cores) which will be spawned. Default value is 1.}
\item{na.tolerance}{a numeric value \eqn{(0.0-1.0)} which indicates the proportion of NA values that will be tolerated to calculate Rao's index in each moving window over \emph{x}. If the relative proportion of NA's in a moving window is bigger than na.tolerance, then the value of the window will be set as NA, otherwise Rao's index will be calculated considering the non-NA values. Default values is 0.0 (i.e., no tolerance for NA's).}
\item{cluster.type}{the type of cluster which will be created. The options are \code{"MPI"} (calls "makeMPIcluster"), "FORK" and "SOCK" (call "makeCluster"). Default type is \code{"SOCK"}.}
\item{debugging}{a boolean variable set to FALSE by default. If TRUE, additional messages will be printed. For debugging only.}
}
\details{
Shannon's index (\eqn{H'}) is calculated on a numerical matrix as \eqn{H' = {\sum_{i=1}^{R} p_i \times log(p_i)}}, where \emph{R} is the total number of categories (i.e., unique numerical values in the considered numerical matrix) and \emph{p} is the relative abundance of each category.}
\value{
A numerical matrix with dimension \code{{dim(x)}}.
}
\references{
Shannon, C.E. (1948). A mathematical theory of communication. Bell System Technical Journal, 27: 379-423, 623-656.
}
\author{
Marcantonio Matteo \email{marcantoniomatteo@gmail.com} \cr
Martina Iannacito \email{martina.iannacito@inria.fr} \cr
Duccio Rocchini \email{duccio.rocchini@unibo.it} \cr}
\note{
Linux users need to install libopenmpi for MPI parallel computing. Linux Ubuntu users may try:
apt-get update; apt-get upgrade; apt-get install mpi; apt-get install libopenmpi-dev; apt-get install r-cran-rmpi

Microsoft Windows users may need some additional work to use "MPI", see:\cr
https://bioinfomagician.wordpress.com/2013/11/18/installing-rmpi-mpi-for-r-on-mac-and-windows/
}
\examples{
#Minimal example; compute Shannon's index  
a <- matrix(c(10,10,10,20,20,20,20,30,30),ncol=3,nrow=3)
shannon <- Shannon(x=a,window=3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}% use one of  RShowDoc("KEYWORDS")