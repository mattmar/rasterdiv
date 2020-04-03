\name{Pielou}
\alias{Pielou}
\title{
Pielou's evenness index (E')
}
\description{
Computes Pielou's evenness index on different classes of numeric matrices using a moving window algorithm.}
\usage{
Pielou(x, window=3, np=1, na.tolerance=1, cluster.type="SOCK",
 debugging=FALSE)
}
\arguments{
\item{x}{input data may be a matrix, a Spatial Grid Data Frame, a RasterLayer or a list of these objects. In the latter case, only the first element of the list will be considered.}
\item{window}{the side of the square moving window, it must be a odd numeric value greater than 1 to ensure that the target pixel is in the center of the moving window. Default value is 3.}
\item{np}{the type of cluster which will be created. The options are \code{"MPI"} (calls "makeMPIcluster"), "FORK" and "SOCK" (call "makeCluster"). Default type is "SOCK".}
\item{na.tolerance}{a numeric value \eqn{(0.0-1.0)} which indicates the proportion of NA values that will be tolerated to calculate Rao's index in each moving window over \emph{x}. If the relative proportion of NA's in a moving window is bigger than na.tolerance, then the value of the window will be set as NA, otherwise Rao's index will be calculated considering the non-NA values. Default values is 0.0 (i.e., no tolerance for NA's).}
\item{cluster.type}{the type of cluster which will be created. The options are "MPI" (calls "makeMPIcluster"), "FORK" and "SOCK" (call "makeCluster"). Default type is \code{"SOCK"}.}
\item{debugging}{a boolean variable set to FALSE by default. If TRUE, additional messages will be printed. For debugging only.}
}
\details{
Pielou evenness's index (\eqn{E'}) is calculated on a numerical matrix as \eqn{E' = {{\sum_{i=1}^{R} p_i \times log(p_i)}\over{log(R)}}}, where \emph{R} is the total number of categories (i.e., unique numerical values in the considered numerical matrix) and \emph{p} is the relative abundance of each category. Pielou's evenness represents the ratio between the observed value of Shannon's Index and the value of Shannon's Index if all categories (R) had the same relative abundance.}
\value{
A numerical matrix with dimension as \code{dim(x)}.
}
\references{
Pielou, E.C. (1966). The measurement of diversity in different types of biological collections. Journal of Theoretical Biology, 13: 131-144.
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{Shannon}}
}

\examples{
#Minimal example; compute Shannon's index  
a <- matrix(c(10,10,10,20,20,20,20,30,30),ncol=3,nrow=3)
renyi <- Pielou(x=a,window=3)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}% use one of  RShowDoc("KEYWORDS")
