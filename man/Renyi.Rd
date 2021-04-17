\name{Renyi}
\alias{Renyi}
\title{Renyi's entropy (\eqn{H})}
\description{
Computes Renyi's entropy (\eqn{{}^qH}) on different classes of numeric matrices using a moving window algorithm.
}
\usage{
Renyi(x, window=3, alpha=1, base=exp(1), rasterOut=TRUE, 
np=1.0, na.tolerance=, cluster.type="SOCK", debugging=FALSE)}
\arguments{
\item{x}{input data may be a matrix, a Spatial Grid Data Frame, a RasterLayer or a list of these objects. In the latter case, only the first element of the list will be considered.}
\item{window}{the side of the square moving window, it must be a odd numeric value greater than 1 to ensure that the target pixel is in the centre of the moving window. Default value is 3.}
\item{alpha}{Order of diversity to compute the index. If \code{alpha} is a vector with length greater than 1, then the index will be calculated over \code{x} for each value in the sequence.}
\item{base}{a numerical value which defines the base of the logarithm in Renyi's entropy formula. Default value is exp(1).}
\item{rasterOut}{Boolean, if TRUE output will be in RasterLayer format with \emph{x} as template.}
\item{np}{the number of processes (cores) which will be spawned. Default value is 1.}
\item{na.tolerance}{a numeric value \eqn{(0.0-1.0)} which indicates the proportion of NA values that will be tolerated to calculate Renyi's index in each moving window over \emph{x}. If the relative proportion of NA's in a moving window is bigger than na.tolerance, then the value of the window will be set as NA, otherwise Rao's index will be calculated considering the non-NA values. Default values is 1.0 (i.e., no tolerance for NA's).}
\item{cluster.type}{the type of cluster which will be created. The options are "MPI" (calls "makeMPIcluster"), "FORK" and "SOCK" (call "makeCluster"). Default type is "SOCK".}
\item{debugging}{a boolean variable set to FALSE by default. If TRUE, additional messages will be printed. For debugging only.}
}
\details{
Renyi's entropy (\eqn{{}^qH}) is calculated on a numerical matrix as \eqn{{}^qH = {1\over(1-q)} \ln(\sum_{i=1}^{R} {p^q}_i)}, where \emph{q} is the considered order of diversity (\code{alpha}), \emph{R} is the total number of categories (i.e., unique numerical values in the considered numerical matrix) and \emph{p} is the relative abundance of each category. If q=1, Shannon.R is called to calculate \eqn{H'} instead of the indefinite \eqn{{}^1D}, if \eqn{p > 2*10^9}, then BerkgerParker.R is called to calculate \eqn{log(1/{}^\infty H)}. Renyi's entropy of low order weight more rare numerical categories, whereas values of higher order weight more dominant categories.
}
\value{
A list of matrices with length equal to the length of "alpha". If length of "alpha" is 1, then a matrix of dimension \code{dim(x)}.
}
\references{
Rényi, A., 1970. Probability Theory. North Holland Publishing Company, Amsterdam.
}
\author{
Matteo Marcantonio \email{marcantoniomatteo@gmail.com} \cr
Martina Iannacito \email{martina.iannacito@inria.fr} \cr
Duccio Rocchini \email{duccio.rocchini@unibo.it} \cr
}
\note{
Linux users need to install libopenmpi for MPI parallel computing. Linux Ubuntu users may try:
apt-get update; apt-get upgrade; apt-get install mpi; apt-get install libopenmpi-dev; apt-get install r-cran-rmpi

Microsoft Windows users may need some additional work to use "MPI", see:\cr
https://bioinfomagician.wordpress.com/2013/11/18/installing-rmpi-mpi-for-r-on-mac-and-windows/
}

\seealso{
\code{\link{Shannon}},
\code{\link{BergerParker}}
}
\examples{
#Minimal example; compute Renyi's index with alpha 1:5 
a <- matrix(c(10,10,10,20,20,20,20,30,30),ncol=3,nrow=3)
renyi <- Renyi(x=a,window=3,alpha=1:5)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}% use one of  RShowDoc("KEYWORDS")