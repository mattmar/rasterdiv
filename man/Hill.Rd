\name{Hill}
\alias{Hill}
\title{Hill's index of diversity - Hill numbers (D)}
\description{
Computes Hill's index of diversity (Hill numbers) on different classes of numeric matrices using a moving window algorithm.
}
\usage{
Hill(x, window = 3, alpha = 1, rasterOut=TRUE, 
	np = 1, na.tolerance=1, cluster.type = "SOCK", 
	debugging = FALSE)
}
\arguments{
\item{x}{input data may be a matrix, a Spatial Grid Data Frame, a RasterLayer or a list of these objects. In the latter case, only the first element of the list will be considered.}
\item{window}{the side of the square moving window, it must be a odd numeric value greater than 1 to ensure that the target pixel is in the center of the moving window. Default value is 3.}
\item{alpha}{Order of the Hill number to compute the index. If \code{"alpha"} is a vector with length greater than 1, then the index will be calculated over \code{x} for each value in the sequence.}
\item{rasterOut}{Boolean, if TRUE output will be in RasterLayer format with \emph{x} as template.}
\item{np}{the number of processes (cores) which will be spawned. Default value is 1.}
\item{na.tolerance}{a numeric value \eqn{(0.0-1.0)} which indicates the proportion of NA values that will be tolerated to calculate Rao's index in each moving window over \emph{x}. If the relative proportion of NA's in a moving window is bigger than na.tolerance, then the value of the window will be set as NA, otherwise Rao's index will be calculated considering the non-NA values. Default values is 0.0 (i.e., no tolerance for NA's).}
\item{cluster.type}{the type of cluster which will be created. The options are "MPI" (calls "makeMPIcluster"), "FORK" and "SOCK" (call "makeCluster"). Default type is "SOCK".}
\item{debugging}{a boolean variable set to FALSE by default. If TRUE, additional messages will be printed. For debugging only.}
}
\details{
Hill numbers (\eqn{{}^qD}) are calculated on a numerical matrices as \eqn{{}^qD = (\sum_{i=1}^{R} {p^q}_i)^{1/(1-q)}},where \emph{q} is the order of the Hill number, \emph{R} is the total number of categories (i.e., unique numerical values in a numerical matrix), \emph{p} is the relative abundance of each category. When q=1, Shannon.R is called to calculate \eqn{exp(H^1)} instead of the indefinite \eqn{{}^1D}. if \eqn{q > 2*10^9}, BerkgerParker.R is called to calculate \eqn{1/{{}^\infty D}}. Hill numbers of low order weight more rare categories, whereas Hill numbers of higher order weight more dominant categories.
}
\value{
A list of matrices of dimension \code{dim(x)} with length equal to the length of \code{alpha}.
}
\references{
Hill, M.O. (1973). Diversity and evenness: a unifying notation and its consequences. Ecology 54, 427-431.
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
%% ~~objects to See Also as , ~~~
\code{\link{BergerParker}}
\code{\link{Shannon}}
}
\examples{
#Minimal example; compute Hill's index with alpha 1:5 
a <- matrix(c(10,10,10,20,20,20,20,30,30),ncol=3,nrow=3)
hill <- Hill(x=a,window=3,alpha=1:5)
}
\keyword{methods}