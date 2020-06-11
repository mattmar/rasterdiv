\name{paRao}
\alias{paRao}
\title{Parametric Rao's index of quadratic entropy (Q)}
\description{
Computes the parametric version of Rao's index of quadratic entropy (Q) on different classes of numeric matrices using a moving window algorithm.
}
\usage{
paRao(x, dist_m="euclidean", window=9, alpha=1,
  method="classic", rasterOut=TRUE, lambda=0,
  na.tolerance=1.0, rescale=FALSE, diag=TRUE, 
  simplify=2, np=1,cluster.type="SOCK", debugging=FALSE)
  }
  \arguments{
  \item{x}{input data may be a matrix, a Spatial Grid Data Frame, a RasterLayer or a list of these objects. In the latter case, if \code{method="classic"} only the first element of the list will be considered.}
  \item{dist_m}{define the type of distance to be calculated between numerical categories. \code{dist_m} can be a character string which defines the name of the distance to derive such as "euclidean". The distance names allowed are the same as for \code{proxy::dist}. Alternatively, \code{dist_m} can be a function which calculates an user defined distance, (i.e., \code{function(x,y) {return(cos(y-x)-sin(y-x))}}) or a matrix of distances. If \code{method="multidimension"} then only "euclidean", "manhattan", "canberra", "minkowski" and "mahalanobis" can be used. Default value is "euclidean".}
  \item{window}{the side of the square moving window, it must be a odd numeric value greater than 1 to ensure that the target pixel is in the center of the moving window. Default value is 3. If \code{proxy::dist} is a matrix then the function will assume that this is the distance matrix, and therefore no distance will be derived.}
  \item{alpha}{weight for the distance matrix. If \code{alpha = 0}, distances will be averaged with an geometric mean, if \code{alpha=1} with an arithmetic mean, if \code{alpha = 2} with a quadratic mean, \code{alpha = 3} with a cubic mean and so on and so forth. if alpha tends to infinite (i.e., higher than the maximum integer allowed in R) then the maximum of the distances will be taken. \emph{"alpha"} can be a vector with length greater than 1, in this case Rao's index will be calculated over \emph{x} for each value in the sequence.}
  \item{method}{Currently, there are two ways to calculate the parameteric version of Rao's index. If \code{method="classic"}, then the normal parameteric Rao's index will be calculated on a single matrix. If \code{method="multidimension"} (experimental!) a list of matrices must be provided as input. In the latter case, the overall distance matrix will be calculated in a multi- or hyper-dimensional system by using the distance measure defined through the function argument \code{dist_m}. Each pairwise distance is then multiplied by the inverse of the squared number of pixels in the considered moving window, and the Rao's Q is finally derived by applying a summation. Default value is \emph{"classic"}}
  \item{rasterOut}{Boolean, if TRUE the output will be a RasterLayer object with \emph{x} as a raster template.}
  \item{lambda}{the value of the lambda of Minkowski's distance. Considered only if \code{dist_m = "minkowski"} and \code{method="multidimension". Default value is 0.}}
  \item{na.tolerance}{Numeric value \eqn{(0.0-1.0)} which indicates the proportion of NA values that will be tolerated to calculate Rao's index in each moving window over \emph{x}. If the relative proportion of NA's in a moving window is bigger than na.tolerance, then the value of the window will be set as NA, otherwise Rao's index will be calculated considering the non-NA values. Default values is 0.0 (i.e., no tolerance for NA's). Defaut value is 1.0.}
  \item{rescale}{Boolean. Considered only if \code{method="multidimension"}. If TRUE, each element of \code{x} is rescaled and centered.}
  \item{diag}{Boolean. If TRUE then the diagonal of the distance matrix is filled with 0's, otherwise with NA's. If \code{diag=TRUE} and \code{alpha=0}, the output matrix will inexorably be composed of 0's.}
  \item{simplify}{Number of decimal places to be retained to calculate distances in Rao's index. Only if \emph{x} is floats.}
  \item{np}{the number of processes (cores) which will be spawned. Default value is 2.}
  \item{cluster.type}{the type of cluster which will be created. The options are \emph{"MPI"} (which calls "makeMPIcluster"), \emph{"FORK"} and \emph{"SOCK"} (which call "makeCluster"). Default type is \emph{"SOCK"}.}
  \item{debugging}{a boolean variable set to FALSE by default. If TRUE, additional messages will be printed. For debugging only.}
  }
  \details{
  The parametric Rao's Index (\eqn{Q}) is an extenstion of Rao's Index which considers a generalised mean between distances. The generalised formula for the parametric Rao's index is \eqn{Q_a = (\sum_{i,j=1}^{N}{1/N^2}\times{d^{\alpha}_{j,j}})^{1/{\alpha}}}. Where \emph{N} is the number of numerical categories, \emph{i} and \emph{j} are pair of numerical categories in the same moving window and \emph{alpha} is a weight given to distances. In the "multidimension" Rao's index, first distances among categories are  calculated considering more than one layer, then the pairwise distance between each pair of numerical categories is multiplied to the square of the size of the moving window (this is somewhat the same as to calculate the variance of the multidimensional distance [1].). \cr
  The theoretical minimum of Rao's Q is 0, when all categories in a window have distance 0. If the distance chosen to calculate Rao's Index ranges between 0 and 1, the maximum value of Rao's Index equals the Simpson Index of Diversity \eqn{1-1/S_i} where \emph{S} is the number of categories in window \emph{i} (given alpha=1).}
  \value{
  A list of matrices of dimension \code{dim(x)} with length equal to the length of \emph{alpha}. If \code{rasterOut=TRUE}, then the output is a list of RasterLayer objects.
  }
  \references{
  [1] Rocchini, D., M. Marcantonio, and C. Ricotta (2017). Measuring Rao’s Q diversity index from remote sensing: An open source solution. Ecological Indicators. 72: 234–238.
  }
  \author{
  Matteo Marcantonio \email{marcantoniomatteo@gmail.com}\cr
  Duccio Rocchini \email{duccio.rocchini@unibo.it}
  }

  \examples{
  #Minimal example; compute classic Rao's index 
  a <- matrix(c(10,10,10,20,20,20,20,30,30),ncol=3,nrow=3)
  out <- paRao(x=a,window=3,dist_m="euclidean",na.tolerance=0,alpha=1)

  #Compute Rao's index with a vector of alpha values 
  out <- paRao(x=a,window=3,dist_m="euclidean",na.tolerance=0,alpha=1:5)

  #Compute multidimensional Rao's index with alpha = 1 
  a <- matrix(c(10,10,10,20,20,20,20,30,30),ncol=3,nrow=3)
  b <- matrix(c(100,100,100,200,200,200,200,300,300),ncol=3,nrow=3)

  out <- paRao(x=list(a,b),window=3,dist_m="euclidean",
    na.tolerance=0,alpha=1,method="multidimensional")
    }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}% use one of  RShowDoc("KEYWORDS")