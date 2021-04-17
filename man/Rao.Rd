\name{Rao}
\alias{Rao}
\title{Rao's index of quadratic entropy (Q)}
\description{
Deprecated, use \code{paRao(..., alpha=1)}. Computes Rao's index of quadratic entropy (Q) on different classes of numeric matrices using a moving window algorithm.
}
\usage{
Rao(x, dist_m="euclidean", window=9, rasterOut = TRUE,
 mode="classic",lambda=0, shannon=FALSE, rescale=FALSE,
 na.tolerance=1.0, simplify=2, np=1, cluster.type="SOCK", 
 debugging=FALSE)
 }
 \arguments{
 \item{x}{input data may be a matrix, a Spatial Grid Data Frame, a RasterLayer or a list of these objects. In the latter case, if \code{mode="classic"} only the first element of the list will be considered.}
 \item{dist_m}{define the type of distance to be calculated between numerical categories. \code{dist_m} can be a character string which defines the name of the distance to derive such as "euclidean". The distance names allowed are the same as for \code{proxy::dist}. Alternatively, \code{dist_m} can be a function which calculates an user defined distance, (i.e., \code{function(x,y) {return(cos(y-x)-sin(y-x))}}) or a matrix of distances. If \code{mode="multidimension"} then only "euclidean", "manhattan", "canberra", "minkowski" and "mahalanobis" can be used. Default value is "euclidean".}
 \item{window}{the side of the square moving window, it must be a odd numeric value greater than 1 to ensure that the target pixel is in the centre of the moving window. Default value is 3. If \code{proxy::dist} is a matrix then the function will assume that this is the distance matrix, and therefore no distance will be derived.}
 \item{rasterOut}{Boolean, if TRUE the output will be in RasterLayer format with \emph{x} as template.}
 \item{mode}{currently, there are two modes to calculate Rao's index. If mode is "classic", then the classic Rao's index will be calculated on one single matrix. If mode is "multidimension" (experimental!) a list of matrices must be provided as input. In this latter case, the overall distance matrix will be calculated in a multi- or hyper-dimensional system by using the measure defined through the function argument \code{dist_m}. Each pairwise distance is then multiplied by the inverse of the squared number of pixels in the considered moving window, and the Rao's Q is finally derived by applying a summation. Default value is "classic".}
 \item{lambda}{the value of the lambda of Minkowski's distance. Considered only if \code{dist_m = "minkowski"} and \code{mode="multidimension". Default value is 0.}}
 \item{shannon}{a boolean variable set to FALSE by default. If TRUE, a matrix with Shannon index will be also calculated.  Default value is FALSE.}
 \item{rescale}{a boolean variable set to FALSE by default. If TRUE, \emph{x} will be scaled and centred to standardise different matrices if mode is "multidimension". Default value is FALSE.}
 \item{na.tolerance}{a numeric value \eqn{(0.0-1.0)} which indicates the proportion of NA values that will be tolerated to calculate Rao's index in each moving window over \emph{x}. If the relative proportion of NA's in a moving window is bigger than na.tolerance, then the value of the window will be set as NA, otherwise Rao's index will be calculated considering the non-NA values. Default values is 1.0 (i.e., full tolerance for NA's).}
 \item{simplify}{Number of decimal places to be retained to calculate distances in Rao's index. Only if \emph{x} is floats.}
 \item{np}{the number of processes (cores) which will be spawned. Default value is 1.}
 \item{cluster.type}{the type of cluster which will be created. The options are \code{"MPI"} (which calls "makeMPIcluster"), \code{"FORK"} and \code{"SOCK"} (which call "makeCluster"). Default type is \code{"SOCK"}.}
 \item{debugging}{a boolean variable set to FALSE by default. If TRUE, additional messages will be printed. For debugging only.}
 }
 \details{
 Classical Rao's Index (\eqn{Q}) is calculated on a numerical matrix as \eqn{Q = \sum_{i=1}^{R}{\sum_{j=1}^{R}{d_{i,j}\times{p_i}\times{p_j}}}} [1]. Where R is the number of categories, whereas \emph{i} and \emph{j} are pair of numerical categories in the same moving window. In the "multidimension" Rao's index, distances among categories are calculated considering more than one layer, then the pairwise distance between each pair of numerical categories is multiplied to the square of the size of the moving window (which is somewhat the same as to calculate the variance of the multidimensional distance [2].). \cr
 The theoretical minimum of Rao's Q is 0, when all categories in a window have distance 0. If the distance chosen to calculate Rao's Index ranges between 0 and 1, the maximum value of Rao's Index equals the Simpson Index of Diversity \eqn{1-1/S_i} where \emph{S} is the number of categories in window \emph{i}.}
 \value{
 If shannon=TRUE, a list of matrices of length two, otherwise a matrix of dimension \code{dim(x)}.
 }
 \references{
 [1] Rao, C.R. (1982). Diversity and dissimilarity coefficients: a unified approach. Theoretical Population Biology, 21: 2443.
 [2] Rocchini, D., M. Marcantonio, and C. Ricotta (2017). Measuring Rao’s Q diversity index from remote sensing: An open source solution. Ecological Indicators. 72: 234–238.
 }
 \author{
 Matteo Marcantonio \email{marcantoniomatteo@gmail.com}\cr
 Duccio Rocchini \email{duccio.rocchini@unibo.it}
 }

 \examples{
# Minimal example; compute Rao's index 
 \dontrun{
  a <- matrix(c(10,10,10,20,20,20,20,30,30),ncol=3,nrow=3)
  out <- Rao(x=a,window=3,dist_m="euclidean",na.tolerance=1.0,shannon=FALSE,mode="classic")

# Compute both Rao and Shannon index 
  out <- Rao(x=a,window=3,dist_m="euclidean",na.tolerance=1.0,shannon=TRUE,mode="classic")

# Compute both Rao and Shannon index multiple windows
  out <- Rao(x=a,window=c(3,5),dist_m="euclidean",na.tolerance=1.0,shannon=TRUE,mode="classic")

# Compute multidimension Rao's index rescaling the input matrices
  a <- matrix(c(10,10,10,20,20,20,20,30,30),ncol=3,nrow=3)
  b <- matrix(c(0.5,0.5,0.1,0.1,0.3,0.3,0.3,0.3,0.3),ncol=3,nrow=3)
  out <- Rao(x=list(a,b),window=3,dist_m="euclidean",na.tolerance=1.0,
    mode="multidimension",rescale=TRUE,debugging=TRUE)
}
}