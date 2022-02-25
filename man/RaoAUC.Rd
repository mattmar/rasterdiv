\name{RaoAUC}
\alias{RaoAUC}
\title{Accumulation function for parametric Rao's index of quadratic entropy (Q)}
\description{
RaoAUC computes the accumulation function (integral or area under the curve) of the parametric version of Rao's index of quadratic entropy (Q) on different classes of numeric matrices using a moving window algorithm.
}
\usage{
RaoAUC(alphas=1:5, x, dist_m="euclidean", window=9,
  method="classic", rasterAUC=TRUE, lambda=0,
  na.tolerance=1.0, rescale=FALSE, diag=TRUE, 
  simplify=2, np=1, cluster.type="SOCK", debugging=FALSE)
  }
  \arguments{
  \item{alphas}{a continuous vector of alphas in the form start:end over which integrated the parametric Rao's index.  Default value is 1:5.}
  \item{x}{input data may be a matrix, a Spatial Grid Data Frame, a RasterLayer or a list of these objects. In the latter case, if \code{method="classic"} only the first element of the list will be considered.}
  \item{dist_m}{define the type of distance to be calculated between numerical categories. \code{dist_m} can be a character string which defines the name of the distance to derive such as "euclidean". The distance names allowed are the same as for \code{proxy::dist}. Alternatively, \code{dist_m} can be a function which calculates an user defined distance, (i.e., \code{function(x,y) {return(cos(y-x)-sin(y-x))}}) or a matrix of distances. If \code{method="multidimension"} then only "euclidean", "manhattan", "canberra", "minkowski" and "mahalanobis" can be used. Default value is "euclidean".}
  \item{window}{the side of the square moving window, it must be a odd numeric value greater than 1 to ensure that the target pixel is in the centre of the moving window. Default value is 3. If \code{proxy::dist} is a matrix then the function will assume that this is the distance matrix, and therefore no distance will be derived.}  
  \item{method}{Currently, there are two ways to calculate the parametric version of Rao's index. If \code{method="classic"}, then the normal parametric Rao's index will be calculated on a single matrix. If \code{method="multidimension"} (experimental!) a list of matrices must be provided as input. In the latter case, the overall distance matrix will be calculated in a multi- or hyper-dimensional system by using the distance measure defined through the function argument \code{dist_m}. Each pairwise distance is then multiplied by the inverse of the squared number of pixels in the considered moving window, and the Rao's Q is finally derived by applying a summation. Default value is \emph{"classic"}}
  \item{rasterAUC}{Boolean, if TRUE the output will be a RasterLayer object with \emph{x} as a raster template.}
  \item{lambda}{the value of the lambda of Minkowski's distance. Considered only if \code{dist_m = "minkowski"} and \code{method="multidimension". Default value is 0.}}
  \item{na.tolerance}{Numeric value \eqn{(0.0-1.0)} which indicates the proportion of NA values that will be tolerated to calculate parametric Rao's index in each moving window over \emph{x}. If the relative proportion of NA's in a moving window is bigger than na.tolerance, then the value of the window will be set as NA, otherwise Rao's index will be calculated considering the non-NA values. Default values is 1.0 (i.e., no tolerance for NA's). Default value is 1.0.}
  \item{rescale}{Boolean. Considered only if \code{method="multidimension"}. If TRUE, each element of \code{x} is rescaled and centred.}
  \item{diag}{Boolean. If TRUE then the diagonal of the distance matrix is filled with 0's, otherwise with NA's. If \code{diag=TRUE} and \code{alpha=0}, the output matrix will inexorably be composed of 0's.}
  \item{simplify}{Number of decimal places to be retained to calculate distances in Rao's index. Only if \emph{x} is floats.}
  \item{np}{the number of processes (cores) which will be spawned. Default value is 2.}
  \item{cluster.type}{the type of cluster which will be created. The options are \emph{"MPI"} (which calls "makeMPIcluster"), \emph{"FORK"} and \emph{"SOCK"} (which call "makeCluster"). Default type is \emph{"SOCK"}.}
  \item{debugging}{a boolean variable set to FALSE by default. If TRUE, additional messages will be printed. For debugging only.}
  }
  \details{
  The accumulation function for the parametric Rao's Index (\eqn{Q}) is calculated integrating numerically over a range of alphas. *RaoAUC* is therefore equal to \eqn{(\int_{a}^{b} {1\over{N^4}}\cdot{d_{i,j}^{\alpha}})^{1\over{\alpha}} dx}. Where \emph{N} is the number of pixels  in a moving window, and \emph{alpha} is a weight assigned to distances.}
  \value{
  A matrix of dimension \code{dim(x)}. If \code{rasterAUC=TRUE}, then the output is a RasterLayer with \emph{x} as template.}
  \references{
  [1] Rocchini, D., M. Marcantonio, and C. Ricotta (2017). Measuring Rao’s Q diversity index from remote sensing: An open source solution. Ecological Indicators. 72: 234–238.
  }
  \author{
  Matteo Marcantonio \email{marcantoniomatteo@gmail.com}
  }
  \seealso{
  \code{\link{paRao}}
  }
  \examples{
  #Minimal example; RaoAUC with alphas ranging from 1 to 10
  a <- matrix(c(10,10,10,20,20,20,20,30,30), ncol=3, nrow=3)
  out <- RaoAUC(alphas=1:10, x=a, window=3, dist_m="euclidean", na.tolerance=1, rasterAUC=TRUE)
  }