% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/paRao.R
\name{paRao}
\alias{paRao}
\title{Parametric Rao's index of quadratic entropy (Q)}
\usage{
paRao(
  x,
  area = NULL,
  field = NULL,
  dist_m = "euclidean",
  window = 9,
  alpha = 1,
  method = "classic",
  rasterOut = TRUE,
  lambda = 0,
  na.tolerance = 1,
  rescale = FALSE,
  diag = TRUE,
  simplify = 0,
  np = 1,
  cluster.type = "SOCK",
  progBar = TRUE,
  debugging = FALSE,
  time_vector = NA,
  stepness = -0.5,
  midpoint = 35,
  cycle_length = "year",
  time_scale = "day"
)
}
\arguments{
\item{x}{Input data may be a matrix, a Spatial Grid Data Frame, a SpatRaster, or a list of these objects.}

\item{area}{Input vector area layer for area-based calculation.}

\item{field}{Column name of the vector area layer to use to calculate the index.}

\item{dist_m}{Define the type of distance to be calculated between numerical categories. `dist_m` can be a character string which defines the name of the distance to derive such as "euclidean". The distance names allowed are the same as for \code{proxy::dist}. Alternatively, `dist_m` can be a function which calculates a user-defined distance, (i.e., \code{function(x,y) {return(cos(y-x)-sin(y-x))}}) or a matrix of distances. If `method="multidimension"` then only "euclidean", "manhattan", "canberra", "minkowski" and "mahalanobis" can be used. Default value is "euclidean". If `dist_m` is a matrix, then the function will assume that the matrix contains the distances. Moreover \emph{"twdtw"} (time weighted dynamic time warping) can be used as a way to calculate distances for time series in the `paRao` multidimensional mode.}

\item{window}{The side of the square moving window, it must be a vector of odd numeric values greater than 1 to ensure that the target pixel is in the centre of the moving window. Default value is 3. `window` can be a vector with length greater than 1, in this case, Rao's index will be calculated over `x` for each value in the vector.}

\item{alpha}{Weight for the distance matrix. If `alpha = 0`, distances will be averaged with a geometric average, if `alpha=1` with an arithmetic mean, if `alpha = 2` with a quadratic mean, `alpha = 3` with a cubic mean, and so on. if `alpha` tends to infinite (i.e., higher than the maximum integer allowed in R) or `alpha=Inf`, then the maximum distance will be taken. `alpha` can be a vector with length greater than 1, in this case, Rao's index will be calculated over `x` for each value in the vector.}

\item{method}{Currently, there are two ways to calculate the parametric version of Rao's index. If `method="classic"`, then the normal parametric Rao's index will be calculated on a single matrix. If `method="multidimension"` (experimental!), a list of matrices must be provided as input. In the latter case, the overall distance matrix will be calculated in a multi- or hyper-dimensional system by using the distance measure defined through the function argument `dist_m`. Each pairwise distance is then multiplied by the inverse of the squared number of pixels in the considered moving window, and the Rao's Q is finally derived by applying a summation. Default value is `"classic"`.}

\item{rasterOut}{Boolean, if TRUE the output will be a SpatRaster object with `x` as a template.}

\item{lambda}{The value of the lambda of Minkowski's distance. Considered only if `dist_m = "minkowski"` and `method="multidimension"`. Default value is 0.}

\item{na.tolerance}{Numeric value (0.0-1.0) which indicates the proportion of NA values that will be tolerated to calculate Rao's index in each moving window over `x`. If the relative proportion of NA's in a moving window is bigger than `na.tolerance`, then the value of the window will be set as NA, otherwise Rao's index will be calculated considering the non-NA values. Default value is 1.0.}

\item{rescale}{Boolean. Considered only if `method="multidimension"`. If TRUE, each element of `x` is rescaled and centred.}

\item{diag}{Boolean. If TRUE then the diagonal of the distance matrix is filled with 0's, otherwise with NA's. If `diag=TRUE` and `alpha=0`, the output matrix will inexorably be 0's.}

\item{simplify}{Number of decimal places to be retained to calculate distances in Rao's index. Default `simplify=0`.}

\item{np}{The number of processes (cores) which will be spawned. Default value is 2.}

\item{cluster.type}{The type of cluster which will be created. The options are `"MPI"` (which calls "makeMPIcluster"), `"FORK"`, and `"SOCK"` (which call "makeCluster"). Default type is `"SOCK"`.}

\item{progBar}{logical. If TRUE a progress bar is shown.}

\item{debugging}{A boolean variable set to FALSE by default. If TRUE, additional messages will be printed. For debugging only.}

\item{time_vector}{time;}

\item{stepness}{numeric; steepness of the logistic function.}

\item{midpoint}{numeric; midpoint of the logistic function}

\item{cycle_length}{string; The length of the cycle. Can be a numeric value or a string specifying the units ('year', 'month', 'day', 'hour', 'minute', 'second'). When numeric, the cycle length is in the same units as time_scale. When a string, it specifies the time unit of the cycle.}

\item{time_scale}{string; Specifies the time scale for the conversion. Must be one of 'year', 'month', 'day', 'hour', 'minute', 'second'. When cycle_length is a string, time_scale changes the unit in which the result is expressed. When cycle_length is numeric, time_scale is used to compute the elapsed time in seconds.}
}
\value{
A list of matrices of dimension `dim(x)` with length equal to the length of `alpha`. If `rasterOut=TRUE` and `x` is a SpatRaster, then the output is a list of SpatRaster objects.
}
\description{
It computes the parametric version of Rao's index of quadratic entropy (Q) on different classes of numeric matrices using a moving window algorithm.
}
\details{
The parametric Rao's Index (Q) is an extension of Rao's Index which considers a generalized mean between distances. The general formula for the parametric Rao's index is Q_a = \deqn{Q = \sum_{i, j} p_i p_j d_{ij}^{\alpha}}. Where `N` is the number of numerical categories, `i` and `j` are pair of numerical categories in the same moving window, and `alpha` is a weight given to distances. In the "multidimension" Rao's index, first the distances among categories are calculated considering more than one feature, and then the overall Rao's Q is derived by using these distances.
}
\examples{
\dontrun{
# loading data
data(volcano)
r <- terra::rast(volcano)

# we want to compute Rao's index on this data using a 3x3 window
res <- paRao(x = r, window = 3, alpha = 2, method = "classic")
terra::plot(res[[1]][[1]])
}

}
\references{
Rao, C. R. (1982). Diversity and dissimilarity coefficients: A unified approach. Theoretical Population Biology, 21(1), 24-43.
}
