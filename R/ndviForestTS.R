#' Simulated NDVI dataset
#'
#' A \code{list} of 8-bit matrices.
#'
#' This list represents a time series of NDVI values of a patch of forest over 3 years. 
#' It is stored as a \code{list}, suitable for explaining how to make helical plots.
#'
#' @format A \code{list} containing matrices:
#' \describe{
#'   \item{ndviForestTS}{List of matrixes of 9 cells simulating NDVI of a patch of forests over 3 years. Each matrix represents a day in the time series.}
#' }
#' @name ndviForestTS
#' @docType data
#' @keywords datasets
#' @examples
#' ndviForestTS <- readRDS(system.file("extdata", "ndviForestTS.rds", package = "rasterdiv"))
NULL